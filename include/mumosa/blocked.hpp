#ifndef MUMOSA_BLOCKED_HPP
#define MUMOSA_BLOCKED_HPP

#include <vector>
#include <algorithm>
#include <cstddef>
#include <optional>

#include "knncolle/knncolle.hpp"
#include "sanisizer/sanisizer.hpp"
#include "scran_blocks/scran_blocks.hpp"
#include "quickstats/quickstats.hpp"

#include "simple.hpp"
#include "utils.hpp"

/**
 * @file blocked.hpp
 * @brief Compute distances to nearest neighbors with blocking.
 */

namespace mumosa {

/**
 * @brief Options for `compute_distance_blocked()`.
 */
struct BlockedOptions {
    /**
     * Number of neighbors for the nearest neighbor search.
     * Larger values improve stability at the risk of including biological heterogeneity into the distance.
     * `num_neighbors + 1` can also be interpreted as the expected minimum size of each subpopulation.
     */
    int num_neighbors = 20;

    /**
     * Policy to use for weighting the contribution from each block when computing the average distance.
     */
    scran_blocks::WeightPolicy block_weight_policy = scran_blocks::WeightPolicy::VARIABLE;

    /**
     * Parameters for the variable block weights, including the threshold at which blocks are considered to be large enough to have equal weight.
     * Only relevant when `BlockedOptions::block_weight_policy = scran_blocks::WeightPolicy::VARIABLE`.
     */
    scran_blocks::VariableWeightParameters variable_block_weight_parameters; 

    /**
     * Number of threads to use.
     * The parallelization mechanism is determined by `knncolle::parallelize()`.
     t*/
    int num_threads = 1;
};

/*
 * We don't apply block-specific scaling factors as we don't want to alter the relative values within the same modality.
 * We shouldn't have to do it in the first place - as it's the same modality! - but more importantly, we could introduce spurious differences between blocks.
 * In the simplest case, two blocks have the same subpopulation structure but the number of cells is different.
 * We would get different distances in each block due to density, causing us to scale each block differently.
 * More generally, we could expect differences in subpopulation structure between blocks, leading to different distances even in the absence of any batch effects.
 * (Mind you, differences in subpopulation structure also interfere with accurate scaling between modalities,
 * but any errors in scaling modalities are much less obvious than those from scaling blocks.)
 */

/**
 * Systematic differences between blocks can artificially inflate the distances to the nearest neighbors within a modality's embedding.
 * Specifically, strong batch effects can reduce the density of the local neighborhood by shifting cells elsewhere.
 * This increases the distance to the nearest neighbors compared to a modality without any batch effects,
 * even if the variance in the local neighborhood is the same between modalities.
 *
 * If the magnitude of the batch effects differ between modalities, this may introduce spurious differences in the median distance-to-neighbor.
 * To improve accuracy in the presence of blocks, this function calls `compute_distance()` on each entry of `prebuilts` separately.
 * It then computes a weighted average of the median distance and RMSDs across blocks (see `scran_blocks::compute_weights()` for details).
 * This ensures that arbitrary shifts in location between blocks have no effect on the distances to the nearest neighbors for each modality.
 *
 * @tparam Index_ Integer type of the number of cells.
 * @tparam Distance_ Floating-point type of the distances.
 *
 * @param prebuilts Vector of length equal to the number of blocks.
 * Each entry contains (1) the number of cells in the block and (2) a pointer to an array of length equal to the number of cells in this block.
 * The latter contains the distance of each cell to its \f$k\f$-nearest neighbor within that block.
 * @param options Further options.
 * 
 * @return Pair containing the weighted average of the median distance to the nearest neighbor (first)
 * and the root-mean-squared distance (second) across blocks.
 * These values can be used in `compute_scale()`.
 * If there are no non-empty blocks, both the median and RMSD are set to zero.
 */
template<typename Index_, typename Distance_>
std::pair<Distance_, Distance_> compute_distance_blocked(const std::vector<std::pair<Index_, Distance_*> >& blocks, const BlockedOptions& options) {
    const auto nblocks = blocks.size();
    auto block_weights = sanisizer::create<std::vector<Distance_> >(nblocks);
    scran_blocks::compute_weights(
        sanisizer::cast<std::size_t>(nblocks),
        [&](std::size_t b) -> Index_ { return blocks[b].first; },
        options.block_weight_policy,
        options.variable_block_weight_parameters,
        [&](std::size_t b, Distance_ w) -> void { block_weights[b] = w; }
    );

    const auto total_weight = [&]{
        quickstats::PairwiseSumWorkspace<Distance_> pswrk;
        quickstats::PairwiseSumOptions psopt;
        return quickstats::pairwise_sum(block_weights.size(), block_weights.data(), pswrk, psopt);
    }();

    auto outputs = sanisizer::create<std::vector<std::pair<Distance_, Distance_> > >(nblocks);
    knncolle::parallelize(options.num_threads, nblocks, [&](const int, I<decltype(nblocks)> start, I<decltype(nblocks)> length) -> void {
        for (I<decltype(nblocks)> b = start, bend = start + length; b < bend; ++b) {
            const auto curweight = block_weights[b];
            const auto curdist = compute_distance(blocks[b].first, blocks[b].second);
            outputs[b].first = curdist.first * curweight;
            outputs[b].second = curdist.second * curweight;
        }
    });

    std::pair<Distance_, Distance_> output{};
    for (I<decltype(nblocks)> b = 0; b < nblocks; ++b) {
        output.first += outputs[b].first;
        output.second += outputs[b].second;
    }

    if (total_weight) {
        output.first /= total_weight;
        output.second /= total_weight;
    }

    return output;
}

/**
 * Overload of `compute_distance_blocked()` that accepts a set of prebuilt neighbor search indices.
 *
 * @tparam Index_ Integer type of the number of cells.
 * @tparam Input_ Numeric type of the input data used to build the search index.
 * This is only required to define the `knncolle::Prebuilt` class and is otherwise ignored.
 * @tparam Distance_ Floating-point type of the distances.
 *
 * @param prebuilts Vector of length equal to the number of blocks.
 * Each entry contains a prebuilt neighbor search index for a single block.
 * A block with no observations may be represented by a null pointer.
 * @param workspace Workspace object, constructed with block sizes that match the number of observations in each entry of `prebuilts`.
 * This can be re-used across multiple `compute_distance_blocked()` calls with the same block sizes.
 * @param options Further options.
 * 
 * @return Pair containing the weighted average of the median distance to the nearest neighbor (first)
 * and the root-mean-squared distance (second) across blocks.
 * These values can be used in `compute_scale()`.
 */
template<typename Index_, typename Input_, typename Distance_>
std::pair<Distance_, Distance_> compute_distance_blocked(
    const std::vector<std::shared_ptr<const knncolle::Prebuilt<Index_, Input_, Distance_> > >& prebuilts,
    Distance_* const buffer,
    const BlockedOptions& options
) {
    const auto nblocks = prebuilts.size();
    std::size_t accumulated = 0;
    std::vector<std::pair<Index_, Distance_*> > blocks;
    blocks.reserve(nblocks);

    for (I<decltype(nblocks)> b = 0; b < nblocks; ++b) {
        const auto nobs = prebuilts[b]->num_observations();
        const auto capped_k = knncolle::cap_k(options.num_neighbors, nobs);

        knncolle::parallelize(options.num_threads, nobs, [&](const int, const Index_ start, const Index_ length) -> void {
            const auto searcher = prebuilts[b]->initialize();
            std::vector<Distance_> cur_distances;
            for (Index_ i = start, end = start + length; i < end; ++i) {
                searcher->search(i, capped_k, NULL, &cur_distances);
                if (cur_distances.size()) {
                    buffer[accumulated + i] = cur_distances.back();
                } else {
                    buffer[accumulated + i] = 0; // i.e., only distance is that to itself.
                }
            }
        });

        blocks.emplace_back(nobs, buffer + accumulated);
        accumulated += nobs;
    }

    return compute_distance_blocked(blocks, options);
}

/**
 * Overload of `compute_distance_blocked()` that accepts an embedding matrix with a block factor.
 *
 * @tparam Index_ Integer type of the number of cells.
 * @tparam Input_ Numeric type of the input data. 
 * @tparam Distance_ Floating-point type of the distances.
 * @tparam Matrix_ Class of the input data matrix for the neighbor search.
 * This should satisfy the `knncolle::Matrix` interface.
 *
 * @param num_dim Number of dimensions in the embedding.
 * @param num_cells Number of cells.
 * @param[in] data Pointer to an array containing the embedding matrix for a modality.
 * This should be stored in column-major layout where each row is a dimension and each column is a cell.
 * The number of rows and columns should be equal to `num_dim` and `num_cells`, respectively.
 * @param[in] blocks Pointer to an array of length equal to `num_cells`, containing the block assignment for each column of `data`.
 * Each value should be a non-negative integer in `[0, num_blocks)`.
 * @param num_blocks Number of blocks.
 * @param builder Algorithm to use for the neighbor search.
 * @param options Further options.
 * 
 * @return Pair containing the weighted average of the median distance to the nearest neighbor (first)
 * and the root-mean-squared distance (second) across blocks.
 * These values can be used in `compute_scale()`.
 */
template<typename Index_, typename Input_, typename Block_, typename Distance_, class Matrix_ = knncolle::Matrix<Index_, Input_> >
std::pair<Distance_, Distance_> compute_distance_blocked(
    const std::size_t num_dim,
    const Index_ num_cells,
    const Input_* const data,
    const Block_* const blocks,
    const std::size_t num_blocks,
    const knncolle::Builder<Index_, Input_, Distance_, Matrix_>& builder,
    Distance_* const buffer,
    const BlockedOptions& options
) {
    // Avoiding allocation of a temporary buffer if we're already dealing with contiguous blocks.
    auto block_details = sanisizer::create<std::vector<std::pair<Index_, Index_> > >(num_blocks);
    Index_ non_contiguous = 0;
    for (Index_ c = 0; c < num_cells; ++c) {
        auto& curblock = block_details[blocks[c]];
        if (curblock.second == 0) {
            curblock.first = c;
            curblock.second = 1;
        } else {
            non_contiguous += (c != curblock.first + curblock.second);
            ++curblock.second;
        }
    }

    const Input_* dataptr = data;
    std::optional<std::vector<Input_> > tmp_data;
    if (non_contiguous) {
        // Otherwise, we reorganize the data so that observations from the same batch are in a single block.
        Index_ accumulated = 0;
        auto offsets = sanisizer::create<std::vector<Index_> >(num_blocks);
        for (std::size_t b = 0; b < num_blocks; ++b) {
            offsets[b] = accumulated;
            block_details[b].first = accumulated;
            accumulated += block_details[b].second; // this won't overflow as we already know that num_cells fits in an Index_.
        }

        tmp_data.emplace(sanisizer::product<typename std::vector<Input_>::size_type>(num_dim, num_cells));
        for (Index_ c = 0; c < num_cells; ++c) {
            auto& off = offsets[blocks[c]];
            std::copy_n(
                data + sanisizer::product_unsafe<std::size_t>(c, num_dim),
                num_dim,
                tmp_data->data() + sanisizer::product_unsafe<std::size_t>(off, num_dim)
            );
            ++off;
        }

        dataptr = tmp_data->data();
    }

    auto prebuilts = sanisizer::create<std::vector<std::shared_ptr<const knncolle::Prebuilt<Index_, Input_, Distance_> > > >(num_blocks);
    knncolle::parallelize(options.num_threads, num_blocks, [&](const int, const std::size_t start, const std::size_t length) -> void {
        for (std::size_t b = start, end = start + length; b < end; ++b) {
            const auto sofar = block_details[b].first;
            const auto cursize = block_details[b].second;
            prebuilts[b] = builder.build_shared(knncolle::SimpleMatrix(num_dim, cursize, dataptr + sanisizer::product_unsafe<std::size_t>(sofar, num_dim)));
        }
    });

    return compute_distance_blocked(prebuilts, buffer, options);
}

}

#endif
