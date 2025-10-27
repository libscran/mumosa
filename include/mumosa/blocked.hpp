#ifndef MUMOSA_BLOCKED_HPP
#define MUMOSA_BLOCKED_HPP

#include <vector>
#include <algorithm>
#include <cstddef>

#include "knncolle/knncolle.hpp"
#include "sanisizer/sanisizer.hpp"
#include "scran_blocks/scran_blocks.hpp"

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

/**
 * @brief Workspace for `compute_distance_blocked()`.
 *
 * @tparam Distance_ Floating-point type of the distances.
 *
 * Instances of this class should typically be created by `prepare_workspace()`.
 */
template<typename Distance_>
struct BlockedWorkspace {
    /**
     * @cond
     */
    std::vector<Distance_> weights;
    Distance_ total_weight;

    std::vector<Distance_> distance_buffer;
    /**
     * @endcond
     */
};

/**
 * @tparam Index_ Integer type of the number of cells.
 * @tparam Distance_ Floating-point type of the distances.
 *
 * @param block_sizes Vector of length equal to the number of blocks, containing the number of observations in each block.
 * @param options Further options.
 * 
 * @return A workspace that can be re-used across multiple `compute_distance_blocked()` calls. 
 */
template<typename Distance_, typename Index_>
BlockedWorkspace<Distance_> create_workspace(const std::vector<Index_>& block_sizes, const BlockedOptions& options) {
    BlockedWorkspace<Distance_> output;
    output.weights = scran_blocks::compute_weights<Distance_>(block_sizes, options.block_weight_policy, options.variable_block_weight_parameters);
    output.total_weight = std::accumulate(output.weights.begin(), output.weights.end(), static_cast<Distance_>(0));

    Index_ max_size = 0;
    if (block_sizes.size()) {
        max_size = *std::max_element(block_sizes.begin(), block_sizes.end());
    }
    sanisizer::resize(output.distance_buffer, max_size);

    return output;
}

/**
 * NOTES:
 *
 * The local neighborhood variance can be considered as the variance within a particular region of the high-dimensional space.
 * The expectation of this variance should not be affected by the number of cells, but the distance to the neighbors will be affected if the density of cells changes.
 *
 * We do not apply block-specific scaling factors as we don't want to alter the relative values within the same modality.
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
 * @tparam Input_ Numeric type of the input data used to build the search index.
 * This is only required to define the `knncolle::Prebuilt` class and is otherwise ignored.
 * @tparam Distance_ Floating-point type of the distances.
 *
 * @param prebuilts Vector of length equal to the number of blocks.
 * Each entry contains a prebuilt neighbor search index for a single block.
 * A block with no observations may be represented by a null pointer.
 * @param workspace Workspace object, constructed with `block_sizes` that contain the number of observations in each entry of `prebuilts`.
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
    BlockedWorkspace<Distance_>& workspace,
    const BlockedOptions& options
) {
    Options simple_opt;
    simple_opt.num_neighbors = options.num_neighbors;
    simple_opt.num_threads = options.num_threads;

    std::pair<Distance_, Distance_> output(0, 0);

    const auto nblocks = prebuilts.size();
    for (I<decltype(nblocks)> b = 0; b < nblocks; ++b) {
        const auto curweight = workspace.weights[b];
        const auto& pbptr = prebuilts[b];
        if (curweight && pbptr && pbptr->num_observations()) {
            const auto curdist = compute_distance(*pbptr, workspace.distance_buffer.data(), simple_opt);
            output.first += curdist.first * curweight;
            output.second += curdist.second * curweight;
        }
    }

    if (workspace.total_weight) {
        output.first /= workspace.total_weight;
        output.second /= workspace.total_weight;
    }

    return output;
}

/**
 * Overload of `compute_distance()` that accepts an embedding matrix with contiguous blocks.
 *
 * @tparam Index_ Integer type of the number of cells.
 * @tparam Input_ Numeric type of the input data. 
 * @tparam Distance_ Floating-point type of the distances.
 * @tparam Matrix_ Class of the input data matrix for the neighbor search.
 * This should satisfy the `knncolle::Matrix` interface.
 *
 * @param num_dim Number of dimensions in the embedding.
 * @param block_sizes Number of cells in each block.
 * @param[in] data Pointer to an array containing the embedding matrix for a modality.
 * This should be stored in column-major layout where each row is a dimension and each column is a cell.
 * The number of rows should be equal to `num_dim` and the number of columns should be equal to the sum of `block_sizes`.
 * Cells from the first block should be stored in the first `block_sizes[0]` columns,
 * cells from the second block should be stored in the next `block_sizes[1]` columns,
 * and so on.
 * @param builder Algorithm to use for the neighbor search.
 * @param options Further options.
 * 
 * @return Pair containing the weighted average of the median distance to the nearest neighbor (first)
 * and the root-mean-squared distance (second) across blocks.
 * These values can be used in `compute_scale()`.
 */
template<typename Index_, typename Input_, typename Distance_, class Matrix_ = knncolle::Matrix<Index_, Input_> >
std::pair<Distance_, Distance_> compute_distance_blocked(
    const std::size_t num_dim,
    const std::vector<Index_> block_sizes,
    const Input_* const data,
    const knncolle::Builder<Index_, Input_, Distance_, Matrix_>& builder,
    const BlockedOptions& options
) {
    const auto num_blocks = block_sizes.size();
    auto prebuilts = sanisizer::create<std::vector<std::shared_ptr<const knncolle::Prebuilt<Index_, Input_, Distance_> > > >(num_blocks);

    Index_ sofar = 0;
    for (I<decltype(num_blocks)> b = 0; b < num_blocks; ++b) {
        const auto cursize = block_sizes[b];
        if (cursize) {
            prebuilts[b] = builder.build_shared(knncolle::SimpleMatrix(num_dim, cursize, data + sanisizer::product_unsafe<std::size_t>(sofar, num_dim)));
        }
        sofar += cursize;
    }

    auto workspace = create_workspace<Distance_>(block_sizes, options);
    return compute_distance_blocked(prebuilts, workspace, options);
}

/**
 * Overload of `compute_distance()` that accepts an embedding matrix with non-contiguous block assignments.
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
 * The number of rows should be equal to `num_dim` and `num_cells`.
 * @param block Pointer to an array of length equal to `num_cells`,
 * containing the block assignment for each column of `data`.
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
    const Block_* const block,
    const knncolle::Builder<Index_, Input_, Distance_, Matrix_>& builder,
    const BlockedOptions& options
) {
    Block_ num_blocks = 0;
    if (num_cells) {
        num_blocks = sanisizer::sum<Block_>(1, *std::max_element(block, block + num_cells));
    }

    auto block_ends = sanisizer::create<std::vector<Index_> >(num_blocks);
    auto block_non_contig = sanisizer::create<std::vector<char> >(num_blocks);
    auto block_sizes = sanisizer::create<std::vector<Index_> >(num_blocks);

    for (Index_ c = 0; c < num_cells; ++c) {
        const auto curb = block[c];
        block_sizes[curb] += 1;

        auto& nc = block_non_contig[curb];
        if (!nc) {
            auto& be = block_ends[curb];
            if (be == 0) {
                be = c + 1;
            } else if (be == c) {
                ++be;
            } else {
                nc = true;
            }
        }
    }

    auto prebuilts = sanisizer::create<std::vector<std::shared_ptr<const knncolle::Prebuilt<Index_, Input_, Distance_> > > >(num_blocks);
    Index_ non_contig_offset = 0;

    for (Block_ b = 0; b < num_blocks; ++b) {
        const auto length = block_sizes[b];
        if (block_non_contig[b]) {
            block_ends[b] = non_contig_offset; // repurposing block_ends to store starting offsets for each non-contiguous block.
            non_contig_offset += length;
        } else if (length) {
            const auto start = block_ends[b] - length;
            const auto ptr = data + sanisizer::product_unsafe<std::size_t>(start, num_dim);
            prebuilts[b] = builder.build_shared(knncolle::SimpleMatrix(num_dim, length, ptr));
        }
    }

    // Handling the non-contiguous blocks.
    if (non_contig_offset > 0) {
        std::vector<Input_> tmp_buffer(sanisizer::product<typename std::vector<Input_>::size_type>(non_contig_offset, num_dim));

        for (Index_ c = 0; c < num_cells; ++c) {
            const auto curb = block[c];
            if (block_non_contig[curb]) {
                auto& curoff = block_ends[curb];
                std::copy_n(
                    data + sanisizer::product_unsafe<std::size_t>(c, num_dim),
                    num_dim,
                    tmp_buffer.data() + sanisizer::product_unsafe<std::size_t>(curoff, num_dim)
                );
                ++curoff;
            }
        }

        non_contig_offset = 0;
        for (Block_ b = 0; b < num_blocks; ++b) {
            if (block_non_contig[b]) {
                const auto length = block_sizes[b];
                const auto ptr = tmp_buffer.data() + sanisizer::product_unsafe<std::size_t>(non_contig_offset, num_dim);
                prebuilts[b] = builder.build_shared(knncolle::SimpleMatrix(num_dim, length, ptr));
                non_contig_offset += length;
            }
        }
    }

    auto workspace = create_workspace<Distance_>(block_sizes, options);
    return compute_distance_blocked(prebuilts, workspace, options);
}

}

#endif
