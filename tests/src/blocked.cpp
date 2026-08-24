#include "scran_tests/scran_tests.hpp"

#include "mumosa/blocked.hpp"
#include "mumosa/compute_scale.hpp"

#include <vector>
#include <random>

class ComputeDistanceBlockedTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        first = scran_tests::simulate_vector(ndim * nobs, scran_tests::SimulateVectorParameters());
        builder.reset(new knncolle::VptreeBuilder<int, double, double>(
            std::make_shared<knncolle::EuclideanDistance<double, double> >()
        ));
    }

    inline static int ndim = 5;
    inline static int nobs = 1234;
    inline static std::vector<double> first;
    inline static std::unique_ptr<knncolle::Builder<int, double, double> > builder;
};

template<typename Block_ = int>
std::vector<Block_> sizes_to_block(const std::vector<int>& sizes) {
    std::vector<Block_> output;
    const std::size_t nblocks = sizes.size();
    for (std::size_t b = 0; b < nblocks; ++b) {
        output.insert(output.end(), sizes[b], b);
    }
    return output;
}

TEST_F(ComputeDistanceBlockedTest, Basic) {
    auto combined = first;
    combined.reserve(first.size() * 2);
    for (auto s : first) {
        combined.push_back(s * 2);
    }

    {
        std::vector<double> buffer(nobs);
        auto ref = mumosa::compute_distance(ndim, nobs, first.data(), *builder, buffer.data(), mumosa::Options());

        auto block = sizes_to_block({ nobs, nobs });
        buffer.resize(nobs * 2);
        auto out = mumosa::compute_distance_blocked(ndim, nobs * 2, combined.data(), block.data(), 2, *builder, buffer.data(), mumosa::BlockedOptions());
        EXPECT_FLOAT_EQ(mumosa::compute_scale(ref, out), 2.0/3);
    }

    // Works for other options.
    {
        mumosa::Options opt;
        opt.num_neighbors = 10;
        std::vector<double> buffer(nobs);
        auto ref = mumosa::compute_distance(ndim, nobs, first.data(), *builder, buffer.data(), opt);

        mumosa::BlockedOptions bopt;
        bopt.num_neighbors = 10;
        auto block = sizes_to_block({ nobs, nobs });
        buffer.resize(nobs * 2);
        auto out = mumosa::compute_distance_blocked(ndim, nobs * 2, combined.data(), block.data(), 2, *builder, buffer.data(), bopt);

        EXPECT_FLOAT_EQ(mumosa::compute_scale(ref, out), 2.0/3);
    }
}

TEST_F(ComputeDistanceBlockedTest, BlockWeights) {
    std::vector<int> sizes { 100, 700, nobs - 800 };
    std::vector<double> buffer(nobs);

    auto ref1 = mumosa::compute_distance(ndim, sizes[0], first.data(), *builder, buffer.data(), mumosa::Options());
    auto ref2 = mumosa::compute_distance(ndim, sizes[1], first.data() + sanisizer::product_unsafe<std::size_t>(sizes[0], ndim), *builder, buffer.data(), mumosa::Options());
    auto ref3 = mumosa::compute_distance(ndim, sizes[2], first.data() + sanisizer::product_unsafe<std::size_t>(sizes[0] + sizes[1], ndim), *builder, buffer.data(), mumosa::Options());

    // By default, weighted by the number of observations if we don't hit 1000.
    const std::size_t num_blocks = sizes.size();
    auto block = sizes_to_block(sizes);
    {
        auto out = mumosa::compute_distance_blocked(ndim, nobs, first.data(), block.data(), num_blocks, *builder, buffer.data(), mumosa::BlockedOptions());
        EXPECT_FLOAT_EQ(out.first, (ref1.first * sizes[0] + ref2.first * sizes[1] + ref3.first * sizes[2]) / nobs); 
        EXPECT_FLOAT_EQ(out.second, (ref1.second * sizes[0] + ref2.second * sizes[1] + ref3.second * sizes[2]) / nobs); 
    }

    // What happens if we reduce the variable threshold?
    {
        mumosa::BlockedOptions bopt;
        bopt.variable_block_weight_parameters.upper_bound = 200;
        auto out = mumosa::compute_distance_blocked(ndim, nobs, first.data(), block.data(), num_blocks, *builder, buffer.data(), bopt);
        EXPECT_FLOAT_EQ(out.first, (ref1.first * 0.5 + ref2.first + ref3.first) / 2.5); 
        EXPECT_FLOAT_EQ(out.second, (ref1.second * 0.5 + ref2.second + ref3.second) / 2.5); 
    }

    // What happens if we set it to EQUAL weights?
    {
        mumosa::BlockedOptions bopt;
        bopt.block_weight_policy = scran_blocks::WeightPolicy::EQUAL;
        auto out = mumosa::compute_distance_blocked(ndim, nobs, first.data(), block.data(), num_blocks, *builder, buffer.data(), bopt);
        EXPECT_FLOAT_EQ(out.first, (ref1.first + ref2.first + ref3.first) / 3); 
        EXPECT_FLOAT_EQ(out.second, (ref1.second + ref2.second + ref3.second) / 3); 
    }
}

TEST_F(ComputeDistanceBlockedTest, EmptyBlocks) {
    std::vector<double> buffer(nobs);

    // Partially empty; we consider three blocks, but all cells are in the second block.
    {
        std::vector<int> block(nobs, 1);
        auto out = mumosa::compute_distance_blocked(ndim, nobs, first.data(), block.data(), 3, *builder, buffer.data(), mumosa::BlockedOptions());
        auto ref = mumosa::compute_distance(ndim, nobs, first.data(), *builder, buffer.data(), mumosa::Options());
        EXPECT_EQ(ref, out);
    }

    // Multiplying block assignments to force the existence of empty (even-numbered) clusters. 
    {
        std::vector<int> sizes { 500, 100, nobs - 600 };
        const std::size_t num_blocks = sizes.size();

        auto block = sizes_to_block<char>(sizes);
        auto ref = mumosa::compute_distance_blocked(ndim, nobs, first.data(), block.data(), num_blocks, *builder, buffer.data(), mumosa::BlockedOptions());

        for (auto& b : block) {
            b = b * 2 + 1;
        }
        auto out = mumosa::compute_distance_blocked(ndim, nobs, first.data(), block.data(), num_blocks * 2 + 1, *builder, buffer.data(), mumosa::BlockedOptions());
        EXPECT_EQ(ref, out);
    }

    // Fully empty.
    {
        auto out = mumosa::compute_distance_blocked(ndim, 0, first.data(), static_cast<int*>(NULL), 0, *builder, buffer.data(), mumosa::BlockedOptions());
        EXPECT_EQ(out.first, 0);
        EXPECT_EQ(out.second, 0);
    }
}

TEST_F(ComputeDistanceBlockedTest, FewPoints) {
    // One point.
    // Check that we avoid indexing the end of an empty distance vector.
    {
        std::vector<int> blocks { 1, 2, 0 };
        std::vector<double> buffer(blocks.size());
        auto alt = mumosa::compute_distance_blocked(ndim, static_cast<int>(blocks.size()), first.data(), blocks.data(), 3, *builder, buffer.data(), {});
        EXPECT_EQ(alt.first, 0);
        EXPECT_EQ(alt.second, 0);
    }

    // Two points.
    // Check that 'k' is properly capped.
    {
        std::vector<int> blocks { 1, 0, 2, 1, 2, 0 };
        std::vector<double> buffer(blocks.size());
        auto alt = mumosa::compute_distance_blocked(ndim, static_cast<int>(blocks.size()), first.data(), blocks.data(), 3, *builder, buffer.data(), {});
        EXPECT_GT(alt.first, 0);
        EXPECT_GT(alt.second, 0);
    }
}

TEST_F(ComputeDistanceBlockedTest, Rearranged) {
    std::vector<int> sizes { 100, 200, 300, 400, nobs - 1000 };
    const std::size_t num_blocks = sizes.size();

    auto block = sizes_to_block<char>(sizes);
    std::vector<double> buffer(nobs);
    auto ref = mumosa::compute_distance_blocked(ndim, nobs, first.data(), block.data(), num_blocks, *builder, buffer.data(), mumosa::BlockedOptions());

    // Swapping a stretch of observations so that the second and fourth blocks are not contiguous.
    {
        auto interspersed_first = first;
        auto interspersed_block = block;

        std::copy_n(
            first.begin() + sanisizer::product_unsafe<std::size_t>(200, ndim),
            sanisizer::product_unsafe<std::size_t>(100, ndim),
            interspersed_first.begin() + sanisizer::product_unsafe<std::size_t>(600, ndim)
        );
        std::copy_n(
            block.begin() + 200,
            100,
            interspersed_block.begin() + 600
        );

        std::copy_n(
            first.begin() + sanisizer::product_unsafe<std::size_t>(600, ndim),
            sanisizer::product_unsafe<std::size_t>(100, ndim),
            interspersed_first.begin() + sanisizer::product_unsafe<std::size_t>(200, ndim)
        );
        std::copy_n(
            block.begin() + 600,
            100,
            interspersed_block.begin() + 200 
        );

        auto out = mumosa::compute_distance_blocked(
            ndim,
            nobs,
            interspersed_first.data(),
            interspersed_block.data(),
            num_blocks,
            *builder,
            buffer.data(),
            mumosa::BlockedOptions()
        );
        EXPECT_EQ(ref, out); // equality assumes that the swap does not change the order of observations within each block.
    }

    // Randomizing all of the observations.
    {
        auto shuffled_block = block;

        std::mt19937_64 rng(23423);
        std::shuffle(shuffled_block.begin(), shuffled_block.end(), rng);

        auto offsets = sizes;
        int cumulative = 0;
        for (auto& o : offsets) {
            auto previous = cumulative;
            cumulative += o;
            o = previous;
        }

        std::vector<double> shuffled_first;
        shuffled_first.reserve(first.size());
        for (auto b : shuffled_block) {
            auto& off = offsets[b];
            auto start = first.begin() + sanisizer::product_unsafe<std::size_t>(off, ndim);
            shuffled_first.insert(shuffled_first.end(), start, start + ndim);
            ++off;
        }

        auto out = mumosa::compute_distance_blocked(
            ndim,
            nobs,
            shuffled_first.data(),
            shuffled_block.data(),
            num_blocks,
            *builder, 
            buffer.data(),
            mumosa::BlockedOptions()
        );
        EXPECT_EQ(ref, out);
    }
}
