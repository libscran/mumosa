#include "scran_tests/scran_tests.hpp"

#include "mumosa/blocked.hpp"
#include "mumosa/compute_scale.hpp"

#include <vector>
#include <random>

class ComputeDistanceBlockedTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        first = scran_tests::simulate_vector(ndim * nobs, scran_tests::SimulationParameters());
        builder.reset(new knncolle::VptreeBuilder<int, double, double>(
            std::make_shared<knncolle::EuclideanDistance<double, double> >()
        ));
    }

    inline static int ndim = 5;
    inline static int nobs = 1234;
    inline static std::vector<double> first;
    inline static std::unique_ptr<knncolle::Builder<int, double, double> > builder;
};

TEST_F(ComputeDistanceBlockedTest, Basic) {
    auto combined = first;
    combined.reserve(first.size() * 2);
    for (auto s : first) {
        combined.push_back(s * 2);
    }

    {
        auto ref = mumosa::compute_distance(ndim, nobs, first.data(), *builder, mumosa::Options());
        auto out = mumosa::compute_distance_blocked(ndim, { nobs, nobs }, combined.data(), *builder, mumosa::BlockedOptions());
        EXPECT_FLOAT_EQ(mumosa::compute_scale(ref, out), 2.0/3);
    }

    // Works for other options.
    {
        mumosa::Options opt;
        opt.num_neighbors = 10;
        auto ref = mumosa::compute_distance(ndim, nobs, first.data(), *builder, opt);

        mumosa::BlockedOptions bopt;
        bopt.num_neighbors = 10;
        auto out = mumosa::compute_distance_blocked(ndim, { nobs, nobs }, combined.data(), *builder, bopt);

        EXPECT_FLOAT_EQ(mumosa::compute_scale(ref, out), 2.0/3);
    }
}

TEST_F(ComputeDistanceBlockedTest, BlockWeights) {
    std::vector<int> sizes { 100, 700, nobs - 800 };

    auto ref1 = mumosa::compute_distance(ndim, sizes[0], first.data(), *builder, mumosa::Options());
    auto ref2 = mumosa::compute_distance(ndim, sizes[1], first.data() + sanisizer::product_unsafe<std::size_t>(sizes[0], ndim), *builder, mumosa::Options());
    auto ref3 = mumosa::compute_distance(ndim, sizes[2], first.data() + sanisizer::product_unsafe<std::size_t>(sizes[0] + sizes[1], ndim), *builder, mumosa::Options());

    // By default, weighted by the number of observations if we don't hit 1000.
    {
        auto out = mumosa::compute_distance_blocked(ndim, sizes, first.data(), *builder, mumosa::BlockedOptions());
        EXPECT_FLOAT_EQ(out.first, (ref1.first * sizes[0] + ref2.first * sizes[1] + ref3.first * sizes[2]) / nobs); 
        EXPECT_FLOAT_EQ(out.second, (ref1.second * sizes[0] + ref2.second * sizes[1] + ref3.second * sizes[2]) / nobs); 
    }

    // What happens if we reduce the variable threshold?
    {
        mumosa::BlockedOptions bopt;
        bopt.variable_block_weight_parameters.upper_bound = 200;
        auto out = mumosa::compute_distance_blocked(ndim, sizes, first.data(), *builder, bopt);
        EXPECT_FLOAT_EQ(out.first, (ref1.first * 0.5 + ref2.first + ref3.first) / 2.5); 
        EXPECT_FLOAT_EQ(out.second, (ref1.second * 0.5 + ref2.second + ref3.second) / 2.5); 
    }

    // What happens if we set it to EQUAL weights?
    {
        mumosa::BlockedOptions bopt;
        bopt.block_weight_policy = scran_blocks::WeightPolicy::EQUAL;
        auto out = mumosa::compute_distance_blocked(ndim, sizes, first.data(), *builder, bopt);
        EXPECT_FLOAT_EQ(out.first, (ref1.first + ref2.first + ref3.first) / 3); 
        EXPECT_FLOAT_EQ(out.second, (ref1.second + ref2.second + ref3.second) / 3); 
    }
}

TEST_F(ComputeDistanceBlockedTest, Empty) {
    // Partially empty.
    {
        auto out = mumosa::compute_distance_blocked(ndim, { 0, nobs, 0 }, first.data(), *builder, mumosa::BlockedOptions());
        auto ref = mumosa::compute_distance(ndim, nobs, first.data(), *builder, mumosa::Options());
        EXPECT_EQ(ref, out);
    }

    // Fully empty.
    {
        auto out = mumosa::compute_distance_blocked(ndim, { 0, 0, 0 }, first.data(), *builder, mumosa::BlockedOptions());
        EXPECT_EQ(out.first, 0);
        EXPECT_EQ(out.second, 0);
    }
}

TEST_F(ComputeDistanceBlockedTest, Rearranged) {
    std::vector<int> sizes { 100, 200, 300, 400, nobs - 1000 };
    const std::size_t num_blocks = sizes.size();

    std::vector<char> blocks;
    int counter = 0;
    for (auto s : sizes) {
        blocks.insert(blocks.end(), s, counter);
        ++counter;
    }
    auto ref = mumosa::compute_distance_blocked(ndim, sizes, first.data(), *builder, mumosa::BlockedOptions());

    {
        auto out = mumosa::compute_distance_blocked(ndim, nobs, first.data(), blocks.data(), num_blocks, *builder, mumosa::BlockedOptions());
        EXPECT_EQ(ref, out);
    }

    // Multiplying block assignments by 2 to force the existence of empty (odd-numbered) clusters. 
    {
        auto blocks2 = blocks;
        for (auto& b : blocks2) {
            b *= 2;
        }
        auto out = mumosa::compute_distance_blocked(ndim, nobs, first.data(), blocks2.data(), num_blocks * 2, *builder, mumosa::BlockedOptions());
        EXPECT_EQ(ref, out);
    }

    // Swapping a stretch of observations so that the second and fourth blocks are not contiguous.
    {
        auto interspersed_first = first;
        auto interspersed_blocks = blocks;

        std::copy_n(
            first.begin() + sanisizer::product_unsafe<std::size_t>(200, ndim),
            sanisizer::product_unsafe<std::size_t>(100, ndim),
            interspersed_first.begin() + sanisizer::product_unsafe<std::size_t>(600, ndim)
        );
        std::copy_n(
            blocks.begin() + 200,
            100,
            interspersed_blocks.begin() + 600
        );

        std::copy_n(
            first.begin() + sanisizer::product_unsafe<std::size_t>(600, ndim),
            sanisizer::product_unsafe<std::size_t>(100, ndim),
            interspersed_first.begin() + sanisizer::product_unsafe<std::size_t>(200, ndim)
        );
        std::copy_n(
            blocks.begin() + 600,
            100,
            interspersed_blocks.begin() + 200 
        );

        auto out = mumosa::compute_distance_blocked(ndim, nobs, interspersed_first.data(), interspersed_blocks.data(), num_blocks, *builder, mumosa::BlockedOptions());
        EXPECT_EQ(ref, out); // equality assumes that the swap does not change the order of observations within each block.
    }

    // Randomizing all of the observations.
    auto shuffled_blocks = blocks;

    std::mt19937_64 rng(23423);
    std::shuffle(shuffled_blocks.begin(), shuffled_blocks.end(), rng);

    auto offsets = sizes;
    int cumulative = 0;
    for (auto& o : offsets) {
        auto previous = cumulative;
        cumulative += o;
        o = previous;
    }

    std::vector<double> shuffled_first;
    shuffled_first.reserve(first.size());
    for (auto b : shuffled_blocks) {
        auto& off = offsets[b];
        auto start = first.begin() + sanisizer::product_unsafe<std::size_t>(off, ndim);
        shuffled_first.insert(shuffled_first.end(), start, start + ndim);
        ++off;
    }

    {
        auto out = mumosa::compute_distance_blocked(ndim, nobs, shuffled_first.data(), shuffled_blocks.data(), num_blocks, *builder, mumosa::BlockedOptions());
        EXPECT_EQ(ref, out);
    }

    // Multiplying block assignments to force the existence of empty (even-numbered) clusters. 
    {
        auto shuffled_blocks2 = shuffled_blocks;
        for (auto& b : shuffled_blocks2) {
            b = b * 2 + 1;
        }
        auto out = mumosa::compute_distance_blocked(ndim, nobs, shuffled_first.data(), shuffled_blocks2.data(), num_blocks * 2, *builder, mumosa::BlockedOptions());
        EXPECT_EQ(ref, out);
    }
}
