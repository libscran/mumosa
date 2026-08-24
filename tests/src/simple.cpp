#include "scran_tests/scran_tests.hpp"

#include "mumosa/simple.hpp"
#include "mumosa/compute_scale.hpp"

#include <vector>

class ComputeDistanceTest : public ::testing::Test {
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

TEST_F(ComputeDistanceTest, Basic) {
    auto second = first;
    for (auto& s : second) {
        s *= 2;
    }

    auto dbuffer = sanisizer::create<std::vector<double> >(nobs);
    auto out1 = mumosa::compute_distance(ndim, nobs, first.data(), *builder, dbuffer.data(), mumosa::Options());
    auto out2 = mumosa::compute_distance(ndim, nobs, second.data(), *builder, dbuffer.data(), mumosa::Options());

    EXPECT_FLOAT_EQ(out1.first / out2.first, 0.5);
    EXPECT_FLOAT_EQ(out1.second / out2.second , 0.5);
    EXPECT_FLOAT_EQ(mumosa::compute_scale(out1, out2), 0.5);

    // Works in parallel.
    {
        mumosa::Options opt;
        opt.num_threads = 3;
        auto pout1 = mumosa::compute_distance(ndim, nobs, first.data(), *builder, dbuffer.data(), opt);
        EXPECT_EQ(out1, pout1);
    }

    // Works for other options.
    {
        mumosa::Options opt;
        opt.num_neighbors = 10;
        auto out10_1 = mumosa::compute_distance(ndim, nobs, first.data(), *builder, dbuffer.data(), opt);
        auto out10_2 = mumosa::compute_distance(ndim, nobs, second.data(), *builder, dbuffer.data(), opt);
        EXPECT_LT(out10_1, out1);
        EXPECT_LT(out10_2, out2);
        EXPECT_FLOAT_EQ(mumosa::compute_scale(out10_1, out10_2), 0.5);
    }
}

TEST_F(ComputeDistanceTest, DifferentlyDimensioned) {
    std::vector<double> second(ndim * 2 * nobs);
    auto fIt = first.begin();
    auto sIt = second.begin();
    for (int o = 0; o < nobs; ++o) {
        std::copy(fIt, fIt + ndim, sIt);
        sIt += ndim;
        std::copy(fIt, fIt + ndim, sIt);
        fIt += ndim;
        sIt += ndim;
    }

    auto dbuffer = sanisizer::create<std::vector<double> >(nobs);
    auto out1 = mumosa::compute_distance(ndim, nobs, first.data(), *builder, dbuffer.data(), mumosa::Options());
    auto out2 = mumosa::compute_distance(ndim * 2, nobs, second.data(), *builder, dbuffer.data(), mumosa::Options());

    const double expected_ratio = 1.0 / std::sqrt(2);
    EXPECT_FLOAT_EQ(out1.first / out2.first, expected_ratio);
    EXPECT_FLOAT_EQ(out1.second / out2.second , expected_ratio);
    EXPECT_FLOAT_EQ(mumosa::compute_scale(out1, out2), expected_ratio);
}

TEST_F(ComputeDistanceTest, Zeros) {
    auto dbuffer = sanisizer::create<std::vector<double> >(nobs);

    // Switches to the RMSD.
    {
        std::vector<double> second(ndim * nobs);
        second[0] = 1;

        auto out1 = mumosa::compute_distance(ndim, nobs, first.data(), *builder, dbuffer.data(), mumosa::Options());
        auto out2 = mumosa::compute_distance(ndim, nobs, second.data(), *builder, dbuffer.data(), mumosa::Options());
        auto scale = mumosa::compute_scale(out1, out2);

        EXPECT_FALSE(std::isinf(scale));
        EXPECT_GT(scale, 0);
    }

    // Falls back to the edge cases.
    {
        std::vector<double> second(ndim * nobs);

        auto out1 = mumosa::compute_distance(ndim, nobs, first.data(), *builder, dbuffer.data(), mumosa::Options());
        auto out2 = mumosa::compute_distance(ndim, nobs, second.data(), *builder, dbuffer.data(), mumosa::Options());

        auto scale = mumosa::compute_scale(out1, out2);
        EXPECT_TRUE(std::isinf(scale));

        scale = mumosa::compute_scale(out2, out1);
        EXPECT_EQ(scale, 0);
    }
}

TEST_F(ComputeDistanceTest, FewPoints) {
    // No points.
    // Check that we avoid returning NaNs.
    {
        auto alt = mumosa::compute_distance(ndim, 0, static_cast<double*>(NULL), *builder, static_cast<double*>(NULL), mumosa::Options());
        EXPECT_EQ(alt.first, 0);
        EXPECT_EQ(alt.second, 0);
    }

    // One point.
    // Check that we avoid indexing the end of an empty distance vector.
    {
        std::vector<double> buffer(1);
        auto alt = mumosa::compute_distance(ndim, 1, first.data(), *builder, buffer.data(), mumosa::Options());
        EXPECT_EQ(alt.first, 0);
        EXPECT_EQ(alt.second, 0);
    }

    // Two points.
    // Check that 'k' is properly capped.
    {
        std::vector<double> buffer(2);
        auto alt = mumosa::compute_distance(ndim, 2, first.data(), *builder, buffer.data(), mumosa::Options());
        EXPECT_GT(alt.first, 0);
        EXPECT_GT(alt.second, 0);
    }
}

TEST(ComputeDistance, ComputeDistances) {
    {
        std::vector<std::pair<double, double> > distances{ {3, 3}, { 2, 2 }, { 1, 1 } };
        auto output = mumosa::compute_scale(distances);
        std::vector<double> expected { 1, 1.5, 3 };
        EXPECT_EQ(output, expected);
    }

    // Skips the first.
    {
        std::vector<std::pair<double, double> > distances{ { 0, 0 }, { 10, 10 }, { 1, 1 } };
        auto output = mumosa::compute_scale(distances);
        EXPECT_TRUE(std::isinf(output[0]));
        EXPECT_EQ(output[1], 1);
        EXPECT_EQ(output[2], 10);
    }

    // Skips all of them.
    {
        std::vector<std::pair<double, double> > distances(3);
        auto output = mumosa::compute_scale(distances);
        std::vector<double> expected(3);
        EXPECT_EQ(output, expected);

    }
}
