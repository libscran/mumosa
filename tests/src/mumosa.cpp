#include "scran_tests/scran_tests.hpp"

#include "mumosa/mumosa.hpp"

#include <vector>

class ScaleByNeighborsTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        first = scran_tests::simulate_vector(ndim * nobs, scran_tests::SimulationParameters());
    }

    inline static int ndim = 5;
    inline static int nobs = 1234;
    inline static std::vector<double> first;
};

TEST_F(ScaleByNeighborsTest, Basic) {
    auto second = first;
    for (auto& s : second) {
        s *= 2;
    }

    auto out1 = mumosa::compute_distance(ndim, nobs, first.data(), knncolle::VptreeBuilder(), mumosa::Options());
    auto out2 = mumosa::compute_distance(ndim, nobs, second.data(), knncolle::VptreeBuilder(), mumosa::Options());
    EXPECT_FLOAT_EQ(mumosa::compute_scale(out1, out2), 0.5);

    // Works in parallel.
    {
        mumosa::Options opt;
        opt.num_threads = 3;
        auto pout1 = mumosa::compute_distance(ndim, nobs, first.data(), knncolle::VptreeBuilder(), opt);
        EXPECT_EQ(out1, pout1);
    }

    // Works for other options.
    {
        mumosa::Options opt;
        opt.num_neighbors = 10;
        auto out10_1 = mumosa::compute_distance(ndim, nobs, first.data(), knncolle::VptreeBuilder(), opt);
        auto out10_2 = mumosa::compute_distance(ndim, nobs, second.data(), knncolle::VptreeBuilder(), opt);
        EXPECT_LT(out10_1, out1);
        EXPECT_LT(out10_2, out2);
        EXPECT_FLOAT_EQ(mumosa::compute_scale(out10_1, out10_2), 0.5);
    }
}

TEST_F(ScaleByNeighborsTest, DifferentlyDimensioned) {
    std::vector<double> second(ndim*2*nobs);
    auto fIt = first.begin();
    auto sIt = second.begin();
    for (int o = 0; o < nobs; ++o) {
        std::copy(fIt, fIt + ndim, sIt);
        sIt += ndim;
        std::copy(fIt, fIt + ndim, sIt);
        fIt += ndim;
        sIt += ndim;
    }

    auto out1 = mumosa::compute_distance(ndim, nobs, first.data(), knncolle::VptreeBuilder(), mumosa::Options());
    auto out2 = mumosa::compute_distance(ndim * 2, nobs, second.data(), knncolle::VptreeBuilder(), mumosa::Options());
    EXPECT_FLOAT_EQ(mumosa::compute_scale(out1, out2), 1.0 / std::sqrt(2));
}

TEST_F(ScaleByNeighborsTest, Zeros) {
    // Switches to the RMSD.
    {
        std::vector<double> second(ndim * nobs);
        second[0] = 1;

        auto out1 = mumosa::compute_distance(ndim, nobs, first.data(), knncolle::VptreeBuilder(), mumosa::Options());
        auto out2 = mumosa::compute_distance(ndim, nobs, second.data(), knncolle::VptreeBuilder(), mumosa::Options());
        auto scale = mumosa::compute_scale(out1, out2);

        EXPECT_FALSE(std::isinf(scale));
        EXPECT_GT(scale, 0);
    }

    // Falls back to the edge cases.
    {
        std::vector<double> second(ndim * nobs);

        auto out1 = mumosa::compute_distance(ndim, nobs, first.data(), knncolle::VptreeBuilder(), mumosa::Options());
        auto out2 = mumosa::compute_distance(ndim, nobs, second.data(), knncolle::VptreeBuilder(), mumosa::Options());

        auto scale = mumosa::compute_scale(out1, out2);
        EXPECT_TRUE(std::isinf(scale));

        scale = mumosa::compute_scale(out2, out1);
        EXPECT_EQ(scale, 0);
    }
}

TEST(ScaleByNeighbors, ComputeDistances) {
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

TEST(ScaleByNeighbors, CombineEmbeddings) {
    size_t nobs = 123;
    auto first = scran_tests::simulate_vector(20 * nobs, [&]{
        scran_tests::SimulationParameters sparams;
        sparams.seed = 1000;
        return sparams;
    }());
    auto second = scran_tests::simulate_vector(5 * nobs, [&]{
        scran_tests::SimulationParameters sparams;
        sparams.seed = 2000;
        return sparams;
    }());

    {
        std::vector<double> output(25 * nobs);
        mumosa::combine_scaled_embeddings<int, int, double, double, double>(
            { 20, 5 },
            nobs,
            std::vector<double*>{ first.data(), second.data() },
            { 0.5, 1.2 },
            output.data()
        );

        // Interleaving is done correctly.
        EXPECT_EQ(output[0], first[0] * 0.5);
        EXPECT_EQ(output[19], first[19] * 0.5);
        EXPECT_EQ(output[25], first[20] * 0.5);
        EXPECT_EQ(output[25 * (nobs - 1)], first[20 * (nobs - 1)] * 0.5);
        EXPECT_EQ(output[25 * (nobs - 1) + 19], first[20 * nobs - 1] * 0.5);

        EXPECT_EQ(output[20], second[0] * 1.2);
        EXPECT_EQ(output[24], second[4] * 1.2);
        EXPECT_EQ(output[45], second[5] * 1.2);
        EXPECT_EQ(output[25 * (nobs - 1) + 20], second[5 * (nobs - 1)] * 1.2);
        EXPECT_EQ(output[25 * nobs - 1], second[5 * nobs - 1] * 1.2);
    }

    // Handles the infinite special case.
    {
        std::vector<double> output(25 * nobs);
        mumosa::combine_scaled_embeddings<int, int, double, double, double>(
            { 20, 5 },
            nobs,
            std::vector<double*>{ first.data(), second.data() },
            { 0.5, std::numeric_limits<double>::infinity() },
            output.data()
        );

        // Interleaving is done correctly.
        EXPECT_EQ(output[0], first[0] * 0.5);
        EXPECT_EQ(output[19], first[19] * 0.5);
        EXPECT_EQ(output[25], first[20] * 0.5);
        EXPECT_EQ(output[25 * (nobs - 1)], first[20 * (nobs - 1)] * 0.5);
        EXPECT_EQ(output[25 * (nobs - 1) + 19], first[20 * nobs - 1] * 0.5);

        EXPECT_EQ(output[20], 0);
        EXPECT_EQ(output[24], 0);
        EXPECT_EQ(output[45], 0);
        EXPECT_EQ(output[25 * (nobs - 1) + 20], 0);
        EXPECT_EQ(output[25 * nobs - 1], 0);
    }
}
