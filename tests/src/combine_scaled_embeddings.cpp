#include "scran_tests/scran_tests.hpp"

#include "mumosa/combine_scaled_embeddings.hpp"

#include <vector>

TEST(CombineScaledEmbeddings, Basic) {
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
        mumosa::combine_scaled_embeddings(
            { static_cast<std::size_t>(20), static_cast<std::size_t>(5) },
            nobs,
            std::vector<double*>{ first.data(), second.data() },
            std::vector<double>{ 0.5, 1.2 },
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
        mumosa::combine_scaled_embeddings(
            { static_cast<std::size_t>(20), static_cast<std::size_t>(5) },
            nobs,
            std::vector<double*>{ first.data(), second.data() },
            std::vector<double>{ 0.5, std::numeric_limits<double>::infinity() },
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
