#include "apply_R3.hpp"

std::array<double,27> apply_R3_single(
    const R3Matrix& R,
    const std::array<double,27>& chi_in)
{
    std::array<double,27> chi_out{};

    for (int i = 0; i < 27; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < 27; j++)
            sum += R.m[i][j] * chi_in[j];

        chi_out[i] = sum;
    }

    return chi_out;
}
