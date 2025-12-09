#include "chi2_matlab.hpp"
#include "apply_R3.hpp"
#include <iostream>
#include <iomanip>

// -----------------------------------------------------------------------------
// MATLAB-equivalent outer product:
// χ_mol = alpha(:) * mu(:)'   → 27×1 vector
//
// MATLAB's linear order is COLUMN-MAJOR:
//
// alpha(:) = [xx; yx; zx; xy; yy; zy; xz; yz; zz]
// mu(:)    = [mx; my; mz]
//
// For each μ-component (slow index b),
// multiply all 9 α components (fast index a)
// -----------------------------------------------------------------------------
static std::array<double,27> outer_alpha_mu(
    const std::array<double,9>& alpha,
    const Vec3& mu)
{
    std::array<double,27> out{};
    double m[3] = { mu.x, mu.y, mu.z };

    int k = 0;
    for (int b = 0; b < 3; b++)      // μ index
        for (int a = 0; a < 9; a++)  // α index
            out[k++] = alpha[a] * m[b];

    return out;
}


// -----------------------------------------------------------------------------
// MATLAB-equivalent χ² computation
// -----------------------------------------------------------------------------
Chi2Result compute_chi2_matlab(
    const HamiltonianEquivResult& H,
    const R3Matrix& R)
{
    Chi2Result out;
    int N = H.N;
    if (N == 0)
    {
        std::cerr << "[chi2_matlab] ERROR: N == 0\n";
        return out;
    }

    out.N = N;
    out.freq.resize(N);
    out.chi_mol.resize(N);
    out.chi_lab.resize(N);



    // Loop through exciton states
    for (int k = 0; k < N; k++)
    {
        // 1. Frequency
        out.freq[k] = H.Sort_Ex_Freq[k];

        // 2. Extract exciton μ and α
        const Vec3& mu_ex = H.mu_ex[k];
        const std::array<double,9>& alpha_ex = H.alpha_ex[k];

        // 3. Compute χ_mol (27×1)
        out.chi_mol[k] = outer_alpha_mu(alpha_ex, mu_ex);

        // 4. Rotate to lab frame
        out.chi_lab[k] = apply_R3_single(R, out.chi_mol[k]);
    }



    return out;
}
