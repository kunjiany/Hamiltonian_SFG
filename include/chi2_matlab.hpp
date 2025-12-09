#ifndef CHI2_MATLAB_HPP
#define CHI2_MATLAB_HPP

#include <array>
#include <vector>
#include "helper_vec3.hpp"
#include "hamiltonian_equiv_matlab.hpp"

struct R3Matrix;
// ---------------------------------------------------------
// Unified χ² output structure
// ---------------------------------------------------------
struct Chi2Result {
    int N = 0;

    std::vector<double> freq;                    // size N
    std::vector<std::array<double,27>> chi_mol;  // χ(mol) size N
    std::vector<std::array<double,27>> chi_lab;  // χ(lab) size N
};

// ---------------------------------------------------------
// Compute χ² using:
//   - H    → Hamiltonian result (mu_ex, alpha_ex, sorted eigenvalues)
//   - R    → 27×27 rotation matrix loaded from HDF5 database
// ---------------------------------------------------------
Chi2Result compute_chi2_matlab(
    const HamiltonianEquivResult& H,
    const R3Matrix& R);

#endif

