#ifndef HAMILTONIAN_EQUIV_MATLAB_HPP
#define HAMILTONIAN_EQUIV_MATLAB_HPP

#include <vector>
#include <array>
#include "helper_vec3.hpp"
#include "get_amideI_geometry.hpp"
#include "get_amideI_properties.hpp"
#include "initialize_amideI_frequency.hpp"
#include "mualphagen.hpp"      

// -----------------------------------------------------------------------------
// Final output that matches MATLAB OneExcitonH.m
// -----------------------------------------------------------------------------
struct HamiltonianEquivResult {
    int N = 0;

    // Sorted exciton frequencies
    std::vector<double> Sort_Ex_Freq;        // size N

    // Sorted eigenvectors (columns sorted with frequencies)
    std::vector<double> Sort_V;              // N×N row-major

    // Rotated site dipoles μ_i (lab frame)
    std::vector<Vec3> mu_rot;                // size N

    // Rotated Raman tensors α_i (lab frame), col-major 3×3 → 9 entries
    std::vector<std::array<double,9>> alpha_rot;  // size N

    // Exciton dipoles μ_k (MATLAB: V' μ_loc V)
    std::vector<Vec3> mu_ex;                      // size N

    // Exciton Raman tensors α_k (MATLAB: V' α_loc V)
    std::vector<std::array<double,9>> alpha_ex;    // size N
};

// -----------------------------------------------------------------------------
// Main MATLAB-equivalent driver
// -----------------------------------------------------------------------------
HamiltonianEquivResult Hamiltonian_equiv_matlab(
    const std::vector<AmideIGeo>& geo,
    const std::vector<AmideIProps>& props,
    const std::vector<AmideIFreq>& freqs,
    double tilt_deg,
    double twist_deg
);

#endif
