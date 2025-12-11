#ifndef COMPUTE_DIPOLE_COUPLING_HPP
#define COMPUTE_DIPOLE_COUPLING_HPP

#include <vector>
#include "helper_vec3.hpp"

double compute_dipole_coupling(
    const Vec3 &mu_i,
    const Vec3 &mu_j,
    const Vec3 &Ri_minus_Rj,
    double distance_cutoff,   // Å
    bool use_cutoff,
    double prefactor          // MATLAB beta prefactor (cm⁻¹·Å³)
);

#endif

