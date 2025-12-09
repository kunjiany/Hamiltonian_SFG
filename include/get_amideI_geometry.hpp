#ifndef GET_AMIDEI_GEOMETRY_HPP
#define GET_AMIDEI_GEOMETRY_HPP

#include <vector>
#include "Extract_Amide_Coordinates.hpp"   // defines AmideIEntry
#include "helper_vec3.hpp"


// Stores geometric information needed for later SFG & Hamiltonian calculations
struct AmideIGeo {
    Vec3 CO_vector;                // normalized CO direction
    Vec3 CN_vector;                // normalized CN direction
    Vec3 vibration_center_coord;   // amide-I vibrational center
};

// Main API: compute geometry for each amide residue
std::vector<AmideIGeo> get_amideI_geometry(
    const std::vector<AmideIEntry>& amideAtoms);

#endif
