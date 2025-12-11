#include "get_amideI_geometry.hpp"
#include "helper_vec3.hpp"
#include <iostream>

// -------------------------
// Compute vibration center
// center = C + 0.665 * CO + 0.256 * CN
// -------------------------
static inline Vec3 compute_vibration_center(
        const Vec3& C_coord,
        const Vec3& CO_vector,
        const Vec3& CN_vector)
{
    return C_coord + CO_vector * 0.665 + CN_vector * 0.256;
}

std::vector<AmideIGeo> get_amideI_geometry(
    const std::vector<AmideIEntry>& amideAtoms)
{
    std::vector<AmideIGeo> geo;
    geo.reserve(amideAtoms.size());

    for (const auto& e : amideAtoms) {

        Vec3 C_coord = {e.C.x, e.C.y, e.C.z};
        Vec3 O_coord = {e.O.x, e.O.y, e.O.z};
        Vec3 N_coord = {e.N.x, e.N.y, e.N.z};

        // Compute normalized CO and CN directions
        Vec3 CO_vector = normalize(O_coord - C_coord);
        Vec3 CN_vector = normalize(N_coord - C_coord);

        // Compute vibrational center of dipole
        Vec3 vibration_center_coord =
            compute_vibration_center(C_coord, CO_vector, CN_vector);

        geo.push_back({
            CO_vector,
            CN_vector,
            vibration_center_coord
        });
    }

    return geo;
}
