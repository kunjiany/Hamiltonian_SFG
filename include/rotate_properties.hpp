#ifndef ROTATE_PROPERTIES_HPP
#define ROTATE_PROPERTIES_HPP

#include <vector>
#include "helper_vec3.hpp"
#include "get_amideI_properties.hpp"  // defines AmideIProps

// Rotated site properties in the rod / lab frame
struct RotatedProps {
    Vec3  dipole_rot;                // rotated dipole μ'

    // Full rotated Raman tensor α' in 3×3 form (row-major)
    double alpha_rot_matrix[3][3];

    // Rotated Raman tensor α' flattened in MATLAB COLUMN-MAJOR order
    // (this is what StrucInfo.alpha would store)
    double alpha_rot_vectorized[9];

    // Backward-compatible alias: same as alpha_rot_vectorized
    double alpha_rot[9];
};

// Rotate intrinsic properties into the rod/lab frame
std::vector<RotatedProps> rotate_properties(
    const std::vector<AmideIProps> &props_intrinsic,
    double tilt_deg,
    double twist_deg
);

#endif
