#ifndef GET_AMIDEI_PROPERTIES_HPP
#define GET_AMIDEI_PROPERTIES_HPP

#include <vector>
#include "helper_vec3.hpp"
#include "get_local_frame.hpp"

// Stores rotated dipole and Raman tensor information
struct AmideIProps {
    Vec3 dipole_sim;
    double alpha_matrix[3][3];    // full tensor
    double alpha_vectorized[9];   // 9 elements (MATLAB reshape order)
    double alpha_reduced[6];      // XX YY ZZ XY YZ XZ
};


// API: compute dipole and Raman properties for all amide modes
std::vector<AmideIProps> get_amideI_properties(
    const std::vector<LocalFrame>& frames);

#endif
