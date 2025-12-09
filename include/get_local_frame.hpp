#ifndef GET_LOCAL_FRAME_HPP
#define GET_LOCAL_FRAME_HPP

#include <vector>
#include "get_amideI_geometry.hpp"   
#include "helper_vec3.hpp"

// Output struct storing coordinate frame axes
struct LocalFrame {
    Vec3 X_axis;
    Vec3 Y_axis;
    Vec3 Z_axis;     // = CO_vector
};

// API: compute X, Y, Z axes for every amide site
std::vector<LocalFrame> get_local_frame(
    const std::vector<AmideIGeo>& geo);

#endif
