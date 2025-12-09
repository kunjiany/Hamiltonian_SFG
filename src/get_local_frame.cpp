#include "get_local_frame.hpp"
#include "helper_vec3.hpp"
#include <cmath>
#include <iostream>


std::vector<LocalFrame> get_local_frame(
    const std::vector<AmideIGeo>& geo)
{
    std::vector<LocalFrame> frames;
    frames.reserve(geo.size());

    for(const auto& g : geo)
    {
        // Z-axis = CO direction
        Vec3 Z_axis = g.CO_vector;

        // X-axis = normalized cross(Z, CN)
        Vec3 X_axis = cross(Z_axis, g.CN_vector);
        X_axis = normalize(X_axis);

        // Y-axis = cross(X, Z)
        Vec3 Y_axis = cross(X_axis, Z_axis);

        frames.push_back({X_axis, Y_axis, Z_axis});
    }

    return frames;
}
