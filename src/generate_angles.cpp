#include "generate_angles.hpp"

std::vector<double> Linspace(double start, double end, int points)
{
    std::vector<double> v;

    if(points <= 0)
        return v;

    // If only one sample is needed:
    if(points == 1) {
        v.push_back(start);
        return v;
    }

    v.resize(points);
    double step = (end - start) / (points - 1);

    for(int i = 0; i < points; i++)
        v[i] = start + step * i;

    return v;
}
