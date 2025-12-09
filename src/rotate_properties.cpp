#include "rotate_properties.hpp"

std::vector<RotatedProps> rotate_properties(
    const std::vector<AmideIProps> &props_intrinsic,
    double tilt_deg,
    double twist_deg)
{
    std::vector<RotatedProps> out(props_intrinsic.size());

    for (size_t i = 0; i < props_intrinsic.size(); i++)
    {
        // NO ROTATION (match MATLAB)
        out[i].dipole_rot = props_intrinsic[i].dipole_sim;

        // α unchanged
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                out[i].alpha_rot_matrix[r][c] = props_intrinsic[i].alpha_matrix[r][c];

        // Flatten (MATLAB column-major)
        int k = 0;
        for (int col = 0; col < 3; col++)
            for (int row = 0; row < 3; row++)
                out[i].alpha_rot_vectorized[k++] =
                    out[i].alpha_rot_matrix[row][col];

        // alias
        for (int t = 0; t < 9; t++)
            out[i].alpha_rot[t] = out[i].alpha_rot_vectorized[t];
    }

    return out;
}
