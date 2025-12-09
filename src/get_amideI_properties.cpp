#include "get_amideI_properties.hpp"
#include <cmath>

// -----------------------------------------------------------
// Define the intrinsic molecular dipole (mu_Mol)
// MATLAB code uses 25°
// -----------------------------------------------------------
static Vec3 get_intrinsic_dipole()
{
    double angle = 25.0 * M_PI / 180.0;
    return {0.0, std::sin(angle), -std::cos(angle)};
}

// -----------------------------------------------------------
// Define intrinsic molecular Raman tensor (3×3)
// MATLAB code scales it by 5, rotated by alpha_angle = 34°
// -----------------------------------------------------------
static void get_intrinsic_raman(double A[3][3])
{
    A[0][0] = 0.05 * 5;
    A[0][1] = 0.0;
    A[0][2] = 0.0;

    A[1][0] = 0.0;
    A[1][1] = 0.20 * 5;
    A[1][2] = 0.0;

    A[2][0] = 0.0;
    A[2][1] = 0.0;
    A[2][2] = 1.0 * 5;
}

// -----------------------------------------------------------
// Rotate Raman tensor by alpha_angle = 34° (same as MATLAB)
// -----------------------------------------------------------
static void rotate_raman_to_molecule(double A[3][3])
{
    double angle = 34.0 * M_PI / 180.0;
    double c = std::cos(angle);
    double s = std::sin(angle);

    // Rotation matrix around x
    double R[3][3] = {
        {1, 0, 0},
        {0, c, -s},
        {0, s,  c}
    };

    double tmp[3][3];
    double Arot[3][3];

    // tmp = R * A
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            tmp[i][j] = R[i][0]*A[0][j] + R[i][1]*A[1][j] + R[i][2]*A[2][j];
        }
    }

    // Arot = tmp * R^T
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            Arot[i][j] = tmp[i][0]*R[j][0] + tmp[i][1]*R[j][1] + tmp[i][2]*R[j][2];
        }
    }

    // copy back
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            A[i][j] = Arot[i][j];
}

// -----------------------------------------------------------
// Rotate dipole into simulation frame: μ_sim = [X,Y,Z]^T * μ
// -----------------------------------------------------------
static Vec3 rotate_dipole(const Vec3& mu, const LocalFrame& f)
{
    return {
        f.X_axis.x*mu.x + f.Y_axis.x*mu.y + f.Z_axis.x*mu.z,
        f.X_axis.y*mu.x + f.Y_axis.y*mu.y + f.Z_axis.y*mu.z,
        f.X_axis.z*mu.x + f.Y_axis.z*mu.y + f.Z_axis.z*mu.z
    };
}

// -----------------------------------------------------------
// Rotate Raman tensor: α_sim = T * α * T^T
// where T is [X_axis; Y_axis; Z_axis]
// -----------------------------------------------------------
static void rotate_raman_to_sim(
    const double A_mol[3][3], const LocalFrame& f,
    double A_out[3][3])
{
    double T[3][3] = {
        {f.X_axis.x, f.Y_axis.x, f.Z_axis.x},
        {f.X_axis.y, f.Y_axis.y, f.Z_axis.y},
        {f.X_axis.z, f.Y_axis.z, f.Z_axis.z}
    };

    double tmp[3][3];

    // tmp = T * A_mol
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            tmp[i][j] = T[i][0]*A_mol[0][j] + T[i][1]*A_mol[1][j] + T[i][2]*A_mol[2][j];
        }
    }

    // A_out = tmp * T^T
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            A_out[i][j] = tmp[i][0]*T[j][0] + tmp[i][1]*T[j][1] + tmp[i][2]*T[j][2];
        }
    }
}

// -----------------------------------------------------------
// Reduce Raman tensor (XX, YY, ZZ, XY, YZ, XZ)
// -----------------------------------------------------------
static void reduce_raman(double A[3][3], double out[6])
{
    out[0] = A[0][0];  // XX
    out[1] = A[1][1];  // YY
    out[2] = A[2][2];  // ZZ
    out[3] = A[0][1];  // XY
    out[4] = A[1][2];  // YZ
    out[5] = A[0][2];  // XZ
}

// -----------------------------------------------------------
// PUBLIC API
// -----------------------------------------------------------
std::vector<AmideIProps> get_amideI_properties(
    const std::vector<LocalFrame>& frames)
{
    std::vector<AmideIProps> properties;
    properties.reserve(frames.size());

    // intrinsic molecular properties
    Vec3 mu_mol = get_intrinsic_dipole();

    double A_mol[3][3];
    get_intrinsic_raman(A_mol);
    rotate_raman_to_molecule(A_mol);

    // rotate for each residue
    for(const auto& f : frames)
    {
        Vec3 mu_sim = rotate_dipole(mu_mol, f);

        double A_sim[3][3];
        rotate_raman_to_sim(A_mol, f, A_sim);

        double A_red[6];
        reduce_raman(A_sim, A_red);

        AmideIProps p;
        p.dipole_sim = mu_sim;

        // full matrix
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                p.alpha_matrix[i][j] = A_sim[i][j];

        // MATLAB reshape(A_sim, 1, 9) column-major
        int k = 0;
        for(int col=0; col<3; col++)
            for(int row=0; row<3; row++)
                p.alpha_vectorized[k++] = A_sim[row][col];

        // reduced
        for(int i=0;i<6;i++)
            p.alpha_reduced[i] = A_red[i];

        properties.push_back(p);
    }

    return properties;
}
