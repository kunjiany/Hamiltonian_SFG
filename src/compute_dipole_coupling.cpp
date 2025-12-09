#include "compute_dipole_coupling.hpp"
#include <cmath>

double compute_dipole_coupling(
    const Vec3 &mu_i,
    const Vec3 &mu_j,
    const Vec3 &Ri_minus_Rj,
    double distance_cutoff,
    bool use_cutoff,
    double prefactor
)
{
    double R2 = Ri_minus_Rj.x * Ri_minus_Rj.x
              + Ri_minus_Rj.y * Ri_minus_Rj.y
              + Ri_minus_Rj.z * Ri_minus_Rj.z;

    if (R2 == 0.0) {
        // diagonal or same-site → no coupling
        return 0.0;
    }

    double R = std::sqrt(R2);

    // cutoff (MATLAB sets Beta=0 outside cutoff)
    if (use_cutoff && R > distance_cutoff)
        return 0.0;

    double R3 = R2 * R;
    double R5 = R3 * R2;

    // dot products
    double muIdotJ = mu_i.x * mu_j.x + mu_i.y * mu_j.y + mu_i.z * mu_j.z;

    double RdotI = Ri_minus_Rj.x * mu_i.x
                 + Ri_minus_Rj.y * mu_i.y
                 + Ri_minus_Rj.z * mu_i.z;

    double RdotJ = Ri_minus_Rj.x * mu_j.x
                 + Ri_minus_Rj.y * mu_j.y
                 + Ri_minus_Rj.z * mu_j.z;

    // MATLAB:
    // Beta = prefactor * (muIdotJ/R^3 - 3*(R·muI)(R·muJ)/R^5)
    double beta = prefactor * (muIdotJ / R3 - 3.0 * (RdotI * RdotJ) / R5);

    return beta;
}
