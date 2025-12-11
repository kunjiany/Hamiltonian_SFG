#include "hamiltonian_equiv_matlab.hpp"
#include "rotate_properties.hpp"
#include "compute_dipole_coupling.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <lapacke.h>


static bool diagonalize(std::vector<double>& H, int N, std::vector<double>& evals)
{
    evals.resize(N);
    int info = LAPACKE_dsyev(
        LAPACK_ROW_MAJOR, 'V', 'U',
        N, H.data(), N, evals.data()
    );
    return (info == 0);
}


HamiltonianEquivResult Hamiltonian_equiv_matlab(
    const std::vector<AmideIGeo>&   geo,
    const std::vector<AmideIProps>& props,
    const std::vector<AmideIFreq>&  freqs,
    double tilt_deg,
    double twist_deg)
{
    HamiltonianEquivResult out;

    int N = static_cast<int>(props.size());
    if (N == 0) return out;
    out.N = N;

    
    auto rotated = rotate_properties(props, tilt_deg, twist_deg);

    out.mu_rot.resize(N);
    out.alpha_rot.resize(N);

    for (int i = 0; i < N; ++i) {
        // μ (3 components)
        out.mu_rot[i] = rotated[i].dipole_rot;  // Vec3

        // α (9 components in column-major order)
        for (int k = 0; k < 9; ++k) {
            out.alpha_rot[i][k] = rotated[i].alpha_rot[k];
        }
    }
    //build hamiltonian
    std::vector<double> H(N * N, 0.0);

    // Diagonal: site frequencies
    for (int i = 0; i < N; ++i) {
        H[i * N + i] = freqs[i].freq;
    }

    // Dipole–dipole coupling prefactor (from MATLAB code)
    const double prefactor =
        5034.0 * std::pow((4.1058 / std::sqrt(1600.0)) * 3.144, 2);

    for (int i = 0; i < N; ++i) {
        const Vec3& Ri  = geo[i].vibration_center_coord;
        const Vec3& mui = out.mu_rot[i];

        for (int j = i + 1; j < N; ++j) {
            const Vec3& Rj  = geo[j].vibration_center_coord;
            const Vec3& muj = out.mu_rot[j];

            Vec3 Rij { Ri.x - Rj.x, Ri.y - Rj.y, Ri.z - Rj.z };

            double Jij = compute_dipole_coupling(
                mui, muj, Rij,
                1000.0,   // dielectric epsilon, same as before
                false,    // no screening flag
                prefactor
            );

            H[i * N + j] = Jij;
            H[j * N + i] = Jij;
        }
    }


    // Diagonalize H → eigenvalues + eigenvectors
    // H is overwritten with eigenvectors (columns) by LAPACKE_dsyev
    std::vector<double> evals;
    if (!diagonalize(H, N, evals)) {
        std::cerr << "[Hamiltonian_equiv_matlab] LAPACK dsyev failed\n";
        return out;
    }


    // 4. Sort eigenvalues and eigenvectors ascending (like MATLAB sort)
    std::vector<int> idx(N);
    for (int i = 0; i < N; ++i) idx[i] = i;

    std::sort(idx.begin(), idx.end(),
              [&](int a, int b) { return evals[a] < evals[b]; });

    out.Sort_Ex_Freq.resize(N);
    out.Sort_V.assign(N * N, 0.0);

    for (int newcol = 0; newcol < N; ++newcol) {
        int oldcol = idx[newcol];
        out.Sort_Ex_Freq[newcol] = evals[oldcol];

        for (int r = 0; r < N; ++r) {
            // H(r,oldcol) in row-major → H[r*N + oldcol]
            out.Sort_V[r * N + newcol] = H[r * N + oldcol];
        }
    }


    // Compute exciton μ_ex and α_ex:
    // μ_ex(k)   = Σ_i V(i,k) * μ_i
    // α_ex(k,:) = Σ_i V(i,k) * α_i(:)
    out.mu_ex.resize(N);
    out.alpha_ex.resize(N);

    for (int k = 0; k < N; ++k) {
        double mx = 0.0, my = 0.0, mz = 0.0;

        for (int i = 0; i < N; ++i) {
            double vik = out.Sort_V[i * N + k];  // eigenvector component

            mx += vik * out.mu_rot[i].x;
            my += vik * out.mu_rot[i].y;
            mz += vik * out.mu_rot[i].z;
        }

        out.mu_ex[k] = Vec3{ mx, my, mz };
    }

    for (int k = 0; k < N; ++k) {
        std::array<double, 9> acc{};
        for (int t = 0; t < 9; ++t) acc[t] = 0.0;

        for (int i = 0; i < N; ++i) {
            double vik = out.Sort_V[i * N + k];

            for (int t = 0; t < 9; ++t) {
                acc[t] += vik * out.alpha_rot[i][t];
            }
        }

        out.alpha_ex[k] = acc;
    }

    return out;
}
