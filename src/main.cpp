#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <filesystem>
#include <omp.h>

#include "Read_Input.hpp"
#include "Fresnel_Read.hpp"          // only for compatibility, not used
#include "Fresnel_Calculation.hpp"   // optional
#include "generate_angles.hpp"

#include "get_amideI_multi.hpp"
#include "hamiltonian_equiv_matlab.hpp"
#include "chi2_matlab.hpp"
#include "load_R3ZXZ1.hpp"
#include "compute_SFG_spectra.hpp"

// ================================================================
// MAIN
// ================================================================
int main()
{
    std::cout << "=== FULL Amide-I + SFG Pipeline (MATLAB Equivalent) ===\n";
    omp_set_dynamic(0);
    // ------------------------------------------------------------
    // 1. Read input.txt (same as old main)
    // ------------------------------------------------------------
    InputParams in = Read_Input("./input/input.txt");

    std::cout << "PDB: " << in.pdbFile << "\n";
    std::cout << "Center Freq: " << in.centerFreq << "\n";
    std::cout << "Tilt range:  " << in.tilt_start
              << " to " << in.tilt_end
              << " (" << in.tilt_points << " points)\n";
    std::cout << "Twist range: " << in.twist_start
              << " to " << in.twist_end
              << " (" << in.twist_points << " points)\n";

    // ------------------------------------------------------------
    // Create output directory
    // ------------------------------------------------------------
    std::filesystem::create_directories(in.SpectraFolder);

    // ------------------------------------------------------------
    // Frequency grid
    // ------------------------------------------------------------
    std::vector<double> freq_grid;
    for (double w = in.spec_range_start; w <= in.spec_range_end; w += in.spec_range_step)
        freq_grid.push_back(w);

    // ------------------------------------------------------------
    // 2. Load site-level amide modes
    // ------------------------------------------------------------
    AmideIMultiOutput M =
        Get_AmideI_Multi(in.centerFreq, 1, 5, in.layer, in.pdbFile);

    int N = M.center.size();
    std::cout << "Total modes = " << N << "\n";

    // Prepare Hamiltonian input arrays
    std::vector<AmideIGeo>  geo(N);
    std::vector<AmideIProps> props(N);
    std::vector<AmideIFreq> freqs(N);

    for (int i = 0; i < N; i++)
    {
        geo[i].vibration_center_coord = M.center[i];
        props[i].dipole_sim           = M.mu_orig[i];

        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                props[i].alpha_matrix[r][c] = M.alpha_matrix[i][r*3 + c];

        freqs[i].freq = M.freq[i];
    }

    // ------------------------------------------------------------
    // 3. Load R3 ZXZ rotation database
    // ------------------------------------------------------------
    R3Database Rdb("./data/R3ZXZ1_database.h5");

    // Generate tilt/twist vectors (same as old MATLAB driver)
    auto tilt_vec  = Linspace(in.tilt_start,  in.tilt_end,  in.tilt_points);
    auto twist_vec = Linspace(in.twist_start, in.twist_end, in.twist_points);

    // ------------------------------------------------------------
    // 4. Orientation loop (same structure as old main)
    // ------------------------------------------------------------

    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int it = 0; it < (int)twist_vec.size(); it++)
    {
        for (int jt = 0; jt < (int)tilt_vec.size(); jt++)
        {
            double twist_deg = twist_vec[it];
            double tilt_deg  = tilt_vec[jt];

            // Each thread prints independently → optional mutex if needed
            #pragma omp critical
            {
                std::cout << "\n>>> tilt = " << tilt_deg
                        << "°, twist = " << twist_deg << "°\n";
            }

            // =============================================
            // A) Exciton Hamiltonian
            // =============================================
            HamiltonianEquivResult H =
                Hamiltonian_equiv_matlab(geo, props, freqs, tilt_deg, twist_deg);

            // =============================================
            // B) χ² using R3 lookup
            // =============================================
            R3Matrix R = Rdb.get_R(twist_deg, tilt_deg);

            Chi2Result chi = compute_chi2_matlab(H, R);

            // =============================================
            // DEBUG OUTPUT (thread safe)
            // =============================================
            {
                std::string tag =
                    "tilt" + std::to_string((int)std::round(tilt_deg)) +
                    "_twist" + std::to_string((int)std::round(twist_deg));

                {
                    std::ofstream fmol("debug1/chi_mol_" + tag + ".txt");
                    std::ofstream flab("debug1/chi_lab_" + tag + ".txt");

                    if (fmol && flab)
                    {
                        fmol << std::setprecision(10);
                        flab << std::setprecision(10);

                        for (int k = 0; k < chi.N; k++) {
                            fmol << "# exciton " << k << "\n";
                            flab << "# exciton " << k << "\n";

                            for (int i = 0; i < 27; i++)
                                fmol << chi.chi_mol[k][i] << (i == 26 ? '\n' : ' ');

                            for (int i = 0; i < 27; i++)
                                flab << chi.chi_lab[k][i] << (i == 26 ? '\n' : ' ');
                        }
                    }
                }
            }

            // =============================================
            // C) SFG calculation
            // =============================================
            SpectrumResult spec =
                compute_SFG_spectra(H, chi, in.width, freq_grid);

            // =============================================
            // D) Write output spectra (thread-safe)
            // =============================================
            {
                std::string fname =
                    in.SpectraFolder + "/" +
                    in.SpectraStorePrefix +
                    "_tilt"  + std::to_string((int)std::round(tilt_deg)) +
                    "_twist" + std::to_string((int)std::round(twist_deg)) +
                    ".txt";

                std::ofstream fout(fname);
                fout << "# freq   SSP(yyz)   PPP(zzz)\n";
                fout << std::setprecision(10);

                for (size_t i = 0; i < spec.freq.size(); i++)
                {
                    fout << spec.freq[i] << " "
                        << spec.I_ssp[i] << " "
                        << spec.I_ppp[i] << "\n";
                }

                #pragma omp critical
                {
                    std::cout << "Wrote: " << fname << "\n";
                }
            }
        }
    }
    std::cout << "\n=== Completed full SFG pipeline ===\n";
    return 0;
}
