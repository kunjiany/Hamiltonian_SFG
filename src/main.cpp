#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <filesystem>

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
    for (double twist_deg : twist_vec)
    {
        for (double tilt_deg : tilt_vec)
        {
            std::cout << "\n>>> tilt = " << tilt_deg
                      << "°, twist = " << twist_deg << "°\n";

            // (A) Exciton Hamiltonian (MATLAB-equivalent)
            HamiltonianEquivResult H =
                Hamiltonian_equiv_matlab(geo, props, freqs, tilt_deg, twist_deg);


            // (B) χ² (rotated to lab frame)
            R3Matrix R = Rdb.get_R(twist_deg, tilt_deg);
            // Print the first row of the rotation matrix for debugging
            std::cout << "R(tilt=" << tilt_deg << ", twist=" << twist_deg << ") first row: "
                    << R.m[0][0] << "  " << R.m[0][1] << "  " << R.m[0][2] << "\n" 
                    << R.m[1][0] << "  " << R.m[1][1] << "  " << R.m[1][2] << "\n" 
                    << R.m[2][0] << "  " << R.m[2][1] << "  " << R.m[2][2] << "\n" <<std::endl;


            Chi2Result chi = compute_chi2_matlab(H, R);

        // =========================================================
        // DEBUG: Print chi_mol and chi_lab for this tilt/twist
        // =========================================================
        {
            std::string tag =
                "tilt" + std::to_string((int)std::round(tilt_deg)) +
                "_twist" + std::to_string((int)std::round(twist_deg));

            std::ofstream fmol("debug1/chi_mol_" + tag + ".txt");
            std::ofstream flab("debug1/chi_lab_" + tag + ".txt");

            if (!fmol || !flab) {
                std::cerr << "ERROR: cannot open debug1 directory or files.\n";
            } else {
                fmol << std::setprecision(10);
                flab << std::setprecision(10);

                for (int k = 0; k < chi.N; k++) {
                    fmol << "# exciton " << k << "\n";
                    flab << "# exciton " << k << "\n";

                    // χ_mol (27 numbers)
                    for (int i = 0; i < 27; i++)
                        fmol << chi.chi_mol[k][i] << (i == 26 ? '\n' : ' ');

                    // χ_lab (27 numbers)
                    for (int i = 0; i < 27; i++)
                        flab << chi.chi_lab[k][i] << (i == 26 ? '\n' : ' ');
                }
            }
        }

            // ==========================================================

            // (C) Compute raw SFG (NO Fresnel)
            SpectrumResult spec =
                compute_SFG_spectra(H, chi, in.width, freq_grid);

            // ----------------------------------------------------
            // (D) Write MATLAB-style output: freq, SSP, PPP
            // ----------------------------------------------------
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

            std::cout << "Wrote: " << fname << "\n";
        }
    }

    std::cout << "\n=== Completed full SFG pipeline ===\n";
    return 0;
}
