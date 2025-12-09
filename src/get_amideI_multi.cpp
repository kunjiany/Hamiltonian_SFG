#include "get_amideI_multi.hpp"
#include <iostream>
#include <stdexcept>

// ================================================================
// Utility: slice a vector (MATLAB helixa:helixb)
// ================================================================
template<typename T>
std::vector<T> slice_vector(const std::vector<T>& v, int a, int b)
{
    if (a < 1 || b > (int)v.size() || a > b)
        throw std::runtime_error("Invalid helix range");

    return std::vector<T>(v.begin() + (a - 1), v.begin() + b);
}

// ================================================================
// Utility: replicate ModeNum → ModeNum * layer
// (MATLAB's layer logic)
// ================================================================
template<typename T>
std::vector<T> replicate_layer(const std::vector<T>& data, int layer)
{
    std::vector<T> out;
    out.reserve(data.size() * layer);

    for (int L = 0; L < layer; L++) {
        out.insert(out.end(), data.begin(), data.end());
    }
    return out;
}



// ============================================================================
// MAIN IMPLEMENTATION — EXACT SEQUENCE OF MATLAB GetAmideIMulti.m
// ============================================================================
AmideIMultiOutput Get_AmideI_Multi(
    double center_freq,
    int helix_a,
    int helix_b,
    int layer,
    const std::string& pdbFile)
{
    AmideIMultiOutput out;

    // --------------------------------------------------------------------
    // 1. Read PDB and extract atoms
    // --------------------------------------------------------------------
    auto atoms = Read_PDB_Atoms(pdbFile);

    // --------------------------------------------------------------------
    // 2. Extract C/O/N triplets (MATLAB’s AmideIAtomSerNo logic)
    // --------------------------------------------------------------------
    auto amide_all = Extract_Amide_Coordinates(atoms);

    // --------------------------------------------------------------------
    // 3. Slice helix segment (helixa:helixb)
    // --------------------------------------------------------------------
    auto amide_seg = slice_vector(amide_all, helix_a, helix_b);
    int ModeNum = amide_seg.size(); //number of mode is the same as the number of amide

    // --------------------------------------------------------------------
    // 4. Build geometry: CO, CN, center
    // --------------------------------------------------------------------
    auto geo = get_amideI_geometry(amide_seg);

    // --------------------------------------------------------------------
    // 5. Local frame X, Y, Z (simulation frame)
    // --------------------------------------------------------------------
    auto frames = get_local_frame(geo);

    // --------------------------------------------------------------------
    // 6. Dipole & Raman properties
    // --------------------------------------------------------------------
    auto props = get_amideI_properties(frames);

    // --------------------------------------------------------------------
    // Layer replication (ModeNum → ModeNum * layer)
    // --------------------------------------------------------------------
    auto geo_L     = replicate_layer(geo, layer);
    auto frames_L  = replicate_layer(frames, layer);
    auto props_L   = replicate_layer(props, layer);
    auto amide_L   = replicate_layer(amide_seg, layer);
    int total_modes = ModeNum * layer;

    // --------------------------------------------------------------------
    // 7. Frequencies
    // --------------------------------------------------------------------
    out.freq.resize(total_modes, center_freq);

    // --------------------------------------------------------------------
    // 8. Anharmonicity = 12 for all
    // --------------------------------------------------------------------
    out.anharm.resize(total_modes, 12.0);

    // --------------------------------------------------------------------
    // 9. Fill output arrays
    // --------------------------------------------------------------------
    out.center.resize(total_modes);
    out.mu_orig.resize(total_modes);
    out.alpha_matrix.resize(total_modes);
    out.alpha_reduced.resize(total_modes);
    out.alpha_vectorized.resize(total_modes);

    out.AtomSerNo.resize(total_modes);
    out.xyz.resize(total_modes);

    out.AtomC.resize(total_modes);
    out.AtomO.resize(total_modes);
    out.AtomN.resize(total_modes);

    for (int i = 0; i < total_modes; i++) {

        // -----------------------------
        // center
        // -----------------------------
        out.center[i] = geo_L[i].vibration_center_coord;

        // -----------------------------
        // dipole
        // -----------------------------
        out.mu_orig[i] = props_L[i].dipole_sim;

        // -----------------------------
        // Raman tensors (must match MATLAB)
        // -----------------------------
        for(int r=0;r<3;r++)
            for(int c=0;c<3;c++)
                out.alpha_matrix[i][r*3 + c] = props_L[i].alpha_matrix[r][c];

        {
            int k = 0;
            for (int col = 0; col < 3; col++) {
                for (int row = 0; row < 3; row++) {
                    out.alpha_vectorized[i][k++] = props_L[i].alpha_matrix[row][col];
                }
            }
        }

        for (int k = 0; k < 6; k++)
            out.alpha_reduced[i][k] = props_L[i].alpha_reduced[k];

        // -----------------------------
        // Raw atom serial numbers
        // -----------------------------
        out.AtomSerNo[i] = {
            amide_L[i].C.serial,
            amide_L[i].O.serial,
            amide_L[i].N.serial
        };

        // -----------------------------
        // Raw xyz coordinates
        // -----------------------------
        out.AtomC[i] = { amide_L[i].C.x, amide_L[i].C.y, amide_L[i].C.z };
        out.AtomO[i] = { amide_L[i].O.x, amide_L[i].O.y, amide_L[i].O.z };
        out.AtomN[i] = { amide_L[i].N.x, amide_L[i].N.y, amide_L[i].N.z };

        out.xyz[i] = { out.AtomC[i], out.AtomO[i], out.AtomN[i] };
    }

    return out;
}
