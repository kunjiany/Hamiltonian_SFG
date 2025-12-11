#include "get_amideI_multi.hpp"
#include <iostream>
#include <stdexcept>


template<typename T>
std::vector<T> slice_vector(const std::vector<T>& v, int a, int b)
{
    if (a < 1 || b > (int)v.size() || a > b)
        throw std::runtime_error("Invalid helix range");

    return std::vector<T>(v.begin() + (a - 1), v.begin() + b);
}

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



AmideIMultiOutput Get_AmideI_Multi(
    double center_freq,
    int helix_a,
    int helix_b,
    int layer,
    const std::string& pdbFile)
{
    AmideIMultiOutput out;
    auto atoms = Read_PDB_Atoms(pdbFile);

    auto amide_all = Extract_Amide_Coordinates(atoms);

    auto amide_seg = slice_vector(amide_all, helix_a, helix_b);
    int ModeNum = amide_seg.size(); //number of mode is the same as the number of amide

    auto geo = get_amideI_geometry(amide_seg);

    auto frames = get_local_frame(geo);

    auto props = get_amideI_properties(frames);
    
    auto geo_L     = replicate_layer(geo, layer);
    auto frames_L  = replicate_layer(frames, layer);
    auto props_L   = replicate_layer(props, layer);
    auto amide_L   = replicate_layer(amide_seg, layer);
    int total_modes = ModeNum * layer;

    out.freq.resize(total_modes, center_freq);


    out.anharm.resize(total_modes, 12.0);

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

        out.center[i] = geo_L[i].vibration_center_coord;
        out.mu_orig[i] = props_L[i].dipole_sim;
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

        out.AtomSerNo[i] = {
            amide_L[i].C.serial,
            amide_L[i].O.serial,
            amide_L[i].N.serial
        };

        out.AtomC[i] = { amide_L[i].C.x, amide_L[i].C.y, amide_L[i].C.z };
        out.AtomO[i] = { amide_L[i].O.x, amide_L[i].O.y, amide_L[i].O.z };
        out.AtomN[i] = { amide_L[i].N.x, amide_L[i].N.y, amide_L[i].N.z };

        out.xyz[i] = { out.AtomC[i], out.AtomO[i], out.AtomN[i] };
    }

    return out;
}
