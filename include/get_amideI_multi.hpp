#ifndef GET_AMIDEI_MULTI_HPP
#define GET_AMIDEI_MULTI_HPP

#include <vector>
#include <string>
#include <array>

#include "helper_vec3.hpp"
#include "Read_PDB_Atoms.hpp"
#include "Extract_Amide_Coordinates.hpp"
#include "get_amideI_geometry.hpp"
#include "get_local_frame.hpp"
#include "get_amideI_properties.hpp"


struct AmideIMultiOutput {

    // Geometry: center of vibration
    std::vector<Vec3> center;  // size = ModeNum * layer

    // Frequency & anharmonicity
    std::vector<double> freq;      // center_freq (or isotopes)
    std::vector<double> anharm;    // always 12.0

    // Dipole (mu_orig)
    std::vector<Vec3> mu_orig;

    // Raman tensors
    std::vector<std::array<double,9>> alpha_matrix;     // full 3×3
    std::vector<std::array<double,9>> alpha_vectorized; // flatten 3×3
    std::vector<std::array<double,6>> alpha_reduced;    // 6 unique

    // Raw atom indices & coordinates
    std::vector<std::array<int,3>> AtomSerNo;          // [C, O, N] serials
    std::vector<std::array<Vec3,3>> xyz;               // xyz[i][0]=C,1=O,2=N

    std::vector<Vec3> AtomC;
    std::vector<Vec3> AtomO;
    std::vector<Vec3> AtomN;
};


AmideIMultiOutput Get_AmideI_Multi(
    double center_freq,
    int helix_a,
    int helix_b,
    int layer,
    const std::string& pdbFile
    // TODO: later isotopes & label frequencies
);

#endif
