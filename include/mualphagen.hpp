#ifndef MUALPHAGEN_HPP
#define MUALPHAGEN_HPP

#include <vector>

// MATLAB equivalent "Trans_Loc" : StatesNum × StatesNum × dim
struct LocOp {
    int StatesNum = 0;     // number of states (N)
    int dim = 0;           // operator dimension (3 for μ, 9 for α)
    std::vector<double> data; // flattened array: (i*StatesNum + j)*dim + c
};

// C++ version of OneDSFG_MuAlphaGen
LocOp OneDSFG_MuAlphaGen(
    int StatesNum,
    const std::vector<std::vector<double>>& T_local
);

#endif
