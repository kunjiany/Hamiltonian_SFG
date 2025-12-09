#include "mualphagen.hpp"
#include <vector>
#include <iostream>

LocOp OneDSFG_MuAlphaGen(
    int StatesNum,
    const std::vector<std::vector<double>>& T_local)
{
    LocOp op;
    op.StatesNum = StatesNum;

    if (StatesNum <= 0) {
        std::cerr << "[OneDSFG_MuAlphaGen] ERROR: StatesNum <= 0\n";
        return op;
    }

    if ((int)T_local.size() != StatesNum) {
        std::cerr << "[OneDSFG_MuAlphaGen] ERROR: size mismatch\n";
        return op;
    }

    // dim = 3 for μ, 9 for α
    op.dim = T_local[0].size();

    const int N = StatesNum;
    const int dim = op.dim;

    // Allocate N×N×dim, zero-initialized
    op.data.assign(N * N * dim, 0.0);

    // MATLAB:
    // A_loc(1,i,:) = T_local(i,:)   % others zero
    //
    // C++ 0-based:
    // A_loc(0,i,:) = T_local[i]

    for (int i = 0; i < N; i++) {
        int base = (0 * N + i) * dim;   // row=0, col=i
        for (int c = 0; c < dim; c++) {
            op.data[base + c] = T_local[i][c];
        }
    }

    return op;
}
