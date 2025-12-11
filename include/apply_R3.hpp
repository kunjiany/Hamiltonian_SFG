#ifndef APPLY_R3_HPP
#define APPLY_R3_HPP

#include <array>
#include "load_R3ZXZ1.hpp"   // for R3Matrix

// Apply R3 (27×27) to a single chi vector (27×1)
std::array<double,27> apply_R3_single(
    const R3Matrix& R,
    const std::array<double,27>& chi_in
);

#endif
