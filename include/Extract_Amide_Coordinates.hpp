#ifndef EXTRACT_AMIDE_COORDINATES_HPP
#define EXTRACT_AMIDE_COORDINATES_HPP

#include <vector>
#include "Read_PDB_Atoms.hpp"

struct AmideIEntry {
    Atom C, O, N;
};

std::vector<AmideIEntry> Extract_Amide_Coordinates(const std::vector<Atom>& atoms);

#endif

