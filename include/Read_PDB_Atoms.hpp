#ifndef READ_PDB_ATOMS_HPP
#define READ_PDB_ATOMS_HPP

#include <vector>
#include <string>

struct Atom {
    std::string name;   // AtomName (C, O, N, CA, etc.)
    int serial;         // NEW: PDB Atom Serial Number
    int resID;          // Residue Number
    double x, y, z;     // Coordinates
};

std::vector<Atom> Read_PDB_Atoms(const std::string &pdbFile);

#endif
