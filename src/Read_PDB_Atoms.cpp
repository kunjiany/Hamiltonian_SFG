#include "Read_PDB_Atoms.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>

std::vector<Atom> Read_PDB_Atoms(const std::string &pdbFile)
{
    std::ifstream fin(pdbFile);
    if(!fin){
        std::cerr << "ERROR: Cannot open PDB file: " << pdbFile << "\n";
        exit(1);
    }

    std::vector<Atom> atoms;
    std::string line;

    while(std::getline(fin, line))
    {
        if (line.size() < 54) continue;

        // Accept both ATOM and HETATM
        if (line.rfind("ATOM", 0) != 0 && line.rfind("HETATM", 0) != 0)
            continue;

        Atom a;

        // ------------------------
        // PDB atom serial (columns 7–11)
        // ------------------------
        std::string serial_str = line.substr(6, 5);
        serial_str.erase(std::remove(serial_str.begin(), serial_str.end(), ' '), serial_str.end());
        a.serial = std::stoi(serial_str);

        // Atom name (columns 12–16)
        a.name = line.substr(12, 4);
        a.name.erase(std::remove(a.name.begin(), a.name.end(), ' '), a.name.end());

        // Residue ID (columns 22–26)
        std::string resnum = line.substr(22, 4);
        resnum.erase(std::remove(resnum.begin(), resnum.end(), ' '), resnum.end());
        a.resID = std::stoi(resnum);

        // Coordinates (columns 30–54)
        a.x = std::stod(line.substr(30, 8));
        a.y = std::stod(line.substr(38, 8));
        a.z = std::stod(line.substr(46, 8));

        atoms.push_back(a);
    }

    return atoms;
}
