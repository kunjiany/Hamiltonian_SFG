#include "Extract_Amide_Coordinates.hpp"
#include <vector>
#include <string>

std::vector<AmideIEntry> Extract_Amide_Coordinates(const std::vector<Atom>& atoms)
{
    std::vector<AmideIEntry> amides;

    int NumAtoms = atoms.size();

    //  Step 1: collect all C atoms in PDB order 
    std::vector<int> C_indices;
    C_indices.reserve(NumAtoms);
    for (int i = 0; i < NumAtoms; i++) {
        if (atoms[i].name == "C") {
            C_indices.push_back(i);
        }
    }

    //  Step 2: for each C, find the next O 
    std::vector<int> O_indices(C_indices.size(), -1);
    for (size_t k = 0; k < C_indices.size(); k++) {
        int idx = C_indices[k] + 1;
        while (idx < NumAtoms) {
            if (atoms[idx].name == "O") {
                O_indices[k] = idx;
                break;
            }
            idx++;
        }
    }

    // Step 3: from each O, find the next N or OXT 
    std::vector<int> N_indices(C_indices.size(), -1);
    for (size_t k = 0; k < C_indices.size(); k++) {
        int oi = O_indices[k];
        if (oi < 0) continue;

        int idx = oi + 1;
        while (idx < NumAtoms) {
            const std::string &nm = atoms[idx].name;
            if (nm == "N" || nm == "OXT") {
                N_indices[k] = idx;
                break;
            }
            idx++;
        }
    }

    // Step 4: assemble valid C/O/N triplets 
    for (size_t k = 0; k < C_indices.size(); k++) {
        if (O_indices[k] < 0 || N_indices[k] < 0) {
            continue;
        }

        amides.push_back({
            atoms[C_indices[k]],
            atoms[O_indices[k]],
            atoms[N_indices[k]]
        });
    }

    return amides;
}
