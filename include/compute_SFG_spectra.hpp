#ifndef COMPUTE_SFG_SPECTRA_HPP
#define COMPUTE_SFG_SPECTRA_HPP

#include <vector>
#include "hamiltonian_equiv_matlab.hpp"
#include "chi2_matlab.hpp"


struct SpectrumResult {
    std::vector<double> freq;
    std::vector<double> I_ssp;   // yyz = index 22
    std::vector<double> I_ppp;   // zzz = index 26
};


SpectrumResult compute_SFG_spectra(
    const HamiltonianEquivResult& H,
    const Chi2Result& chi,
    double width,
    const std::vector<double>& freq_grid
);

#endif

