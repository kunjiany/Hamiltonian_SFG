#include "compute_SFG_spectra.hpp"
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>

using cplx = std::complex<double>;

// Complex Lorentzian: 1 / (Δω + i Γ)
static inline cplx lorentz(double dw, double width) {
    return cplx(1.0, 0.0) / cplx(dw, width);
}

SpectrumResult compute_SFG_spectra(
    const HamiltonianEquivResult& H,
    const Chi2Result& chi,
    double width,
    const std::vector<double>& freq_grid)
{
    SpectrumResult out;
    out.freq = freq_grid;

    int N = H.N;
    if (N <= 0) {
        std::cerr << "[SFG] Error: N = 0\n";
        return out;
    }

    if ((int)chi.chi_lab.size() != N) {
        std::cerr << "[SFG] Error: chi_lab size mismatch\n";
        return out;
    }

    size_t nf = freq_grid.size();
    out.I_ssp.assign(nf, 0.0);
    out.I_ppp.assign(nf, 0.0);

    // SSP  MATLAB column 23   C++ index 22
    // PPP   MATLAB column 27   C++ index 26
    const int IDX_SSP = 22;
    const int IDX_PPP = 26;

    for (size_t i = 0; i < nf; i++)
    {
        double w = freq_grid[i];
        cplx sum_ssp(0, 0);
        cplx sum_ppp(0, 0);

        for (int k = 0; k < N; k++)
        {
            double wk = H.Sort_Ex_Freq[k];
            double dw = w - wk;

            double chi_ssp_k = chi.chi_lab[k][IDX_SSP];
            double chi_ppp_k = chi.chi_lab[k][IDX_PPP];

            cplx Lk = lorentz(dw, width);

            sum_ssp += cplx(chi_ssp_k, 0.0) * Lk;
            sum_ppp += cplx(chi_ppp_k, 0.0) * Lk;
        }

        // RAW intensities (MATLAB SFG_Calculation.m)
        out.I_ssp[i] = std::norm(sum_ssp);
        out.I_ppp[i] = std::norm(sum_ppp);
    }

    // Write MATLAB-format output
    {
        std::ofstream fs("cpp_sfg_ssp.txt");
        fs << std::setprecision(10);
        for (size_t i = 0; i < nf; i++)
            fs << out.freq[i] << " " << out.I_ssp[i] << "\n";
    }
    {
        std::ofstream fp("cpp_sfg_ppp.txt");
        fp << std::setprecision(10);
        for (size_t i = 0; i < nf; i++)
            fp << out.freq[i] << " " << out.I_ppp[i] << "\n";
    }

    std::cout << "[SFG] Wrote cpp_sfg_ssp.txt and cpp_sfg_ppp.txt (RAW, no Fresnel)\n";
    return out;
}
