#ifndef FRESNEL_CALCULATION_HPP
#define FRESNEL_CALCULATION_HPP

#include "Fresnel_Read.hpp"

// Fresnel prefactors used in SFG scoring, matching MATLAB Fresnel_Prism.m
// at the polymer–water interface (PoWa).
struct FresnelFactors {
    double F_ssp_yyz = 0.0;  // SSP, χ_yyz
    double F_ppp_xxz = 0.0;  // PPP, χ_xxz
    double F_ppp_xzx = 0.0;  // PPP, χ_xzx
    double F_ppp_zxx = 0.0;  // PPP, χ_zxx
    double F_ppp_zzz = 0.0;  // PPP, χ_zzz
};

// Compute Fresnel prefactors exactly as in MATLAB Fresnel_Prism.m
// for the air–CaF2–polymer–water system, using only the polymer–water
// local field factors (PoWa) that appear in the scoring code.
FresnelFactors Fresnel_Calculation(const FresnelParams& fp);

#endif // FRESNEL_CALCULATION_HPP


