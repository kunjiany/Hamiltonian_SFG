#include "Fresnel_Calculation.hpp"
#include <cmath>
#include <complex>
#include <iostream>


//Not use now
using cd = std::complex<double>;

// -----------------------------------------------------------
// Helpers
// -----------------------------------------------------------
static inline cd cexp_i(double x) {
    return cd(std::cos(x), std::sin(x));
}

// Clamp to [-1, 1] to avoid asin() domain errors caused by
// small floating-point overshoots, just like MATLAB effectively does.
static inline double clamp_sin(double x) {
    if (x > 1.0)  return 1.0;
    if (x < -1.0) return -1.0;
    return x;
}

FresnelFactors Fresnel_Calculation(const FresnelParams& fp)
{
    FresnelFactors FF;

    // Short names mapping directly to MATLAB variables
    const double d         = fp.polymer_thickness_nm;   // nm
    const double A_Vis     = fp.A_vis_deg * M_PI / 180.0; // A_Vis (rad)
    const double A_Pump    = fp.A_ir_deg  * M_PI / 180.0; // A_Pump (rad)

    // Visible indices + wavelength (vi)
    const double n0vi      = fp.n0_vis;
    const double n1vi      = fp.n1_vis;
    const double n2vi      = fp.n2_vis;
    const double n3vi      = fp.n3_vis;
    const double Lambda_vi = fp.lambda_vis_nm;

    // IR indices + wavelength (in)
    const double n0in      = fp.n0_ir;
    const double n1in      = fp.n1_ir;
    const double n2in      = fp.n2_ir;
    const double n3in      = fp.n3_ir;
    const double Lambda_in = fp.lambda_ir_nm;

    // SFG wavelength + indices (su)
    const double Lambda_su = fp.lambda_sfg_nm;
    const double n0su      = fp.n0_sfg;
    const double n1su      = fp.n1_sfg;
    const double n2su      = fp.n2_sfg;
    const double n3su      = fp.n3_sfg;

    // ============================================================
    // 1) Visible beam (vi): angles via Snell's law
    //    Direct translation of MATLAB:
    //
    // Sigma0vi = pi/2 - A_Vis;
    // Sin_Sigma0vi = sin(Sigma0vi);
    // Sin_Sigma1vi = (n0vi * Sin_Sigma0vi)/n1vi;
    // Sigma1vi = asin(Sin_Sigma1vi);
    // ...
    // ============================================================
    double Sigma0vi    = M_PI/2.0 - A_Vis;
    double Sin_Sigma0vi = std::sin(Sigma0vi);
    double Cos_Sigma0vi = std::cos(Sigma0vi);

    double Sin_Sigma1vi = clamp_sin((n0vi * Sin_Sigma0vi) / n1vi);
    double Sigma1vi     = std::asin(Sin_Sigma1vi);
    double Cos_Sigma1vi = std::cos(Sigma1vi);

    double Phi1vi       = M_PI/2.0 - Sigma1vi;
    double Sin_Phi1vi   = std::sin(Phi1vi);
    double Cos_Phi1vi   = std::cos(Phi1vi);

    double Sin_Phi2vi = clamp_sin((n1vi * Sin_Phi1vi) / n2vi);
    double Phi2vi     = std::asin(Sin_Phi2vi);
    double Cos_Phi2vi = std::cos(Phi2vi);

    double Sin_Phi3vi = clamp_sin((n1vi * Sin_Phi1vi) / n3vi);
    double Phi3vi     = std::asin(Sin_Phi3vi);
    double Cos_Phi3vi = std::cos(Phi3vi);

    // Fresnel coefficients (vi), exactly as MATLAB
    double rp12vi = ( n2vi*Cos_Phi1vi - n1vi*Cos_Phi2vi )
                  / ( n2vi*Cos_Phi1vi + n1vi*Cos_Phi2vi );
    double tp12vi = ( 2.0*n1vi*Cos_Phi1vi )
                  / ( n2vi*Cos_Phi1vi + n1vi*Cos_Phi2vi );

    double rs12vi = ( n1vi*Cos_Phi1vi - n2vi*Cos_Phi2vi )
                  / ( n1vi*Cos_Phi1vi + n2vi*Cos_Phi2vi );
    double ts12vi = ( 2.0*n1vi*Cos_Phi1vi )
                  / ( n1vi*Cos_Phi1vi + n2vi*Cos_Phi2vi );

    double rp23vi = ( n2vi*Cos_Phi2vi - n3vi*Cos_Phi3vi )
                  / ( n3vi*Cos_Phi2vi + n2vi*Cos_Phi3vi );
    double tp23vi = ( 2.0*n2vi*Cos_Phi2vi )
                  / ( n3vi*Cos_Phi2vi + n2vi*Cos_Phi3vi );

    double rs23vi = ( n2vi*Cos_Phi2vi - n3vi*Cos_Phi3vi )
                  / ( n2vi*Cos_Phi2vi + n3vi*Cos_Phi3vi );
    double ts23vi = ( 2.0*n3vi*Cos_Phi3vi )
                  / ( n2vi*Cos_Phi2vi + n3vi*Cos_Phi3vi );

    // Thetavi = (2*pi/Lambda_vi) * n2vi * d * Cos_Phi2vi;
    double Thetavi = (2.0 * M_PI / Lambda_vi) * n2vi * d * Cos_Phi2vi;
    cd Evi = cexp_i(2.0 * Thetavi);

    // tpvi, tsvi: field transmitted from 0 to d for 3 layers
    cd tpvi = tp12vi / ( cd(1.0,0.0) + rp12vi*rp23vi * Evi );
    cd tsvi = ts12vi / ( cd(1.0,0.0) + rs12vi*rs23vi * Evi );

    // ============================================================
    // 2) IR beam (in)
    //    Direct translation of MATLAB block
    // ============================================================
    double Sigma0in    = M_PI/2.0 - A_Pump;
    double Sin_Sigma0in = std::sin(Sigma0in);
    double Cos_Sigma0in = std::cos(Sigma0in);

    double Sin_Sigma1in = clamp_sin((n0in * Sin_Sigma0in) / n1in);
    double Sigma1in     = std::asin(Sin_Sigma1in);
    double Cos_Sigma1in = std::cos(Sigma1in);

    double Phi1in       = M_PI/2.0 - Sigma1in;
    double Sin_Phi1in   = std::sin(Phi1in);
    double Cos_Phi1in   = std::cos(Phi1in);

    double Sin_Phi2in = clamp_sin((n1in * Sin_Phi1in) / n2in);
    double Phi2in     = std::asin(Sin_Phi2in);
    double Cos_Phi2in = std::cos(Phi2in);

    double Sin_Phi3in = clamp_sin((n1in * Sin_Phi1in) / n3in);
    double Phi3in     = std::asin(Sin_Phi3in);
    double Cos_Phi3in = std::cos(Phi3in);

    double rp12in = ( n2in*Cos_Phi1in - n1in*Cos_Phi2in )
                  / ( n2in*Cos_Phi1in + n1in*Cos_Phi2in );
    double tp12in = ( 2.0*n1in*Cos_Phi1in )
                  / ( n2in*Cos_Phi1in + n1in*Cos_Phi2in );

    double rs12in = ( n1in*Cos_Phi1in - n2in*Cos_Phi2in )
                  / ( n1in*Cos_Phi1in + n2in*Cos_Phi2in );
    double ts12in = ( 2.0*n1in*Cos_Phi1in )
                  / ( n1in*Cos_Phi1in + n2in*Cos_Phi2in );

    double rp23in = ( n2in*Cos_Phi2in - n3in*Cos_Phi3in )
                  / ( n3in*Cos_Phi2in + n2in*Cos_Phi3in );
    double tp23in = ( 2.0*n2in*Cos_Phi2in )
                  / ( n3in*Cos_Phi2in + n2in*Cos_Phi3in );

    double rs23in = ( n2in*Cos_Phi2in - n3in*Cos_Phi3in )
                  / ( n2in*Cos_Phi2in + n3in*Cos_Phi3in );
    double ts23in = ( 2.0*n3in*Cos_Phi3in )
                  / ( n2in*Cos_Phi2in + n3in*Cos_Phi3in );

    double Thetain = (2.0 * M_PI / Lambda_in) * n2in * d * Cos_Phi2in;
    cd Ein = cexp_i(2.0 * Thetain);

    cd tpin = tp12in / ( cd(1.0,0.0) + rp12in*rp23in * Ein );
    cd tsin = ts12in / ( cd(1.0,0.0) + rs12in*rs23in * Ein );

    // ============================================================
    // 3) SFG beam (su)
    //    Direct translation of the "sum generation beam" section
    // ============================================================
    double Sin_Phi1su =
        Lambda_su / n1su *
        ( n1vi*Sin_Phi1vi / Lambda_vi + n1in*Sin_Phi1in / Lambda_in );
    Sin_Phi1su = clamp_sin(Sin_Phi1su);
    double Phi1su     = std::asin(Sin_Phi1su);
    double Cos_Phi1su = std::cos(Phi1su);

    double Sin_Phi2su = clamp_sin( (n1su * Sin_Phi1su) / n2su );
    double Phi2su     = std::asin(Sin_Phi2su);
    double Cos_Phi2su = std::cos(Phi2su);

    double Sin_Phi3su = clamp_sin( (n1su * Sin_Phi1su) / n3su );
    double Phi3su     = std::asin(Sin_Phi3su);
    double Cos_Phi3su = std::cos(Phi3su);

    double Sigma1su    = Phi1su - M_PI/4.0;
    double Sin_Sigma1su = std::sin(Sigma1su);
    double Cos_Sigma1su = std::cos(Sigma1su);

    double Sin_Sigma0su = clamp_sin( (n1su * Sin_Sigma1su) / n0su );
    double Sigma0su     = std::asin(Sin_Sigma0su);
    double Cos_Sigma0su = std::cos(Sigma0su);

    // Fresnel (su)
    double rp12su = ( n2su*Cos_Phi1su - n1su*Cos_Phi2su )
                  / ( n2su*Cos_Phi1su + n1su*Cos_Phi2su );
    double tp12su = ( 2.0*n1su*Cos_Phi1su )
                  / ( n2su*Cos_Phi1su + n1su*Cos_Phi2su );

    double rs12su = ( n1su*Cos_Phi1su - n2su*Cos_Phi2su )
                  / ( n1su*Cos_Phi1su + n2su*Cos_Phi2su );
    double ts12su = ( 2.0*n1su*Cos_Phi1su )
                  / ( n1su*Cos_Phi1su + n2su*Cos_Phi2su );

    double rp23su = ( n2su*Cos_Phi2su - n3su*Cos_Phi3su )
                  / ( n3su*Cos_Phi2su + n2su*Cos_Phi3su );
    double tp23su = ( 2.0*n2su*Cos_Phi2su )
                  / ( n3su*Cos_Phi2su + n2su*Cos_Phi3su );

    double rs23su = ( n2su*Cos_Phi2su - n3su*Cos_Phi3su )
                  / ( n2su*Cos_Phi2su + n3su*Cos_Phi3su );
    double ts23su = ( 2.0*n3su*Cos_Phi3su )
                  / ( n2su*Cos_Phi2su + n3su*Cos_Phi3su );

    double Thetasu = (2.0 * M_PI / Lambda_su) * n2su * d * Cos_Phi2su;
    cd Esu = cexp_i(2.0 * Thetasu);

    cd tpsu = tp12su / ( cd(1.0,0.0) + rp12su*rp23su * Esu );
    cd tssu = ts12su / ( cd(1.0,0.0) + rs12su*rs23su * Esu );

    // ============================================================
    // 4) Polymer–Water interface (PoWa) local fields
    //    Only these are used in the scoring (MATLAB)
    // ============================================================
    double nmPoWasu = (n2su + n3su) / 2.0;
    double nmPoWavi = (n2vi + n3vi) / 2.0;
    double nmPoWain = (n2in + n3in) / 2.0;

    cd e_iThetasu = cexp_i(Thetasu);
    cd e_iThetavi = cexp_i(Thetavi);
    cd e_iThetain = cexp_i(Thetain);

    // SFG (su) local fields at PoWa
    cd LxxPoWasu = tpsu * e_iThetasu * (1.0 - rp23su) * (Cos_Phi2su / Cos_Phi1su);
    cd LyyPoWasu = tssu * e_iThetasu * (1.0 + rs23su);
    cd LzzPoWasu = tpsu * e_iThetasu * (1.0 + rp23su)
                 * (n1su * n2su) / (nmPoWasu * nmPoWasu);

    // Visible (vi) local fields at PoWa
    cd LxxPoWavi = tpvi * e_iThetavi * (1.0 - rp23vi) * (Cos_Phi2vi / Cos_Phi1vi);
    cd LyyPoWavi = tsvi * e_iThetavi * (1.0 + rs23vi);
    cd LzzPoWavi = tpvi * e_iThetavi * (1.0 + rp23vi)
                 * (n1vi * n2vi) / (nmPoWavi * nmPoWavi);

    // IR (in) local fields at PoWa
    cd LxxPoWain = tpin * e_iThetain * (1.0 - rp23in) * (Cos_Phi2in / Cos_Phi1in);
    cd LyyPoWain = tsin * e_iThetain * (1.0 + rs23in);
    cd LzzPoWain = tpin * e_iThetain * (1.0 + rp23in)
                 * (n1in * n2in) / (nmPoWain * nmPoWain);

    // ============================================================
    // 5) Input/output at air–prism interfaces (0–1 and 1–0)
    // ============================================================
    cd tp01in = 2.0 * n0in * Cos_Sigma0in
              / ( n1in*Cos_Sigma0in + n0in*Cos_Sigma1in );
    cd tp01vi = 2.0 * n0vi * Cos_Sigma0vi
              / ( n1vi*Cos_Sigma0vi + n0vi*Cos_Sigma1vi );
    cd ts01in = 2.0 * n0in * Cos_Sigma0in
              / ( n0in*Cos_Sigma0in + n1in*Cos_Sigma1in );
    cd ts01vi = 2.0 * n0vi * Cos_Sigma0vi
              / ( n0vi*Cos_Sigma0vi + n1vi*Cos_Sigma1vi );

    cd tp10su = 2.0 * n1su * Cos_Sigma1su
              / ( n0su*Cos_Sigma1su + n1su*Cos_Sigma0su );
    cd ts10su = 2.0 * n1su * Cos_Sigma1su
              / ( n1su*Cos_Sigma1su + n0su*Cos_Sigma0su );

    // ============================================================
    // 6) Fresnel prefactors exactly as in MATLAB scoring
    // ============================================================
    cd Fsspyyz1 =
        ts10su * LyyPoWasu *
        ts01vi * LyyPoWavi *
        tp01in * LzzPoWain *
        Sin_Phi1in;

    cd Fpppxxz1 =
       -tp10su * LxxPoWasu * Cos_Phi1su *
        tp01vi * LxxPoWavi * Cos_Phi1vi *
        tp01in * LzzPoWain * Sin_Phi1in;

    cd Fpppxzx1 =
       -tp10su * LxxPoWasu * Cos_Phi1su *
        tp01vi * LzzPoWavi * Sin_Phi1vi *
        tp01in * LxxPoWain * Cos_Phi1in;

    cd Fpppzxx1 =
        tp10su * LzzPoWasu * Sin_Phi1su *
        tp01vi * LxxPoWavi * Cos_Phi1vi *
        tp01in * LxxPoWain * Cos_Phi1in;

    cd Fpppzzz1 =
        tp10su * LzzPoWasu * Sin_Phi1su *
        tp01vi * LzzPoWavi * Sin_Phi1vi *
        tp01in * LzzPoWain * Sin_Phi1in;

    FF.F_ssp_yyz = std::abs(Fsspyyz1);
    FF.F_ppp_xxz = std::abs(Fpppxxz1);
    FF.F_ppp_xzx = std::abs(Fpppxzx1);
    FF.F_ppp_zxx = std::abs(Fpppzxx1);
    FF.F_ppp_zzz = std::abs(Fpppzzz1);

    return FF;
}
