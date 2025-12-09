#ifndef FRESNEL_READ_HPP
#define FRESNEL_READ_HPP

#include <string>
#include <unordered_map>

// All physical parameters needed by Fresnel_Calculation
struct FresnelParams {
    std::string geometry;

    double A_vis_deg = 0.0;
    double A_ir_deg  = 0.0;

    double n0_vis, n1_vis, n2_vis, n3_vis;
    double n0_ir,  n1_ir, n2_ir, n3_ir;
    double n0_sfg, n1_sfg, n2_sfg, n3_sfg;

    double lambda_vis_nm;
    double lambda_ir_nm;
    double lambda_sfg_nm;

    double polymer_thickness_nm;
};

FresnelParams Read_Fresnel_File(const std::string& filename);

#endif
