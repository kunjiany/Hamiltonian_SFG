#include "Fresnel_Read.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>

// ---------------- Trim helpers ----------------
static inline std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// Remove comments (# and ;)
static inline std::string remove_comments(const std::string& s) {
    size_t p1 = s.find('#');
    size_t p2 = s.find(';');
    size_t p = std::min(
        (p1 == std::string::npos ? s.size() : p1),
        (p2 == std::string::npos ? s.size() : p2)
    );
    return s.substr(0, p);
}

// --------------------------------------------------
//      MAIN FRESNEL READER
// --------------------------------------------------
FresnelParams Read_Fresnel_File(const std::string& filename)
{
    FresnelParams fp;
    std::unordered_map<std::string, std::string> kv;

    std::ifstream fin(filename);
    if (!fin) {
        std::cerr << "ERROR: Cannot open Fresnel file: " << filename << "\n";
        exit(1);
    }

    std::string line;
    while (std::getline(fin, line)) {

        line = remove_comments(line);
        line = trim(line);
        if (line.empty()) continue;

        size_t eq = line.find('=');
        if (eq == std::string::npos) continue;

        std::string key = trim(line.substr(0, eq));
        std::string val = trim(line.substr(eq + 1));
        kv[key] = val;
    }

    auto get = [&](const std::string& key) -> std::string {
        if (!kv.count(key)) {
            std::cerr << "ERROR: Missing Fresnel parameter: " << key << "\n";
            exit(1);
        }
        return kv[key];
    };

    // Required fields
    fp.geometry   = get("geometry");

    fp.A_vis_deg  = std::stod(get("A_vis"));
    fp.A_ir_deg   = std::stod(get("A_ir"));

    fp.n0_vis = std::stod(get("n0_vis"));
    fp.n1_vis = std::stod(get("n1_vis"));
    fp.n2_vis = std::stod(get("n2_vis"));
    fp.n3_vis = std::stod(get("n3_vis"));

    fp.n0_ir  = std::stod(get("n0_ir"));
    fp.n1_ir  = std::stod(get("n1_ir"));
    fp.n2_ir  = std::stod(get("n2_ir"));
    fp.n3_ir  = std::stod(get("n3_ir"));

    fp.n0_sfg = std::stod(get("n0_sfg"));
    fp.n1_sfg = std::stod(get("n1_sfg"));
    fp.n2_sfg = std::stod(get("n2_sfg"));
    fp.n3_sfg = std::stod(get("n3_sfg"));

    fp.lambda_vis_nm = std::stod(get("lambda_vis_nm"));
    fp.lambda_ir_nm  = std::stod(get("lambda_ir_nm"));
    fp.lambda_sfg_nm = std::stod(get("lambda_sfg_nm"));
    fp.polymer_thickness_nm = std::stod(get("polymer_thickness_nm"));

    return fp;
}
