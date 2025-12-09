#include "Read_Input.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <cctype>

// Trim helper
static inline std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if(start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// Remove comments starting with ';' or '#'
static inline std::string remove_comments(const std::string& line) {
    size_t pos_semi = line.find(';');
    size_t pos_hash = line.find('#');
    size_t pos = std::min(
        pos_semi == std::string::npos ? line.size() : pos_semi,
        pos_hash == std::string::npos ? line.size() : pos_hash
    );
    return line.substr(0, pos);
}

InputParams Read_Input(const std::string &filename)
{
    InputParams p;
    std::unordered_map<std::string, std::string> kv;

    std::ifstream fin(filename);
    if(!fin){
        std::cerr << "ERROR: Cannot open " << filename << "\n";
        exit(1);
    }

    std::string line;
    while(std::getline(fin, line)) {

        // Remove inline comments like GROMACS/LAMMPS
        line = remove_comments(line);

        // Trim whitespace
        line = trim(line);
        if(line.empty()) continue;

        // Split at '='
        size_t eq = line.find('=');
        if(eq == std::string::npos) continue;

        std::string key = trim(line.substr(0, eq));
        std::string val = trim(line.substr(eq + 1));

        kv[key] = val;
    }

    // ---------- Assign required parameters ----------
    if(kv.count("PDB_file"))
        p.pdbFile = kv["PDB_file"];
    else {
        std::cerr << "ERROR: Missing required input: PDB_file\n";
        exit(1);
    }

    if(kv.count("center_freq"))
        p.centerFreq = std::stod(kv["center_freq"]);
    else {
        std::cerr << "ERROR: Missing required input: center_freq\n";
        exit(1);
    }

    if(kv.count("layer"))
        p.layer = std::stoi(kv["layer"]);
    else {
        std::cerr << "ERROR: Missing required input: layer\n";
        exit(1);
    }

    if(kv.count("tilt_start"))
        p.tilt_start = std::stod(kv["tilt_start"]);
    else {
        std::cerr << "ERROR: Missing required input: tilt_start\n";
        exit(1);
    }

    if(kv.count("tilt_end"))
        p.tilt_end = std::stod(kv["tilt_end"]);
    else {
        std::cerr << "ERROR: Missing required input: tilt_end\n";
        exit(1);
    }

    if(kv.count("tilt_points"))
        p.tilt_points = std::stoi(kv["tilt_points"]);
    else {
        std::cerr << "ERROR: Missing required input: tilt_points\n";
        exit(1);
    }

    if(kv.count("twist_start"))
        p.twist_start = std::stod(kv["twist_start"]);
    else {
        std::cerr << "ERROR: Missing required input: twist_start\n";
        exit(1);
    }

    
    if(kv.count("twist_end"))
        p.twist_end= std::stod(kv["twist_end"]);
    else {
        std::cerr << "ERROR: Missing required input: twist_end\n";
        exit(1);
    }

    if(kv.count("twist_points"))
        p.twist_points = std::stoi(kv["twist_points"]);
    else {
        std::cerr << "ERROR: Missing required input: twist_points\n";
        exit(1);
    }

    if(kv.count("use_cutoff")) {
        std::string cutoff = kv["use_cutoff"];
        if(cutoff == "yes") {
            p.use_cutoff = true;
            if(kv.count("cutoff_distance")){
                p.cutoff_distance = std::stod(kv["cutoff_distance"]);
            }
            else {
                std::cerr << "ERROR: Missing required input: cutoff_distance\n";
                exit(1);
            }
        }
        else if(cutoff == "no") {
            p.use_cutoff = false; 
            p.cutoff_distance = 10000000; // very large dummy variable
        }
        else {
            std::cerr << "ERROR: Missing or wrong required input: cutoff_distance\n";
                exit(1);
        }

    }
    else {
        std::cerr << "ERROR: Missing required input: use_cutoff\n";
        exit(1);
    }

    if(kv.count("width"))
        p.width = std::stod(kv["width"]);
    else {
        std::cerr << "ERROR: Missing required input: width\n";
        exit(1);
    }

    if(kv.count("spec_range_start"))
        p.spec_range_start= std::stod(kv["spec_range_start"]);
    else {
        std::cerr << "ERROR: Missing required input: spec_range_start\n";
        exit(1);
    }

    if(kv.count("spec_range_end"))
        p.spec_range_end = std::stod(kv["spec_range_end"]);
    else {
        std::cerr << "ERROR: Missing required input: spec_range_end\n";
        exit(1);
    }

    if(kv.count("spec_range_step"))
        p.spec_range_step = std::stod(kv["spec_range_step"]);
    else {
        std::cerr << "ERROR: Missing required input:spec_range_step\n";
        exit(1);
    }

    if(kv.count("SpectraFolder"))
        p.SpectraFolder = kv["SpectraFolder"];
    else {
        std::cerr << "ERROR: Missing required input: SpectraFolder\n";
        exit(1);
    }

    if(kv.count("SpectraStorePrefix"))
        p.SpectraStorePrefix = kv["SpectraStorePrefix"];
    else {
        std::cerr << "ERROR: Missing required input: SpectraStorePrefix\n";
        exit(1);
    }

    // ---------- Validation ----------
    if(p.centerFreq <= 0) {
        std::cerr << "ERROR: center_freq must be positive.\n";
        exit(1);
    }
    if(p.layer <= 0) {
        std::cerr << "ERROR: layer must be positive.\n";
        exit(1);
    }
    if(p.tilt_start < 0 || p.tilt_start > 180) {
        std::cerr << "ERROR: tilt_start must in range [0,180].\n";
        exit(1);
    }
    if(p.tilt_end < p.tilt_start || p.tilt_end > 180) {
        std::cerr << "ERROR: tilt_end must be [tilt_start,180].\n";
        exit(1);
    }
    if(p.tilt_points<=0) {
        std::cerr << "ERROR: tilt_point need to >0 integer.\n";
        exit(1);
    }
    else {
        if (p.tilt_start == p.tilt_end && p.tilt_points!=1) {
            std::cerr << "ERROR: tilt_point must be 1 if only 1 angle point is given.\n";
            exit(1);           
        }
    }
    if(p.twist_start < 0 || p.twist_start > 360) {
        std::cerr << "ERROR: twist_start must in range [0,360].\n";
        exit(1);
    }
    if(p.twist_end < p.twist_start || p.twist_end > 360) {
        std::cerr << "ERROR: twist_end must be [twist_start,360].\n";
        exit(1);
    }
    if(p.twist_points<=0) {
        std::cerr << "ERROR: twist_point need to >0 integer.\n";
        exit(1);
    }
    else {
        if (p.twist_start == p.twist_end && p.twist_points!=1) {
            std::cerr << "ERROR: twist_point must be 1 if only 1 angle point is given.\n";
            exit(1);           
        }
    }
    if(p.cutoff_distance < 5) {
        std::cerr << "ERROR: cutoff_distance must be greater than or equal to 5.\n";
        exit(1);
    }
    if(p.width <= 0) {
        std::cerr << "ERROR: width must be greater than 0.\n";
        exit(1);
    }
    if(p.spec_range_start<=0) {
        std::cerr << "ERROR: spec_range_start must be greater than 0.\n";
        exit(1);
    }
    if(p.spec_range_start>p.spec_range_end) {
        std::cerr << "ERROR: spec_range_start must be smaller than p.spec_range_end.\n";
        exit(1);
    }
    if(p.spec_range_step<=0) {
        std::cerr << "ERROR: spec_range_step need to >0.\n";
        exit(1);
    }
    return p;
}


