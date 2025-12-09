#ifndef READ_INPUT_HPP
#define READ_INPUT_HPP

#include <string>

struct InputParams {
    std::string pdbFile;
    double centerFreq;
    int layer;
    double tilt_start; 
    double tilt_end; 
    int tilt_points; 
    double twist_start;
    double twist_end;
    int twist_points; 
    bool use_cutoff;
    double cutoff_distance;
    double width;           
    double spec_range_start; 
    double spec_range_end;   
    double spec_range_step;  
    std::string SpectraFolder;      
    std::string SpectraStorePrefix;  
};

InputParams Read_Input(const std::string &filename);

#endif

