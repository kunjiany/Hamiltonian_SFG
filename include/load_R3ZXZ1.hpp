#ifndef LOAD_R3ZXZ1_HPP
#define LOAD_R3ZXZ1_HPP

#include <H5Cpp.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>
#include <iostream>
#include <mutex>


static std::mutex h5_mutex;

struct R3Matrix {
    double m[27][27];
};

class R3Database {
private:
    H5::H5File file;
    H5::DataSet dset_R3;

    bool matlab_mode = false;
    bool alt_mode    = false;

    hsize_t dims[4];

public:
    size_t nTheta;
    size_t nPsi;

    R3Database(const std::string& fname)
        : file(fname, H5F_ACC_RDONLY)
    {
        dset_R3 = file.openDataSet("R3");
        H5::DataSpace dsp = dset_R3.getSpace();
        dsp.getSimpleExtentDims(dims);

        if (dims[0] == 181 && dims[1] == 361 && dims[2] == 27 && dims[3] == 27) {
            matlab_mode = true;
            nTheta = 181; nPsi = 361;
            std::cout << "[R3 loader] MATLAB mode detected.\n";
        }
        else if (dims[0] == 19 && dims[1] == 37 && dims[2] == 27 && dims[3] == 27) {
            matlab_mode = true;
            nTheta = 19; nPsi = 37;
            std::cout << "[R3 loader] 10-degree grid detected.\n";
        }
        else if (dims[0] == 27 && dims[1] == 27 && dims[2] == 361 && dims[3] == 181) {
            alt_mode = true;
            nTheta = 181; nPsi = 361;
            std::cout << "[R3 loader] ALT mode detected.\n";
        }
        else {
            throw std::runtime_error("ERROR: Unrecognized R3 dimensions.");
        }
    }

    R3Matrix get_R(double psi_deg, double theta_deg)
    {
        std::lock_guard<std::mutex> lock(h5_mutex);
        
        psi_deg   = fmod(fmod(psi_deg, 360.0) + 360.0, 360.0);
        theta_deg = fmod(fmod(theta_deg,180.0) + 180.0,180.0);

        int iPsi   = round(psi_deg   / 10.0);
        int iTheta = round(theta_deg / 10.0);

        iPsi   = std::max(0, std::min((int)nPsi-1,   iPsi));
        iTheta = std::max(0, std::min((int)nTheta-1, iTheta));

        std::vector<double> buf(27 * 27);

        H5::DataSpace full = dset_R3.getSpace();

        if (matlab_mode) {
            hsize_t offset[4] = { (hsize_t)iTheta, (hsize_t)iPsi, 0, 0 };
            hsize_t count[4]  = { 1, 1, 27, 27 };

            full.selectHyperslab(H5S_SELECT_SET, count, offset);

            hsize_t memdims[4] = { 1,1,27,27 };
            H5::DataSpace mem(4, memdims);

            dset_R3.read(buf.data(), H5::PredType::NATIVE_DOUBLE, mem, full);
        }
        else if (alt_mode) {
            hsize_t offset[4] = { 0, 0, (hsize_t)iPsi, (hsize_t)iTheta };
            hsize_t count[4]  = { 27, 27, 1, 1 };

            full.selectHyperslab(H5S_SELECT_SET, count, offset);

            hsize_t memdims[4] = { 27,27,1,1 };
            H5::DataSpace mem(4, memdims);

            dset_R3.read(buf.data(), H5::PredType::NATIVE_DOUBLE, mem, full);
        }

        R3Matrix R;

        for (int r = 0; r < 27; r++)
            for (int c = 0; c < 27; c++)
                R.m[r][c] = buf[c * 27 + r];

        return R;   // ★ FIXED: always return
    }
};

#endif
