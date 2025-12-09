#include "initialize_amideI_frequency.hpp"

std::vector<AmideIFreq> initialize_amideI_frequency(
    std::size_t n_modes,
    double center_freq,
    double anharm)
{
    std::vector<AmideIFreq> freq_list;
    freq_list.reserve(n_modes);

    for(std::size_t i = 0; i < n_modes; ++i) {
        AmideIFreq fi;
        fi.freq   = center_freq;
        fi.anharm = anharm;
        freq_list.push_back(fi);
    }

    return freq_list;
}
