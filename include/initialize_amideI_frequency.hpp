#ifndef INITIALIZE_AMIDEI_FREQUENCY_HPP
#define INITIALIZE_AMIDEI_FREQUENCY_HPP

#include <vector>

// Holds per–amide mode vibrational parameters
struct AmideIFreq {
    double freq;    // fundamental frequency (cm^-1)
    double anharm;  // anharmonicity (cm^-1)
};

// Initialize frequency + anharmonicity for all modes
std::vector<AmideIFreq> initialize_amideI_frequency(
    std::size_t n_modes,
    double center_freq,
    double anharm = 12.0
);

#endif
