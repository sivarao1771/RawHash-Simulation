#include "memristor_crossbar.h"
#include <cmath>
#include <algorithm>

MemristorCrossbar::MemristorCrossbar(int input_dim, int hash_bits, double conductance_std)
    : input_dim_(input_dim), hash_bits_(hash_bits), rng_(std::random_device{}()) {

    conductances_.resize(input_dim_);
    std::lognormal_distribution<double> dist(0.0, conductance_std);

    for (int i = 0; i < input_dim_; ++i) {
        conductances_[i].resize(hash_bits_);
        for (int j = 0; j < hash_bits_; ++j) {
            conductances_[i][j] = dist(rng_);
        }
    }
}

std::vector<bool> MemristorCrossbar::hash(const std::vector<float>& events) const {
    if ((int)events.size() != input_dim_) {
        throw std::runtime_error("Event vector size does not match crossbar input dimension");
    }

    std::vector<bool> result(hash_bits_);
    for (int j = 0; j < hash_bits_; ++j) {
        double dot = 0.0;
        for (int i = 0; i < input_dim_; ++i) {
            dot += events[i] * conductances_[i][j];
        }
        // Threshold at zero (since we use zero‑mean conductances)
        result[j] = (dot > 0.0);
    }
    return result;
}

uint64_t MemristorCrossbar::hashToU64(const std::vector<float>& events) const {
    auto bits = hash(events);
    uint64_t val = 0;
    for (size_t i = 0; i < bits.size() && i < 64; ++i) {
        if (bits[i]) val |= (1ULL << i);
    }
    return val;
}