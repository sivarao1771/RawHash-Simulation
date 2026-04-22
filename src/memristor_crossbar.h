#ifndef MEMRISTOR_CROSSBAR_H
#define MEMRISTOR_CROSSBAR_H

#include <vector>
#include <random>
#include <bitset>

class MemristorCrossbar {
public:
    // Create a crossbar with 'input_dim' rows and 'hash_bits' columns.
    // The conductance matrix is initialised with lognormal noise to mimic device variability.
    MemristorCrossbar(int input_dim, int hash_bits, double conductance_std = 5e-6);

    // Hash a vector of analog events into a binary fingerprint.
    // Returns a vector<bool> where each element is one bit of the hash.
    std::vector<bool> hash(const std::vector<float>& events) const;

    // Convenience: return the hash as a 64‑bit integer (compatible with RawHash's mm128_t::x).
    uint64_t hashToU64(const std::vector<float>& events) const;

    // Accessors
    int inputDim() const { return input_dim_; }
    int hashBits() const { return hash_bits_; }

private:
    int input_dim_;
    int hash_bits_;
    std::vector<std::vector<double>> conductances_; // [input_dim][hash_bits]
    mutable std::mt19937 rng_;                       // for reproducibility
};

#endif