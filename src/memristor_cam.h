#ifndef MEMRISTOR_CAM_H
#define MEMRISTOR_CAM_H

#include <vector>
#include <bitset>
#include <cstdint>

class MemristorCAM {
public:
    static constexpr size_t HASH_BITS = 64;

    struct Row {
        std::bitset<HASH_BITS> hash;
        uint64_t metadata;   // Encoded: id, position, strand (matches mm128_t::y)
    };

    // threshold = allowed Hamming distance for a match
    MemristorCAM(int mismatch_threshold = 7);

    // Store a reference hash and its associated metadata
    void program(const std::bitset<HASH_BITS>& hash, uint64_t metadata);

    // Search for all rows with Hamming distance <= threshold.
    // Returns a vector of the stored metadata for each matching row.
    std::vector<uint64_t> search(const std::bitset<HASH_BITS>& query) const;

    // Clear all stored rows
    void clear();

    size_t size() const { return rows_.size(); }

private:
    std::vector<Row> rows_;
    int threshold_;
};

#endif

#ifdef __cplusplus
extern "C" {
#endif
int ri_idx_cam_search(uint64_t query_hash, uint64_t* out_metadata, int max_results);
#ifdef __cplusplus
}
#endif