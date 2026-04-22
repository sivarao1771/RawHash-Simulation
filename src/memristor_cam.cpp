#include "memristor_cam.h"

MemristorCAM::MemristorCAM(int mismatch_threshold)
    : threshold_(mismatch_threshold) {}

void MemristorCAM::program(const std::bitset<HASH_BITS>& hash, uint64_t metadata) {
    rows_.push_back({hash, metadata});
}

std::vector<uint64_t> MemristorCAM::search(const std::bitset<HASH_BITS>& query) const {
    std::vector<uint64_t> matches;
    for (const auto& row : rows_) {
        auto diff = (query ^ row.hash).count();
        if (diff <= static_cast<size_t>(threshold_)) {
            matches.push_back(row.metadata);
        }
    }
    return matches;
}

void MemristorCAM::clear() {
    rows_.clear();
}

// C wrapper for the memristor CAM (used by C code via extern "C")
static MemristorCAM* g_cam = nullptr;

extern "C" {

void ri_idx_enable_memristor(int threshold) {
    if (!g_cam) g_cam = new MemristorCAM(threshold);
}

void ri_idx_disable_memristor() {
    delete g_cam;
    g_cam = nullptr;
}
int ri_idx_cam_search(uint64_t query_hash, uint64_t* out_metadata, int max_results) {
    if (!g_cam) return 0;
    std::bitset<64> query(query_hash);
    auto matches = g_cam->search(query);
    int n = 0;
    for (auto meta : matches) {
        if (n >= max_results) break;
        out_metadata[n++] = meta;
    }
    return n;
}

// Optional: if we want to restore CAM query later, we can add a C‑callable function here

} // extern "C"