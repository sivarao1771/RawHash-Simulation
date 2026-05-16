#include "memristor_cam.h"

MemristorCAM::MemristorCAM(int mismatch_threshold)
    : threshold_(mismatch_threshold) {}

void MemristorCAM::program(const std::bitset<HASH_BITS>& hash, uint64_t metadata) {
    #ifdef USE_MEMRISTOR
        uint64_t ones = static_cast<uint64_t>(hash.count());
        g_mem_stats.cam_set_cycles.fetch_add(ones,           std::memory_order_relaxed);
        g_mem_stats.cam_reset_cycles.fetch_add(HASH_BITS - ones, std::memory_order_relaxed);
        g_mem_stats.cam_rows_written.fetch_add(1,            std::memory_order_relaxed);
    #endif
    rows_.push_back({hash, metadata});
}

std::vector<uint64_t> MemristorCAM::search(const std::bitset<HASH_BITS>& query) const {
    #ifdef USE_MEMRISTOR
        g_mem_stats.cam_search_cycles.fetch_add(1, std::memory_order_relaxed);
        g_mem_stats.cam_sa_firings.fetch_add(
            static_cast<uint64_t>(rows_.size()), std::memory_order_relaxed);
        g_mem_stats.cam_xor_decisions.fetch_add(
            static_cast<uint64_t>(rows_.size()) * HASH_BITS, std::memory_order_relaxed);
    #endif
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

void ri_idx_cam_program(uint64_t hash_val, uint64_t metadata) {
    if (!g_cam) return;
    std::bitset<64> hash_bits(hash_val);
    g_cam->program(hash_bits, metadata);
}

// Optional: if we want to restore CAM query later, we can add a C‑callable function here

} // extern "C"