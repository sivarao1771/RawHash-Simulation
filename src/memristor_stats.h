#ifndef MEMRISTOR_STATS_H
#define MEMRISTOR_STATS_H
#ifdef USE_MEMRISTOR

#include <cstdint>
#include <cstdio>
#include <atomic>
#include "memristor_device_params.h"

struct MemristorStats {
    std::atomic<uint64_t> crossbar_init_cycles{0};
    std::atomic<uint64_t> vmm_read_cycles{0};
    std::atomic<uint64_t> crossbar_comparator_decisions{0};
    std::atomic<uint64_t> cam_xor_decisions{0};
    std::atomic<uint64_t> cam_set_cycles{0};
    std::atomic<uint64_t> cam_reset_cycles{0};
    std::atomic<uint64_t> cam_rows_written{0};
    std::atomic<uint64_t> cam_search_cycles{0};
    std::atomic<uint64_t> cam_sa_firings{0};

    void print(const char* phase) const {
        using namespace MemristorParams;
        double t_init   = crossbar_init_cycles.load() * T_CROSSBAR_INIT_NS;
        double t_vmm    = vmm_read_cycles.load()      * T_VMM_NS;
        double t_write  = cam_set_cycles.load()       * T_SET_NS
                        + cam_reset_cycles.load()     * T_RESET_NS;
        double t_search = cam_search_cycles.load()    * T_CAM_SEARCH_NS;
        double t_total  = t_init + t_vmm + t_write + t_search;

        fprintf(stderr, "\n=== Memristor Cycle Counts [%s] ===\n", phase);
        fprintf(stderr, "CROSSBAR\n");
        fprintf(stderr, "  Init write cycles (512 cells)  : %llu\n",
                (unsigned long long)crossbar_init_cycles.load());
        fprintf(stderr, "  VMM read cycles                : %llu\n",
                (unsigned long long)vmm_read_cycles.load());
        fprintf(stderr, "  Comparator decisions (VMM x 64): %llu\n",
                (unsigned long long)crossbar_comparator_decisions.load());
        fprintf(stderr, "CAM WRITES\n");
        fprintf(stderr, "  SET   cycles (1-bits, HRS->LRS): %llu\n",
                (unsigned long long)cam_set_cycles.load());
        fprintf(stderr, "  RESET cycles (0-bits, LRS->HRS): %llu\n",
                (unsigned long long)cam_reset_cycles.load());
        fprintf(stderr, "  Rows programmed                : %llu\n",
                (unsigned long long)cam_rows_written.load());
        fprintf(stderr, "CAM SEARCHES\n");
        fprintf(stderr, "  Search cycles (parallel)       : %llu\n",
                (unsigned long long)cam_search_cycles.load());
        fprintf(stderr, "  XOR decisions (rows x 64)      : %llu\n",
                (unsigned long long)cam_xor_decisions.load());
        fprintf(stderr, "  SA firings (rows x searches)   : %llu\n",
                (unsigned long long)cam_sa_firings.load());
        fprintf(stderr, "PROJECTED TIME (He et al. 2025)\n");
        fprintf(stderr, "  Crossbar init  : %.2f ns (%.4f ms)\n", t_init,   t_init/1e6);
        fprintf(stderr, "  VMM total      : %.2f ns (%.4f ms)\n", t_vmm,    t_vmm/1e6);
        fprintf(stderr, "  CAM write      : %.2f ns (%.4f ms)\n", t_write,  t_write/1e6);
        fprintf(stderr, "  CAM search     : %.2f ns (%.4f ms)\n", t_search, t_search/1e6);
        fprintf(stderr, "  TOTAL          : %.2f ns (%.4f ms)\n", t_total,  t_total/1e6);
        fprintf(stderr, "===================================\n\n");
    }

    void reset() {
        crossbar_init_cycles.store(0); vmm_read_cycles.store(0);
        crossbar_comparator_decisions.store(0);
        cam_xor_decisions.store(0);
        cam_set_cycles.store(0); cam_reset_cycles.store(0);     cam_rows_written.store(0);
        cam_search_cycles.store(0);    cam_sa_firings.store(0);
    }
};

extern MemristorStats g_mem_stats;

#ifdef __cplusplus
extern "C" {
#endif
    void mem_stats_print(const char* phase);
    void mem_stats_reset(void);
#ifdef __cplusplus
}
#endif

#endif // USE_MEMRISTOR
#endif // MEMRISTOR_STATS_H