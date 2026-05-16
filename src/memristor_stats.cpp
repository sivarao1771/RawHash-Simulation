#include "memristor_stats.h"
#ifdef USE_MEMRISTOR
MemristorStats g_mem_stats;
extern "C" {
    void mem_stats_print(const char* phase) { g_mem_stats.print(phase); }
    void mem_stats_reset(void)              { g_mem_stats.reset();       }
}
#endif