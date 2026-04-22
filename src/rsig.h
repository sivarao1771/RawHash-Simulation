#ifndef RSIG_H
#define RSIG_H

#ifndef NPOD5RH
#define NPOD5RH   // Trick the IDE into skipping POD5 includes
#endif

#include "rutils.h"

#ifdef __cplusplus
#ifndef NHDF5RH
#include "hdf5_tools.hpp"
#endif
#endif

#ifndef NPOD5RH
#include "pod5_format/c_api.h"
#endif

#ifndef NSLOW5RH
#include <slow5/slow5.h>
#include <slow5/slow5_mt.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ri_sig_s {
    uint32_t rid, l_sig;
    char *name;
    uint64_t offset;
    float* sig;
} ri_sig_t;

typedef struct { size_t n, m; ri_sig_t **a; } rhsig_v;

typedef struct ri_sig_file_s {
    int num_read;
    int cur_read;

    char** raw_path;
    char** ch_path;
#ifdef __cplusplus
#ifndef NHDF5RH
    hdf5_tools::File* fp;
#endif
#endif

    unsigned long int pod5_batch_count;
    unsigned long int pod5_batch;
    unsigned long int pod5_row_count;
    unsigned long int pod5_row;
#ifndef NPOD5RH
    Pod5ReadRecordBatch_t* batch;
    Pod5FileReader_t* pp;
#endif

#ifndef NSLOW5RH
    slow5_mt_t* slow5_mt;
    slow5_batch_t* slow5_batch;
    slow5_file_t* sp;
#endif
} ri_sig_file_t;

void ri_sig_close(ri_sig_file_t *fp);
ri_sig_file_t *open_sig(const char *fn, int io_n_threads);
ri_sig_file_t **open_sigs(int n, const char **fn, int io_n_threads);

void ri_seq_to_sig(const char *str,
                   int len,
                   const ri_pore_t* pore,
                   const int pore_kmer,
                   const int strand,
                   uint32_t* s_len,
                   float* s_values);

void ri_read_sig(ri_sig_file_t* fp, ri_sig_t* s, int io_n_threads);
void find_sfiles(const char *A, ri_char_v *fnames);

#ifdef __cplusplus
}
#endif
#endif