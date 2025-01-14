/*
Description: profiling related data

Copyright : All right reserved by 

Author : 
Date : 2024/04/06
*/
#ifndef PROFILING_H_
#define PROFILING_H_

#include <stdint.h>

#define USE_RDTSC 1

#define LIM_THREAD 128
#define LIM_THREAD_PROF_TYPE 128
#define LIM_GLOBAL_PROF_TYPE 128
#define LIM_THREAD_DATA_TYPE 128
#define LIM_GLOBAL_DATA_TYPE 128

#ifdef SHOW_PERF
extern uint64_t proc_freq;
extern uint64_t tprof[LIM_THREAD_PROF_TYPE][LIM_THREAD];
extern uint64_t gprof[LIM_GLOBAL_PROF_TYPE];
#endif

#ifdef SHOW_DATA_PERF
extern int64_t tdat[LIM_THREAD_DATA_TYPE][LIM_THREAD];
extern int64_t gdat[LIM_GLOBAL_DATA_TYPE];
#endif

#ifdef SHOW_PERF
#if USE_RDTSC
#define PROF_START(tmp_time) const uint64_t prof_tmp_##tmp_time = __rdtsc()
#define PROF_END(result, tmp_time) result += __rdtsc() - prof_tmp_##tmp_time
#else
#define PROF_START(tmp_time) uint64_t prof_tmp_##tmp_time = realtime_msec()
#define PROF_END(result, tmp_time) result += realtime_msec() - prof_tmp_##tmp_time
#endif
#else
#define PROF_START(tmp_time)
#define PROF_END(result, tmp_time)
#endif

extern uint64_t calc_num;
extern uint64_t more_num;
extern uint64_t not_more_num;

// GLOBAL
enum {
    G_ALL = 0,
    G_PIPELINE,
    G_READ,
    G_WRITE,
    G_COMPUTE,
    G_PREPARE,
    G_LOAD_IDX,
    G_MEM_PREPARE,
    G_MEM_KERNEL,
    G_MEM_PESTAT,
    G_MEM_SAM,
    G_KSW_LOOP,
    G_KSW_END_LOOP,
    G_parse_seq,
    G_read_seq,
    G_read_wait_1,
    G_read_wait_2,
    G_comp_wait_1,
    G_comp_wait_2,
    G_write_wait_1,
    G_write_wait_2
};

// THREAD
enum
{
    T_SEED_ALL = 0,
    T_CHAIN_ALL,
    T_ALN_ALL,
    T_INS_SIZE,
    T_SEED_1,
    T_SEED_2,
    T_SEED_3,
    T_SAL,
    T_GEN_CHAIN,
    T_FLT_CHAIN,
    T_FLT_CHANNED_SEEDS,
    T_READ_SA,
    T_BSW,
    T_BSW_ALL,
    T_SAM_MATESW,
    T_SAM_REG2ALN,
    T_SEED_3_1,
    T_SEED_3_2,
    T_SEED_3_3,
    T_SEEDING,
    T_EXTENSION,
    T_SAM,
};

int display_stats(int);

#endif