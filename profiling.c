/*
Description: profiling related data

Copyright : All right reserved by 

Author : 
Date : 2024/04/06
*/

#include <stdio.h>
#include "utils.h"
#include "profiling.h"

uint64_t calc_num = 0;
uint64_t more_num = 0;
uint64_t not_more_num = 0;

#ifdef SHOW_PERF
uint64_t tprof[LIM_THREAD_PROF_TYPE][LIM_THREAD] = {0};
uint64_t proc_freq = 1000;
uint64_t gprof[LIM_GLOBAL_PROF_TYPE] = {0};
#endif

#ifdef SHOW_DATA_PERF
/*
tdat[0]: read nums
tdat[1]: seed-1 full match nums
*/
int64_t tdat[LIM_THREAD_DATA_TYPE][LIM_THREAD] = {0};
int64_t gdat[LIM_GLOBAL_DATA_TYPE] = {0};

#endif

int find_opt(uint64_t *a, int len, double *max, double *min, double *avg)
{
    int i = 0;
    uint64_t umax = 0, umin = UINT64_MAX, uavg = 0;
    for (i = 0; i < len; i++)
    {
        if (a[i] > umax) umax = a[i];
        if (a[i] < umin) umin = a[i];
        uavg += a[i];
    }
    *avg = uavg * 1.0 / len / proc_freq;
    *max = umax * 1.0 / proc_freq;
    *min = umin * 1.0 / proc_freq;
    return 1;
}

int display_stats(int nthreads)
{
#ifdef SHOW_PERF
    double avg, max, min;

    fprintf(stderr, "[steps in main_mem]\n");
    fprintf(stderr, "time_parse_arg: %0.2lf s\n", gprof[G_PREPARE] * 1.0 / proc_freq);
    fprintf(stderr, "time_load_idx:  %0.2lf s\n", gprof[G_LOAD_IDX] * 1.0 / proc_freq);
    fprintf(stderr, "time_pipeline:  %0.2lf s\n", gprof[G_PIPELINE] * 1.0 / proc_freq);
    fprintf(stderr, "time_all:       %0.2lf s\n", gprof[G_ALL] * 1.0 / proc_freq);
    
    fprintf(stderr, "\n[steps in pipeline]\n");
    fprintf(stderr, "time_read:    %0.2lf s\n", gprof[G_READ] * 1.0 / proc_freq);
    fprintf(stderr, "time_compute: %0.2lf s\n", gprof[G_COMPUTE] * 1.0 / proc_freq);
    fprintf(stderr, "time_write:   %0.2lf s\n", gprof[G_WRITE] * 1.0 / proc_freq);

    fprintf(stderr, "\n[steps in mem_process_seqs]\n");
    fprintf(stderr, "time_mem_prepare: %0.2lf s\n", gprof[G_MEM_PREPARE] * 1.0 / proc_freq);
    fprintf(stderr, "time_mem_kernel:  %0.2lf s\n", gprof[G_MEM_KERNEL] * 1.0 / proc_freq);
    fprintf(stderr, "time_mem_pestat:  %0.2lf s\n", gprof[G_MEM_PESTAT] * 1.0 / proc_freq);
    fprintf(stderr, "time_mem_sam:     %0.2lf s\n", gprof[G_MEM_SAM] * 1.0 / proc_freq);

    fprintf(stderr, "\n[steps in kernel]\n");
    find_opt(tprof[T_SEED_ALL], nthreads, &max, &min, &avg);
    fprintf(stderr, "time_seed_all:     %0.2lf (%0.2lf, %0.2lf) s\n", avg, max, min);
    find_opt(tprof[T_CHAIN_ALL], nthreads, &max, &min, &avg);
    fprintf(stderr, "time_chain_all:    %0.2lf (%0.2lf, %0.2lf) s\n", avg, max, min);
    find_opt(tprof[T_ALN_ALL], nthreads, &max, &min, &avg);
    fprintf(stderr, "time_aln_all:      %0.2lf (%0.2lf, %0.2lf) s\n", avg, max, min);
    find_opt(tprof[T_INS_SIZE], nthreads, &max, &min, &avg);
    fprintf(stderr, "time_ins_size_all: %0.2lf (%0.2lf, %0.2lf) s\n", avg, max, min);

    fprintf(stderr, "\n[steps in seeding]\n");
    find_opt(tprof[T_SEED_1], nthreads, &max, &min, &avg);
    fprintf(stderr, "time_seed_1: %0.2lf s\n", avg);
    find_opt(tprof[T_SEED_2], nthreads, &max, &min, &avg);
    fprintf(stderr, "time_seed_2: %0.2lf s\n", avg);
    find_opt(tprof[T_SEED_3], nthreads, &max, &min, &avg);
    fprintf(stderr, "time_seed_3: %0.2lf s\n", avg);
    
    fprintf(stderr, "\n[steps in chain]\n");
    find_opt(tprof[T_GEN_CHAIN], nthreads, &max, &min, &avg);
    fprintf(stderr, "time_gen_chain:         %0.2lf s\n", avg);
    find_opt(tprof[T_FLT_CHAIN], nthreads, &max, &min, &avg);
    fprintf(stderr, "time_flt_chain:         %0.2lf s\n", avg);
    find_opt(tprof[T_FLT_CHANNED_SEEDS], nthreads, &max, &min, &avg);
    fprintf(stderr, "time_flt_chained_seeds: %0.2lf s\n", avg);
    find_opt(tprof[T_SAL], nthreads, &max, &min, &avg);
    fprintf(stderr, "time_sal:               %0.2lf s\n", avg);
    find_opt(tprof[T_BSW], nthreads, &max, &min, &avg);
    fprintf(stderr, "time_bsw:               %0.2lf s\n", avg);

    fprintf(stderr, "\n[steps in gen sam]\n");
    find_opt(tprof[T_SAM_MATESW], nthreads, &max, &min, &avg);
    fprintf(stderr, "time_mate_sw: %0.2lf s\n", avg);
    find_opt(tprof[T_SAM_REG2ALN], nthreads, &max, &min, &avg);
    fprintf(stderr, "time_reg2aln: %0.2lf s\n", avg);

    fprintf(stderr, "time_ksw_loop:     %0.2lf s\n", gprof[G_KSW_LOOP] * 1.0 / proc_freq);
    fprintf(stderr, "time_ksw_end_loop: %0.2lf s\n", gprof[G_KSW_END_LOOP] * 1.0 / proc_freq);

    fprintf(stderr, "parse seq: %0.2lf s\n", gprof[G_parse_seq] * 1.0 / proc_freq);
    fprintf(stderr, "read seq: %0.2lf s\n", gprof[G_read_seq] * 1.0 / proc_freq);

//    fprintf(stderr, "read_wait_1: %0.2lf s\n", gprof[G_read_wait_1] * 1.0 / proc_freq);
//    fprintf(stderr, "read_wait_2: %0.2lf s\n", gprof[G_read_wait_2] * 1.0 / proc_freq);
//    fprintf(stderr, "comp_wait_1: %0.2lf s\n", gprof[G_comp_wait_1] * 1.0 / proc_freq);
//    fprintf(stderr, "comp_wait_2: %0.2lf s\n", gprof[G_comp_wait_2] * 1.0 / proc_freq);
//    fprintf(stderr, "write_wait_1: %0.2lf s\n", gprof[G_write_wait_1] * 1.0 / proc_freq);
//    fprintf(stderr, "write_wait_2: %0.2lf s\n", gprof[G_write_wait_2] * 1.0 / proc_freq);

    // fprintf(stderr, "real_cal num: %ld\n", calc_num);

    fprintf(stderr, "more num: %ld\n", more_num);
    fprintf(stderr, "no more num: %ld\n", not_more_num);

    fprintf(stderr, "\n");
#endif

    return 1;
}