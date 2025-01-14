/*
Description: fmt-idxseed (fm-index twice extend in one search step)

Copyright : All right reserved by 

Author : 
Date : 2023/12/24
*/
#ifndef FMT_INDEX_H_
#define FMT_INDEX_H_

#include <stdint.h>
#include <stddef.h>
#include "bwt.h"
#include "utils.h"

#define FMT_OCC_INTV_SHIFT 8
#define FMT_OCC_INTERVAL (1 << FMT_OCC_INTV_SHIFT)
#define FMT_OCC_INTV_MASK (FMT_OCC_INTERVAL - 1)

#define FMT_MID_INTV_SHIFT 6
#define FMT_MID_INTERVAL (1 << FMT_MID_INTV_SHIFT)
#define FMT_MID_INTV_MASK (FMT_MID_INTERVAL - 1)

// c（），
#define fmt_set_intv(fmt, c, ik) ((ik).x[0] = (fmt)->L2[(int)(c)] + 1, (ik).x[2] = (fmt)->L2[(int)(c) + 1] - (fmt)->L2[(int)(c)], (ik).x[1] = (fmt)->L2[3 - (c)] + 1, (ik).info = 0)
// k（bwt str（$））check point occ（kOCC_INTERVAL）

#if FMT_MID_INTERVAL == 8
#define fmt_occ_intv(b, k) ((b)->bwt + (k) / FMT_OCC_INTERVAL * (FMT_OCC_INTERVAL / 8 + 144))
#elif FMT_MID_INTERVAL == 16
#define fmt_occ_intv(b, k) ((b)->bwt + (k) / FMT_OCC_INTERVAL * (FMT_OCC_INTERVAL / 8 + 80))
#elif FMT_MID_INTERVAL == 32
//#define fmt_occ_intv(b, k) ((b)->bwt + (k) / FMT_OCC_INTERVAL * (FMT_OCC_INTERVAL / 8 + 48))
#define fmt_occ_intv(b, k) ((b)->bwt + ((k) >> 8) * 80)
#elif FMT_MID_INTERVAL == 64
//#define fmt_occ_intv(b, k) ((b)->bwt + (k) / FMT_OCC_INTERVAL * (FMT_OCC_INTERVAL / 8 + 32))
#define fmt_occ_intv(b, k) ((b)->bwt + ((k) >> 8 << 6))
#elif FMT_MID_INTERVAL == 128
#define fmt_occ_intv(b, k) ((b)->bwt + (k) / FMT_OCC_INTERVAL * (FMT_OCC_INTERVAL / 8 + 24))
#else
#define fmt_occ_intv(b, k) ((b)->bwt + (k) / FMT_OCC_INTERVAL * (FMT_OCC_INTERVAL / 8 + 20))
#endif

// valbwt basebpre-bwtT G C A（32（8bit）），occocc
#define __fmt_occ_e2_aux2(fmt, b, val) \
    ((fmt)->cnt_occ[(b)][(val) & 0xff] + (fmt)->cnt_occ[b][(val) >> 8 & 0xff] + (fmt)->cnt_occ[b][(val) >> 16 & 0xff] + (fmt)->cnt_occ[b][(val) >> 24])

#define __fmt_mid_sum(x) \
    ((x) >> 24 & 0xff) + ((x) >> 16 & 0xff) + ((x) >> 8 & 0xff) + ((x) & 0xff)

// sa
#define SA_INTV 4

// fm-index, extend twice in one search step (one memory access)
typedef struct
{
    bwtint_t primary;     // S^{-1}(0), or the primary index of BWT
    bwtint_t sec_primary; // second primary line
    bwtint_t L2[5];       // C(), cumulative count
    bwtint_t seq_len;     // sequence length
    bwtint_t bwt_size;    // size of bwt, about seq_len/4
    uint32_t *bwt;        // BWT
    // occurance array, separated to two parts
    uint32_t cnt_occ[16][256]; // 16-24b（）occ，8-16bocc，0-8aocc，ba
    uint8_t sec_bcp;           // base couple for sec primary line, AA=>0, AC=>1 ... TT=>15
    uint8_t first_base;        // 2bitint，0,1,2,3
    uint8_t last_base;         // dollarbase
    // ref pac
    bwtint_t l_pac; // 
    uint8_t *pac; // 2bit
    // kmerfmt
    KmerHash kmer_hash;
    // suffix array
    int sa_intv;
    bwtint_t n_sa;
    uint8_t *sa;
} FMTIndex;

// fmt
void dump_fmt(const char *fn, const FMTIndex *fmt);
// fmt
FMTIndex *fmt_restore_fmt(const char *fn);
void fmt_create_kmer_index(FMTIndex *fmt);
// kmer hash
void fmt_dump_kmer_idx(const char *fn, const KmerHash *kh);
// kmer hash
KmerHash fmt_restore_kmer_idx(const char *fn);
// sa
void fmt_restore_sa(const char *fn, FMTIndex *fmt);
// interval-bwtfmt-index
FMTIndex *create_fmt_from_bwt(bwt_t *bwt);
// ，bwt basebpre-bwt strocc
void fmt_direct_e2_occ(const FMTIndex *fmt, bwtint_t k, int b1, int b2, bwtint_t cnt[4]);
void fmt_direct2_e2_occ(const FMTIndex *fmt, bwtint_t k, int b1, int b2, bwtint_t cnt[4]);
void fmt_e2_occ(const FMTIndex *fmt, bwtint_t k, int b1, int b2, bwtint_t cnt[4]);
// 
void fmt_direct2_extend2(const FMTIndex *fmt, bwtintv_t *ik, bwtintv_t *ok1, bwtintv_t *ok2, int is_back, int b1, int b2);
void fmt_direct_extend2(const FMTIndex *fmt, bwtintv_t *ik, bwtintv_t *ok1, bwtintv_t *ok2, int is_back, int b1, int b2);
void fmt_extend2(const FMTIndex *fmt, bwtintv_t *ik, bwtintv_t *ok1, bwtintv_t *ok2, int is_back, int b1, int b2);
// 
void fmt_direct_extend1(const FMTIndex *fmt, bwtintv_t *ik, bwtintv_t *ok, int is_back, int b1);
void fmt_extend1(const FMTIndex *fmt, bwtintv_t *ik, bwtintv_t *ok, int is_back, int b1);
// KMER_LEN，
// void gen_all_seq(char **seq_arr, int kmer_len);
// kmerposfmt
void kmer_setval_at(uint8_t *mem_addr, bwtintv_t ik, int pos);
// kmerfmt
void kmer_getval_at(uint8_t *mem_addr, bwtintv_t *ok, int pos);
void fmt_kmer_get(const FMTIndex *fmt, bwtintv_t *ok, uint32_t qbit, int pos);
// smem（seed）
int fmt_smem(const FMTIndex *fmt, int len, const uint8_t *q, int x, int min_intv, int min_seed_len, bwtintv_t *lm, bwtintv_v *mem, bwtintv_v *tmpvec);
int fmt_smem_2(const FMTIndex *fmt, int len, const uint8_t *q, int x, int min_intv, int min_seed_len, bwtintv_t *lm, bwtintv_v *mem, bwtintv_v *tmpvec, uint32_v qev);

int fmt_seed_strategy1(const FMTIndex *fmt, int len, const uint8_t *q, int x, int min_len, int max_intv, bwtintv_t *mem);

bwtint_t fmt_sa(const FMTIndex *fmt, bwtint_t k);
bwtint_t fmt_sa_offset(const FMTIndex *fmt, bwtint_t k);

void fmt_gen_cnt_occ(FMTIndex *fmt);
#endif