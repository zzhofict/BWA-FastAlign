/*
Description: fmt-idxseed（fm-index twice search in one time）

Copyright : All right reserved by 

Author : 
Date : 2023/12/24
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fmt_idx.h"
#include "utils.h"
#include "bntseq.h"
#include "kvec.h"
#include "kstring.h"
#include "debug.h"

// occ，pattern
void fmt_gen_cnt_occ(FMTIndex *fmt)
{
    // 0-8：aocc，8-16：bocc，16-24：bocc
    int i, a, b, ti;
    uint32_t oa, ooa, ob, oob;
    for (i = 0; i != 256; ++i) // 
    {
        for (a = 0; a < 4; ++a) // ba
        {
            oa = 0;
            ooa = 0;
            oa += ((i >> 4 & 3) == a) + ((i & 3) == a);
            ooa += ((i >> 4 & 3) > a) + ((i & 3) > a);
            for (b = 0; b < 4; ++b)
            {
                oob = ob = 0;
                oob += ((i >> 6 & 3) > b && (i >> 4 & 3) == a) + ((i >> 2 & 3) > b && (i & 3) == a);
                ob += ((i >> 6 & 3) == b && (i >> 4 & 3) == a) + ((i >> 2 & 3) == b && (i & 3) == a);
                ti = a << 2 | b;
                fmt->cnt_occ[ti][i] = ob << 24 | oob << 16 | oa << 8 | ooa;
            }
        }
    }
}

// fmt-indexcount table，4bwt，0,1,2,3bwtA,C,G,T，pre-bwt
void fmt_gen_cnt_table(uint32_t cnt_table[4][256])
{
    int i, j, k;
    uint32_t x = 0;
    for (i = 0; i != 256; ++i) // 
    {
        for (k = 0; k < 4; ++k) // bwt
        {
            x = 0;                   // for [A,C,G,T][A,C,G,T]
            for (j = 0; j != 4; ++j) // pre-bwt
                x |= (((i >> 6 & 3) == j && (i >> 4 & 3) == k) + ((i >> 2 & 3) == j && (i & 3) == k)) << (j << 3);
            cnt_table[k][i] = x;
        }
    }
}

// fmt
void dump_fmt(const char *fn, const FMTIndex *fmt)
{
    FILE *fp;
    fp = xopen(fn, "wb");
    err_fwrite(&fmt->primary, sizeof(bwtint_t), 1, fp);
    err_fwrite(&fmt->sec_primary, sizeof(bwtint_t), 1, fp);
    err_fwrite(&fmt->sec_bcp, sizeof(uint8_t), 1, fp);
    err_fwrite(&fmt->first_base, sizeof(uint8_t), 1, fp);
    err_fwrite(&fmt->last_base, sizeof(uint8_t), 1, fp);
    err_fwrite(fmt->L2 + 1, sizeof(bwtint_t), 4, fp);
    err_fwrite(fmt->bwt, 4, fmt->bwt_size, fp);
    err_fflush(fp);
    err_fclose(fp);
}

// fmt
FMTIndex *fmt_restore_fmt(const char *fn)
{
    FMTIndex *fmt;
    fmt = (FMTIndex *)calloc(1, sizeof(FMTIndex));
    FILE *fp = xopen(fn, "rb");

    fseek(fp, 0, SEEK_END);
    fmt->bwt_size = (ftell(fp) - sizeof(bwtint_t) * 6 - 3) >> 2; // 32wordsize
    fmt->bwt = (uint32_t *)calloc(fmt->bwt_size, 4);
    fseek(fp, 0, SEEK_SET);
    err_fread_noeof(&fmt->primary, sizeof(bwtint_t), 1, fp);
    err_fread_noeof(&fmt->sec_primary, sizeof(bwtint_t), 1, fp);
    err_fread_noeof(&fmt->sec_bcp, sizeof(uint8_t), 1, fp);
    err_fread_noeof(&fmt->first_base, sizeof(uint8_t), 1, fp);
    err_fread_noeof(&fmt->last_base, sizeof(uint8_t), 1, fp);
    err_fread_noeof(fmt->L2 + 1, sizeof(bwtint_t), 4, fp);
    fread_fix(fp, fmt->bwt_size << 2, fmt->bwt);
    fmt->seq_len = fmt->L2[4];
    err_fclose(fp);
    fmt_gen_cnt_occ(fmt); // ，
    return fmt;
}

// kmer hash
void fmt_dump_kmer_idx(const char *fn, const KmerHash *kh)
{
    FILE *fp;
    fp = xopen(fn, "wb");
    err_fwrite(kh->ke10, 1, (1 << (10 << 1)) * sizeof(KmerEntryArr), fp);
    err_fwrite(kh->ke11, 1, (1 << (11 << 1)) * sizeof(KmerEntry), fp);
    err_fwrite(kh->ke12, 1, (1 << (12 << 1)) * sizeof(KmerEntry), fp);
#if HASH_KMER_LEN > 12
    err_fwrite(kh->ke13, 1, (1 << (13 << 1)) * sizeof(KmerEntry), fp);
#endif
#if HASH_KMER_LEN > 13
    err_fwrite(kh->ke14, 1, (1 << (14 << 1)) * sizeof(KmerEntry), fp);
#endif
    err_fflush(fp);
    err_fclose(fp);
}

// kmer hash
KmerHash fmt_restore_kmer_idx(const char *fn)
{
    FILE *fp = xopen(fn, "rb");
    KmerHash khash;
    KmerHash *kh = &khash;
    int len = 1 << (10 << 1);
    kh->ke10 = (KmerEntryArr *)malloc(len * sizeof(KmerEntryArr));
    fread_fix(fp, len * sizeof(KmerEntryArr), kh->ke10);
    len = 1 << (11 << 1);
    kh->ke11 = (KmerEntry *)malloc(len * sizeof(KmerEntry));
    fread_fix(fp, len * sizeof(KmerEntry), kh->ke11);
    len = 1 << (12 << 1);
    kh->ke12 = (KmerEntry *)malloc(len * sizeof(KmerEntry));
    fread_fix(fp, len * sizeof(KmerEntry), kh->ke12);
#if HASH_KMER_LEN > 12
    len = 1 << (13 << 1);
    kh->ke13 = (KmerEntry *)malloc(len * sizeof(KmerEntry));
    fread_fix(fp, len * sizeof(KmerEntry), kh->ke13);
#endif
#if HASH_KMER_LEN > 13
    len = 1 << (14 << 1);
    kh->ke14 = (KmerEntry *)malloc(len * sizeof(KmerEntry));
    fread_fix(fp, len * sizeof(KmerEntry), kh->ke14);
#endif
    err_fclose(fp);
    return khash;
}

// KMER_LEN，
void gen_kmer_base(uint8_t **seq_arr, uint64_t *kmer_arr_size, int kmer_len)
{
    uint64_t i;
    uint64_t seq_up_val = (1 << (kmer_len << 1));
    *kmer_arr_size = (uint64_t)seq_up_val;
    *seq_arr = realloc(*seq_arr, seq_up_val * (uint64_t)kmer_len);
    for (i = 0; i < seq_up_val; ++i)
    {
        const uint64_t base_idx = i * kmer_len;
        for (int j = kmer_len - 1; j >= 0; --j)
        {
            (*seq_arr)[base_idx + kmer_len - 1 - j] = (i >> (j << 1)) & 3;
        }
    }
}

uint64_t global_num = 0;

// fmtseed，，
bwtintv_t fmt_search(FMTIndex *fmt, const uint8_t *q, int qlen)
{
    bwtintv_t ik;
    bwtintv_t ok1;
    bwtintv_t ok2;
    int i, c1, c2, x = 0;

    fmt_set_intv(fmt, q[x], ik);
    ik.info = x + 1;
    for (i = x + 1; i + 1 < qlen; i += 2)
    {
        if (q[i] < 4 && q[i + 1] < 4)
        {
            c1 = 3 - q[i];
            c2 = 3 - q[i + 1];
            fmt_extend2(fmt, &ik, &ok1, &ok2, 0, c1, c2);
            ik = ok2;
            ik.info = i + 1;
        }
        else
        {
            break;
        }
    }

    if (i < qlen && q[i] < 4)
    { // 
        c1 = 3 - q[i];
        fmt_extend1(fmt, &ik, &ok1, 0, c1);
        //if (qlen == 14) fprintf(stderr, "%ld %ld %ld\n", ok1.x[0], ok1.x[1], ok1.x[2]);
        //if (qlen == 14) ++global_num;
        //if (qlen == 14 && global_num % 10000 == 0) fprintf(stderr, "%ld\n", global_num);
        ik = ok1;
        ik.info = i + 1;
    }
    return ik;
}

// xmer hash
void fmt_create_kmer_index(FMTIndex *fmt) {
    uint64_t kmer_arr_size = 0;
    uint8_t *seq_arr = 0;
    gen_kmer_base(&seq_arr, &kmer_arr_size, 10);
    bwtintv_t ik;
    uint64_t j;
    int i, c1, c2;
    int qlen = 10;
    bwtint_t tk[4], tl[4];
    KmerHash *kh = &fmt->kmer_hash;

    kh->ke10 = (KmerEntryArr *)malloc((1 << (10 << 1)) * sizeof(KmerEntryArr));
    kh->ke11 = (KmerEntry *)malloc((1 << (11 << 1)) * sizeof(KmerEntry));
    kh->ke12 = (KmerEntry *)malloc((1 << (12 << 1)) * sizeof(KmerEntry));
    kh->ke13 = (KmerEntry *)malloc((1 << (13 << 1)) * sizeof(KmerEntry));
    kh->ke14 = (KmerEntry *)malloc((1 << (14 << 1)) * sizeof(KmerEntry));

    // 10kmer
    for (j = 0; j < kmer_arr_size; ++j)
    {
        uint8_t *q = &seq_arr[j * 10];
        uint8_t *mem_addr = kh->ke10[j].intv_arr;
        fmt_set_intv(fmt, q[0], ik);
        kmer_setval_at(mem_addr, ik, 0);

        // 
        for (i = 1; i + 1 < qlen; i += 2)
        {
            c1 = 3 - q[i];
            c2 = 3 - q[i + 1];

            fmt_e2_occ(fmt, ik.x[1] - 1, c1, c2, tk);
            fmt_e2_occ(fmt, ik.x[1] - 1 + ik.x[2], c1, c2, tl);
            // 
            ik.x[0] = ik.x[0] + (ik.x[1] <= fmt->primary && ik.x[1] + ik.x[2] - 1 >= fmt->primary) + tl[0] - tk[0];
            ik.x[1] = fmt->L2[c1] + 1 + tk[1];
            ik.x[2] = tl[1] - tk[1];
            kmer_setval_at(mem_addr, ik, i);

            // 
            ik.x[0] = ik.x[0] + (ik.x[1] <= fmt->primary && ik.x[1] + ik.x[2] - 1 >= fmt->primary) + tl[2] - tk[2];
            ik.x[1] = fmt->L2[c2] + 1 + tk[3];
            ik.x[2] = tl[3] - tk[3];
            kmer_setval_at(mem_addr, ik, i + 1);
        }
        { // 
            c1 = 3 - q[i];
            c2 = 3;
            fmt_e2_occ(fmt, ik.x[1] - 1, c1, c2, tk);
            fmt_e2_occ(fmt, ik.x[1] - 1 + ik.x[2], c1, c2, tl);
            // 
            ik.x[0] = ik.x[0] + (ik.x[1] <= fmt->primary && ik.x[1] + ik.x[2] - 1 >= fmt->primary) + tl[0] - tk[0];
            ik.x[1] = fmt->L2[c1] + 1 + tk[1];
            ik.x[2] = tl[1] - tk[1];
            kmer_setval_at(mem_addr, ik, i);
        }
    }
    // 11kmer
    gen_kmer_base(&seq_arr, &kmer_arr_size, 11);
    for (j = 0; j < kmer_arr_size; ++j)
    {
        bwtintv_t p = fmt_search(fmt, &seq_arr[j * 11], 11);
        kmer_setval_at(fmt->kmer_hash.ke11[j].intv_arr, p, 0);
    }

    // 12kmer
    gen_kmer_base(&seq_arr, &kmer_arr_size, 12);
    for (j = 0; j < kmer_arr_size; ++j)
    {
        bwtintv_t p = fmt_search(fmt, &seq_arr[j * 12], 12);
        kmer_setval_at(fmt->kmer_hash.ke12[j].intv_arr, p, 0);
    }
    gen_kmer_base(&seq_arr, &kmer_arr_size, 13);
    for (j = 0; j < kmer_arr_size; ++j)
    {
        bwtintv_t p = fmt_search(fmt, &seq_arr[j * 13], 13);
        kmer_setval_at(fmt->kmer_hash.ke13[j].intv_arr, p, 0);
    }

    // 14kmer
    gen_kmer_base(&seq_arr, &kmer_arr_size, 14);
    //fprintf(stderr, "14-size:%ld\n", kmer_arr_size);
    for (j = 0; j < kmer_arr_size; ++j)
    {
        //if (j % 10000 == 0) fprintf(stderr, "arr_size: %ld, %ld\n", j, kmer_arr_size);
        bwtintv_t p = fmt_search(fmt, &seq_arr[j * 14], 14);
        kmer_setval_at(fmt->kmer_hash.ke14[j].intv_arr, p, 0);
    }
    free(seq_arr);
}

// sa
void fmt_restore_sa(const char *fn, FMTIndex *fmt)
{
    char skipped[256];
    FILE *fp;
    bwtint_t primary;
    fp = xopen(fn, "rb");
    err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
    xassert(primary == fmt->primary, "SA-BWT inconsistency: primary is not the same.");
    err_fread_noeof(skipped, sizeof(bwtint_t), 4, fp); // skip
    err_fread_noeof(&fmt->sa_intv, sizeof(bwtint_t), 1, fp);
    err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
    xassert(primary == fmt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");
    fmt->n_sa = (fmt->seq_len + fmt->sa_intv) / fmt->sa_intv;
    fmt->sa = (uint8_t *)malloc(SA_BYTES(fmt->n_sa));
    fread_fix(fp, SA_BYTES(fmt->n_sa), fmt->sa);
    err_fclose(fp);
}

// interval-bwtfmt-index
FMTIndex *create_fmt_from_bwt(bwt_t *bwt)
{
    // FILE *fmt_out = fopen("fmt.txt", "w");
    FMTIndex *fmt = (FMTIndex *)calloc(1, sizeof(FMTIndex));
    fmt_gen_cnt_occ(fmt);

    bwtint_t i, j, k, m, n, n_occ, cnt[4], cnt2[4];
    uint32_t c[4], c2[16]; /*cbwt，c2pre-bwtbwt，AA..TT*/
    uint32_t *buf;         /* fmtbwt */

#ifdef FMT_MID_INTERVAL
    // check point
    uint32_t mc[4] = {0};
    uint32_t cnt_table[4][256]; // 4cnt_table，0,1,2,3,
    fmt_gen_cnt_table(cnt_table);
#endif

    fmt->seq_len = bwt->seq_len; // bwt，$，bwt matrix1
    for (i = 0; i < 5; ++i)
        fmt->L2[i] = bwt->L2[i]; // 
    fmt->primary = bwt->primary; // $，bwt matrix

    n_occ = (bwt->seq_len + FMT_OCC_INTERVAL - 1) / FMT_OCC_INTERVAL + 1; // check point 
    fmt->bwt_size = (fmt->seq_len * 2 + 15) >> 4;                         // 
    fmt->bwt_size += n_occ * 20;                                          // A,C,G,TAA,AC.....TG,TT20

#ifdef FMT_MID_INTERVAL
    uint32_t s1;
    bwtint_t mn_occ = (bwt->seq_len >> FMT_OCC_INTV_SHIFT) * (FMT_OCC_INTERVAL / FMT_MID_INTERVAL - 1);
    bwtint_t last_seq_len = bwt->seq_len % FMT_OCC_INTERVAL;
    mn_occ += (last_seq_len + FMT_MID_INTERVAL - 1) / FMT_MID_INTERVAL - 1;
    fmt->bwt_size += mn_occ * 4;
    i = 0;
#endif

    buf = (uint32_t *)calloc(fmt->bwt_size, 4); // fmt
    c[0] = c[1] = c[2] = c[3] = 0;
    // c2，ACGT，1occ
    for (i = 0; i < 4; ++i)
    {
        bwtint_t before_first_line = fmt->L2[i];
        bwt_occ4(bwt, before_first_line, cnt);
        for (j = i * 4, k = 0; k < 4; ++j, ++k)
            c2[j] = cnt[k];
    }
    // kbuf
    for (i = k = 0; i < bwt->seq_len; ++i)
    {
        // occ
        if (i % FMT_OCC_INTERVAL == 0)
        {
            memcpy(buf + k, c, sizeof(uint32_t) * 4); // bwt strocc
            k += 4;
            memcpy(buf + k, c2, sizeof(uint32_t) * 16); // pre-bwt:bwtocc
            k += 16;
#ifdef FMT_MID_INTERVAL
            mc[0] = mc[1] = mc[2] = mc[3] = 0;
#endif
        }
        // 328（pre-bwt）8(bwt)
        if (i % 16 == 0) // 3216，16，16
        {
            uint32_t pre_bwt_16_seq = 0;                   // 16pre-bwt
            uint32_t *bwt_addr = bwt_occ_intv(bwt, i) + 8;//4; // 48occuint32uint64，bwti，bwt-cp（check point）4uint32_t(8uint32_t)occ
            int offset = (i % OCC_INTERVAL) / 16;          // OCC_INTERVAL，16uint32，
            uint32_t bwt_16_seq = *(bwt_addr + offset);    // 16
            for (j = 0; j < 16; ++j)                       // bwt，
            {
                bwtint_t cur_str_line = i + j;   // bwt str
                if (cur_str_line < bwt->seq_len) // bwt str（bwt strbwt matrix1，bwt str$）
                {
                    uint8_t bwt_base = bwt_B0(bwt, cur_str_line); // bwt
                    // （bwt matrix）
                    bwtint_t cur_mtx_line = cur_str_line;
                    if (cur_str_line >= bwt->primary) // bwt$，，$，seq，bwt matrix
                        cur_mtx_line += 1;
                    bwt_occ4(bwt, cur_mtx_line, cnt); // bwt-checkpointocc
                    for (m = 0; m < 4; ++m)
                        c[m] = (uint32_t)cnt[m]; // mcur_bwt_mtx_line()，bwtocc

                    cnt[bwt_base] -= 1;                                                 // cur_bwt_mtx_line()，bwt_occ4(bwt, cur_bwt_mtx_line-1, cnt)
                    bwtint_t bwt_base_mtx_line = bwt->L2[bwt_base] + 1 + cnt[bwt_base]; // bwt_basebwt matrix（LF）

                    bwt_occ4(bwt, bwt_base_mtx_line, cnt2); // bwt_base_mtx_lineocc
                    for (n = 0; n < 4; ++n)
                    {
                        int c2_idx = bwt_base << 2 | n; // bwt base
                        c2[c2_idx] = (uint32_t)cnt2[n]; // pre-bwt:bwt 
                    }
                    bwtint_t bwt_base_str_line = bwt_base_mtx_line;         // bwt-str
                    if (bwt_base_str_line >= bwt->primary)                  // base_linebwt str，$，1
                        bwt_base_str_line -= 1;                             // bwt（$）
                    uint32_t pre_bwt_base = bwt_B0(bwt, bwt_base_str_line); // bwtpre-bwt
                    // ，bwt_basebwt matrix，$，bwt_base，
                    // pre_bwt_baseprimarybwt base，$，
                    if (bwt_base_mtx_line == bwt->primary)
                    {
                        // sec_bcp
                        fmt->sec_bcp = pre_bwt_base << 2 | bwt_base; // $A
                        fmt->sec_primary = cur_mtx_line;             // pre-bwt base$（bwt-matrix）
                        fmt->first_base = bwt_base;                  // 
                        fmt->last_base = pre_bwt_base;               // $（primarybwt base）
                    }
                    //  pre-bwt
                    pre_bwt_16_seq = pre_bwt_16_seq | (pre_bwt_base << (15 - j) * 2); // uint32_t
                }
                else
                    break;
            }
            // bwtpre_bwt
            uint32_t pre_and_bwt_seq = 0;
            uint32_t pre_and_bwt_seq_2 = 0;
            for (m = 0; m < 8; ++m)
            {
                int lshift_bit = 30 - 2 * m;
                pre_and_bwt_seq |= (((pre_bwt_16_seq & (3 << lshift_bit)) >> (m * 2)) | ((bwt_16_seq & (3 << lshift_bit)) >> ((m * 2) + 2)));
            }
            buf[k++] = pre_and_bwt_seq;

            if (j > 8)
            {
                for (m = 8; m > 0; --m)
                {
                    int lshift_bit = 2 * m - 2;
                    pre_and_bwt_seq_2 |= (((pre_bwt_16_seq & (3 << lshift_bit)) << (m * 2)) | ((bwt_16_seq & (3 << lshift_bit)) << (m * 2 - 2)));
                }

#ifdef FMT_MID_INTERVAL // 8+8mid interval occ
                s1 = pre_and_bwt_seq;
                for (m = 0; m < 4; ++m)
                    mc[m] += cnt_table[m][s1 & 0xff] + cnt_table[m][s1 >> 8 & 0xff] + cnt_table[m][s1 >> 16 & 0xff] + cnt_table[m][s1 >> 24];
#endif

#if FMT_MID_INTERVAL == 8 // mid interval8，
                for (m = 0; m < 4; ++m)
                    buf[k++] = mc[m];
#endif

                buf[k++] = pre_and_bwt_seq_2;

#ifdef FMT_MID_INTERVAL
                s1 = pre_and_bwt_seq_2;
                for (m = 0; m < 4; ++m)
                    mc[m] += cnt_table[m][s1 & 0xff] + cnt_table[m][s1 >> 8 & 0xff] + cnt_table[m][s1 >> 16 & 0xff] + cnt_table[m][s1 >> 24];

                if ((i + 16) % FMT_OCC_INTERVAL != 0 && j == 16 && ((i + 16) & FMT_MID_INTV_MASK) == 0)
                    for (m = 0; m < 4; ++m)
                        buf[k++] = mc[m];
#endif
            }
        }
    }
    // the last element
    memcpy(buf + k, c, sizeof(uint32_t) * 4);
    k += 4;
    memcpy(buf + k, c2, sizeof(uint32_t) * 16);
    k += 16;
    xassert(k == fmt->bwt_size, "inconsistent fmt_size");
    // update fmt
    fmt->bwt = buf;
    return fmt;
}

// ，bwt basebpre-bwt strocc
inline void fmt_e2_occ(const FMTIndex *fmt, bwtint_t k, int b1, int b2, bwtint_t cnt[4])
{
    uint32_t x = 0;
    uint32_t *p, *q, tmp;
    bwtint_t str_line = k, cp_line = k & (~FMT_OCC_INTV_MASK); // cp = check point
    int i, ti = b1 << 2 | b2;
    cnt[0] = 0;
    cnt[1] = 0;
    cnt[2] = 0;
    if (k == (bwtint_t)(-1))
    {
        p = fmt->bwt + 4 + b1 * 4;
        for (i = b2 + 1; i < 4; ++i)
            cnt[2] += p[i];
        cnt[3] = p[b2];
        return;
    }
    k -= (k >= fmt->primary); // kbwtbwt（$，$，1）
    p = fmt_occ_intv(fmt, k);
    // fprintf(stderr, "k: %ld\n", k);
    for (i = b1 + 1; i < 4; ++i)
        cnt[0] += p[i]; // b1occ
    cnt[1] = p[b1];     // b1occ
    q = p + 4 + b1 * 4;
    for (i = b2 + 1; i < 4; ++i)
        cnt[2] += q[i]; // b2occ
    cnt[3] = q[b2];     // b2occ
    p += 20;

    // mid interval
    int mk = k % FMT_OCC_INTERVAL;
    int n_mintv = mk >> FMT_MID_INTV_SHIFT;
    if (n_mintv > 0) // mid interval
    {
        p += n_mintv * (4 + (FMT_MID_INTERVAL >> 3)) - 4; // mid interval check point，A C G T
        q = p + b1;
        for (i = b1 + 1; i < 4; ++i)
            x += p[i]; // b1occ
        cnt[0] += __fmt_mid_sum(x);
        x = *q;
        cnt[1] += __fmt_mid_sum(x); // b1occ
        for (i = 3; i > b2; --i)
            cnt[2] += x >> (i << 3) & 0xff; // b2occ
        cnt[3] += x >> (b2 << 3) & 0xff;    // b2occ
        x = 0;
        p += 4;
    }
    uint32_t *end = p + ((k >> 3) - ((k & ~FMT_MID_INTV_MASK) >> 3));
    for (; p < end; ++p)
    {
        x += __fmt_occ_e2_aux2(fmt, ti, *p);
    }

    tmp = *p & ~((1U << ((~k & 7) << 2)) - 1);
    x += __fmt_occ_e2_aux2(fmt, ti, tmp);

    if (b1 == 0)
    {
        x -= (~k & 7) << 8;
        if (b2 == 0)
            x -= (~k & 7) << 24;
    }
    // second_primary,
    if (b1 == fmt->first_base && cp_line < fmt->sec_primary && str_line >= fmt->sec_primary)
    {
        if (b2 < fmt->last_base)
            cnt[2] -= 1;
        else if (b2 == fmt->last_base)
            cnt[3] -= 1;
    }
    cnt[0] += x & 0xff;
    cnt[1] += x >> 8 & 0xff;
    cnt[2] += x >> 16 & 0xff;
    cnt[3] += x >> 24 & 0xff;
}

// ，bwt basebpre-bwt strocc
inline void fmt_direct_e2_occ(const FMTIndex *fmt, bwtint_t k, int b1, int b2, bwtint_t cnt[4])
{
    uint32_t x = 0;
    uint32_t *p, *q, tmp;
    bwtint_t str_line = k, cp_line = k & (~FMT_OCC_INTV_MASK); // cp = check point
    int ti = b1 << 2 | b2;
    if (k == (bwtint_t)(-1))
    {
        p = fmt->bwt + 4 + b1 * 4;
        cnt[3] = p[b2];
        return;
    }
    k -= (k >= fmt->primary); // kbwtbwt（$，$，1）
    p = fmt_occ_intv(fmt, k);
    q = p + 4 + b1 * 4;
    cnt[3] = q[b2]; // b2occ
    p += 20;

    // mid interval
    int mk = k % FMT_OCC_INTERVAL;
    int n_mintv = mk >> FMT_MID_INTV_SHIFT;
    if (n_mintv > 0) // mid interval
    {
        p += n_mintv * (4 + (FMT_MID_INTERVAL >> 3)) - 4; // mid interval check point，A C G T
        q = p + b1;
        x = *q;
        cnt[3] += x >> (b2 << 3) & 0xff; // b2occ
        x = 0;
        p += 4;
    }
    uint32_t *end = p + ((k >> 3) - ((k & ~FMT_MID_INTV_MASK) >> 3));
    for (; p < end; ++p)
    {
        x += __fmt_occ_e2_aux2(fmt, ti, *p);
    }
    tmp = *p & ~((1U << ((~k & 7) << 2)) - 1);
    x += __fmt_occ_e2_aux2(fmt, ti, tmp);
    if (b1 == 0)
    {
        x -= (~k & 7) << 8;
        if (b2 == 0)
            x -= (~k & 7) << 24;
    }
    // second_primary,
    if (b1 == fmt->first_base && cp_line < fmt->sec_primary && str_line >= fmt->sec_primary && b2 == fmt->last_base)
        cnt[3] -= 1;
    cnt[3] += x >> 24 & 0xff;
}

inline void fmt_direct2_e2_occ(const FMTIndex *fmt, bwtint_t k, int b1, int b2, bwtint_t cnt[4])
{
    uint32_t x = 0;
    uint32_t *p, *q, tmp;
    bwtint_t str_line = k, cp_line = k & (~FMT_OCC_INTV_MASK); // cp = check point
    int ti = b1 << 2 | b2;
    cnt[1] = 0;
    if (k == (bwtint_t)(-1))
    {
        p = fmt->bwt + 4 + b1 * 4;
        cnt[3] = p[b2];
        return;
    }
    k -= (k >= fmt->primary); // kbwtbwt（$，$，1）
    p = fmt_occ_intv(fmt, k);
    cnt[1] = p[b1]; // b1occ
    q = p + 4 + b1 * 4;
    cnt[3] = q[b2]; // b2occ
    p += 20;

    // mid interval
    int mk = k % FMT_OCC_INTERVAL;
    int n_mintv = mk >> FMT_MID_INTV_SHIFT;
    if (n_mintv > 0) // mid interval
    {
        p += n_mintv * (4 + (FMT_MID_INTERVAL >> 3)) - 4; // mid interval check point，A C G T
        q = p + b1;
        x = *q;
        cnt[1] += __fmt_mid_sum(x);      // b1occ
        cnt[3] += x >> (b2 << 3) & 0xff; // b2occ
        x = 0;
        p += 4;
    }
    uint32_t *end = p + ((k >> 3) - ((k & ~FMT_MID_INTV_MASK) >> 3));
    for (; p < end; ++p)
    {
        x += __fmt_occ_e2_aux2(fmt, ti, *p);
    }
    tmp = *p & ~((1U << ((~k & 7) << 2)) - 1);
    x += __fmt_occ_e2_aux2(fmt, ti, tmp);
    if (b1 == 0)
    {
        x -= (~k & 7) << 8;
        if (b2 == 0)
            x -= (~k & 7) << 24;
    }
    // second_primary,
    if (b1 == fmt->first_base && cp_line < fmt->sec_primary && str_line >= fmt->sec_primary && b2 == fmt->last_base)
        cnt[3] -= 1;
    cnt[1] += x >> 8 & 0xff;
    cnt[3] += x >> 24 & 0xff;
}

// 
inline void fmt_extend2(const FMTIndex *fmt, bwtintv_t *ik, bwtintv_t *ok1, bwtintv_t *ok2, int is_back, int b1, int b2)
{
    bwtint_t tk[4], tl[4];
    bwtintv_t intv = {0};
    // tkk，tll
    fmt_e2_occ(fmt, ik->x[!is_back] - 1, b1, b2, tk);
    fmt_e2_occ(fmt, ik->x[!is_back] - 1 + ik->x[2], b1, b2, tl);
    // 
    intv.x[!is_back] = fmt->L2[b1] + 1 + tk[1];
    intv.x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= fmt->primary && ik->x[!is_back] + ik->x[2] - 1 >= fmt->primary) + tl[0] - tk[0];
    intv.x[2] = tl[1] - tk[1];
    *ok1 = intv;
    // 
    intv.x[is_back] = intv.x[is_back] + (intv.x[!is_back] <= fmt->primary && intv.x[!is_back] + intv.x[2] - 1 >= fmt->primary) + tl[2] - tk[2];
    intv.x[!is_back] = fmt->L2[b2] + 1 + tk[3];
    intv.x[2] = tl[3] - tk[3];
    *ok2 = intv;
}

// 
inline void fmt_direct_extend2(const FMTIndex *fmt, bwtintv_t *ik, bwtintv_t *ok1, bwtintv_t *ok2, int is_back, int b1, int b2)
{
    bwtint_t tk[4], tl[4];
    // tkk，tll
    fmt_direct_e2_occ(fmt, ik->x[!is_back] - 1, b1, b2, tk);
    fmt_direct_e2_occ(fmt, ik->x[!is_back] - 1 + ik->x[2], b1, b2, tl);
    ok2->x[!is_back] = fmt->L2[b2] + 1 + tk[3];
    ok2->x[2] = tl[3] - tk[3];
    ok2->x[is_back] = 0;
}

// 
inline void fmt_direct2_extend2(const FMTIndex *fmt, bwtintv_t *ik, bwtintv_t *ok1, bwtintv_t *ok2, int is_back, int b1, int b2)
{
    bwtint_t tk[4], tl[4];
    // tkk，tll
    fmt_direct2_e2_occ(fmt, ik->x[!is_back] - 1, b1, b2, tk);
    fmt_direct2_e2_occ(fmt, ik->x[!is_back] - 1 + ik->x[2], b1, b2, tl);
    ok1->x[!is_back] = fmt->L2[b1] + 1 + tk[1];
    ok1->x[2] = tl[1] - tk[1];
    ok1->x[is_back] = 0;
    ok2->x[!is_back] = fmt->L2[b2] + 1 + tk[3];
    ok2->x[2] = tl[3] - tk[3];
    ok2->x[is_back] = 0;
}

// 
inline void fmt_extend1(const FMTIndex *fmt, bwtintv_t *ik, bwtintv_t *ok, int is_back, int b1)
{
    bwtint_t tk[4], tl[4];
    int b2 = 3; // ，T，，b2
    // tkk，tll
    fmt_e2_occ(fmt, ik->x[!is_back] - 1, b1, b2, tk);
    fmt_e2_occ(fmt, ik->x[!is_back] - 1 + ik->x[2], b1, b2, tl);
    // 
    ok->x[!is_back] = fmt->L2[b1] + 1 + tk[1];
    ok->x[2] = tl[1] - tk[1];
    // 
    ok->x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= fmt->primary && ik->x[!is_back] + ik->x[2] - 1 >= fmt->primary) + tl[0] - tk[0];
}
inline void fmt_direct_extend1(const FMTIndex *fmt, bwtintv_t *ik, bwtintv_t *ok, int is_back, int b1)
{
    bwtint_t tk[4], tl[4];
    int b2 = 3; // ，T，，b2
    // tkk，tll
    fmt_direct2_e2_occ(fmt, ik->x[!is_back] - 1, b1, b2, tk);
    fmt_direct2_e2_occ(fmt, ik->x[!is_back] - 1 + ik->x[2], b1, b2, tl);
    // 
    ok->x[!is_back] = fmt->L2[b1] + 1 + tk[1];
    ok->x[2] = tl[1] - tk[1];
}

// 
static void direct_extend(const FMTIndex *fmt, int len, const uint8_t *q, int left_pos, int right_pos, bwtint_t mtx_line, bwtint_t ref_pos, bwtintv_t *mt, int mtidx)
{
#define PAC_BASE(pac, l) ((pac)[(l) >> 2] >> ((~(l) & 3) << 1) & 3)
#define EXTEND_BASE_LOOP(qcond, rcond, qstep, rstep) \
    while (k != qcond && r != rcond)                 \
    {                                                \
        const int base = PAC_BASE(fmt->pac, r);      \
        if (q[k] != base)                            \
            break;                                   \
        k += qstep;                                  \
        r += rstep;                                  \
    }
#define EXTEND_BASE_LOOP_COMP(qcond, rcond, qstep, rstep) \
    while (k != qcond && r != rcond)                      \
    {                                                     \
        const int base = 3 - PAC_BASE(fmt->pac, r);       \
        if (q[k] != base)                                 \
            break;                                        \
        k += qstep;                                       \
        r += rstep;                                       \
    }

    int k;
    int64_t r, rp;
    mt->num_match = 1;
    // rp = fmt_sa(fmt, mtx_line);
    rp = ref_pos;
    r = rp >= fmt->l_pac ? (fmt->l_pac << 1) - 1 - rp : rp;
    k = right_pos;
    if (rp < fmt->l_pac) // 
    {
        // 
        r += right_pos - left_pos;
        EXTEND_BASE_LOOP(len, fmt->l_pac, 1, 1);
        mt->rm[mtidx].qe = k;
        mt->rm[mtidx].reverse = 0;
        // ，x
        r -= k - left_pos + 1;
        k = left_pos - 1;
        EXTEND_BASE_LOOP(-1, -1, -1, -1);
        mt->rm[mtidx].qs = k + 1;
        mt->rm[mtidx].rs = r + 1;
    }
    else // 
    {
        r -= right_pos - left_pos;
        EXTEND_BASE_LOOP_COMP(len, -1, 1, -1);
        mt->rm[mtidx].qe = k;
        mt->rm[mtidx].reverse = 1;
        // x
        r += k - left_pos + 1;
        k = left_pos - 1;
        EXTEND_BASE_LOOP_COMP(-1, fmt->l_pac, -1, 1);
        mt->rm[mtidx].qs = k + 1;
        mt->rm[mtidx].rs = r - 1;
    }
    mt->info = mt->rm[mtidx].qs;
    mt->info = mt->info << 32 | mt->rm[mtidx].qe;
    mt->x[2] = 1;
}

static inline void fmt_reverse_intvs(bwtintv_v *p)
{
    if (p->n > 1)
    {
        int j;
        for (j = 0; j < p->n >> 1; ++j)
        {
            bwtintv_t tmp = p->a[p->n - 1 - j];
            p->a[p->n - 1 - j] = p->a[j];
            p->a[j] = tmp;
        }
    }
}

int fmt_smem_forward(const FMTIndex *fmt, int len, const uint8_t *q, int x, int min_intv, int min_seed_len, bwtintv_v *mem)
{
    int i, ret, kmer_len;
    bwtintv_t ik = {0}, ok1 = {0}, ok2 = {0};
    bwtintv_t mt = {0};
    uint32_t qbit = 0;
    mem->n = 0;
    if (q[x] > 3) return x + 1;

    if (min_intv < 1) min_intv = 1; // the interval size should be at least 1

    qbit = build_forward_kmer(&q[x], len - x, HASH_KMER_LEN, &kmer_len);
    bwt_kmer_get(&fmt->kmer_hash, &ik, qbit, kmer_len-1); // 
    ik.info = x + kmer_len;

// check change of the interval size and whether the interval size is too small to be extended further
#define PUSH_VAL_AND_SKIP_FORWARD(iv)         \
    do                                \
    {                                 \
        kv_push(bwtintv_t, *mem, iv); \
        goto fmt_smem_forward_end;    \
    } while (0)

    if (kmer_len != HASH_KMER_LEN) // N
        PUSH_VAL_AND_SKIP_FORWARD(ik);

    // kmer
    for (i = (int)ik.info; i + 1 < len; i += 2)
    { // forward search
        if (q[i] < 4 && q[i + 1] < 4)
        {
            fmt_extend2(fmt, &ik, &ok1, &ok2, 0, 3 - q[i], 3 - q[i + 1]);
            __builtin_prefetch(fmt_occ_intv(fmt, ok2.x[1] - 1), 0, 2);
            __builtin_prefetch(fmt_occ_intv(fmt, ok2.x[1] - 1 + ok2.x[2]), 0, 2);

#if 1
            if (min_intv == 1 && ok2.x[2] == min_intv)
            {
                direct_extend(fmt, len, q, x, i + 2, ok2.x[0], fmt_sa(fmt, ok2.x[0]), &mt, 0);
                kv_push(bwtintv_t, *mem, mt);
                ret = (uint32_t)mt.info;
                goto fmt_smem_forward_end;
            }
#endif
            if (ok2.x[2] < min_intv) {
                if (ok1.x[2] < min_intv) {
                    PUSH_VAL_AND_SKIP_FORWARD(ik);
                } else {
                    ok1.info = i + 1;
                    PUSH_VAL_AND_SKIP_FORWARD(ok1);
                }
            } else {
                ik = ok2;
                ik.info = i + 2;
            }
        }
        else if (q[i] < 4) // q[i+1] >= 4
        {
            fmt_extend1(fmt, &ik, &ok1, 0, 3 - q[i]);
            if (ok1.x[2] < min_intv) {
                PUSH_VAL_AND_SKIP_FORWARD(ik);
            } else {
                ok1.info = i + 1;
                PUSH_VAL_AND_SKIP_FORWARD(ok1);
            }
        }
        else // q[i] >= 4
        {
            PUSH_VAL_AND_SKIP_FORWARD(ik);
        }
    }
    for (; i == len - 1; ++i) // 
    {
        if (q[i] < 4)
        {
            fmt_extend1(fmt, &ik, &ok1, 0, 3 - q[i]);
            if (ok1.x[2] < min_intv) {
                PUSH_VAL_AND_SKIP_FORWARD(ik);
            } else {
                ok1.info = i + 1;
                PUSH_VAL_AND_SKIP_FORWARD(ok1);
            }
        }
        else
            PUSH_VAL_AND_SKIP_FORWARD(ik);
    }

fmt_smem_forward_end:
    if (mem->n == 0)
        kv_push(bwtintv_t, *mem, ik);
    ret = mem->a[0].info;
    mem->a[0].info |= (uint64_t)(x) << 32;
    if (mt.num_match == 0)
        ret = (uint32_t)mem->a[0].info; // this will be the returned value，
    return ret;
}

// smem（seed）(lm: long_smem)
int fmt_smem(const FMTIndex *fmt, int len, const uint8_t *q, int x, int min_intv, int min_seed_len, bwtintv_t *lm, bwtintv_v *mem, bwtintv_v *tmpvec)
{
    if (x == 0 || q[x - 1] > 3)
        return fmt_smem_forward(fmt, len, q, x, min_intv, min_seed_len, mem);

    // int flag = 0;
    int i, j, ret, kmer_len;
    bwtintv_t ik = {0}, ok1 = {0}, ok2 = {0};
    bwtintv_t mt = {0};
    bwtintv_v *curr;
    uint32_t qbit = 0;
    mem->n = 0;
    if (q[x] > 3) return x + 1;

    if (min_intv < 1) min_intv = 1; // the interval size should be at least 1
    curr = tmpvec;    // use the temporary vector if provided

    qbit = build_forward_kmer(&q[x], len - x, HASH_KMER_LEN, &kmer_len);
    bwt_kmer_get(&fmt->kmer_hash, &ik, qbit, 0); // 
    ik.info = x + 1;
// check change of the interval size and whether the interval size is too small to be extended further
#define CHECK_INTV_CHANGE(iv, ov, end_pos) \
    if (ov.x[2] != iv.x[2])                \
    {                                      \
        kv_push(bwtintv_t, *curr, iv);     \
        if (ov.x[2] < min_intv)            \
            break;                         \
    }                                      \
    iv = ov;                               \
    iv.info = end_pos

#define PUSH_VAL_AND_SKIP(iv)          \
    do                                 \
    {                                  \
        kv_push(bwtintv_t, *curr, iv); \
        goto backward_search;          \
    } while (0)

    // kmer
    for (curr->n = 0, j = 1; j < kmer_len; ++j)
    {
        bwt_kmer_get(&fmt->kmer_hash, &ok1, qbit, j);
        CHECK_INTV_CHANGE(ik, ok1, x + j + 1);
    }
    if (kmer_len != HASH_KMER_LEN) // N
        PUSH_VAL_AND_SKIP(ik);

    // kmer
    for (i = (int)ik.info; i + 1 < len; i += 2)
    { // forward search
        if (q[i] < 4 && q[i + 1] < 4)
        {
            fmt_extend2(fmt, &ik, &ok1, &ok2, 0, 3 - q[i], 3 - q[i + 1]);
            __builtin_prefetch(fmt_occ_intv(fmt, ok2.x[1] - 1), 0, 2);
            __builtin_prefetch(fmt_occ_intv(fmt, ok2.x[1] - 1 + ok2.x[2]), 0, 2);
            CHECK_INTV_CHANGE(ik, ok1, i + 1);
            CHECK_INTV_CHANGE(ik, ok2, i + 2);
#if 1
            // 
            if (min_intv == 1 && ok2.x[2] == min_intv)
            {
                direct_extend(fmt, len, q, x, i + 2, ok2.x[0], fmt_sa(fmt, ok2.x[0]), & mt, 0);
                kv_push(bwtintv_t, *mem, mt); // min_seed_len
                ret = (uint32_t)mt.info;
                if (mt.rm[0].qs == 0 || q[mt.rm[0].qs - 1] > 3)
                    goto fmt_smem_end;
                goto backward_search;
            }
#endif
        }
        else if (q[i] < 4) // q[i+1] >= 4
        {
            fmt_extend1(fmt, &ik, &ok1, 0, 3 - q[i]);
            CHECK_INTV_CHANGE(ik, ok1, i + 1);
            PUSH_VAL_AND_SKIP(ik);
        }
        else // q[i] >= 4
        {
            PUSH_VAL_AND_SKIP(ik);
        }
    }
    for (; i == len - 1; ++i) // 
    {
        if (q[i] < 4)
        {
            fmt_extend1(fmt, &ik, &ok1, 0, 3 - q[i]);
            CHECK_INTV_CHANGE(ik, ok1, i + 1);
        }
        else
            PUSH_VAL_AND_SKIP(ik);
    }
    if (i == len)
        kv_push(bwtintv_t, *curr, ik); // push the last interval if we reach the end

backward_search:
    fmt_reverse_intvs(curr); // s.t. smaller intervals (i.e. longer matches) visited first
    if (mt.num_match == 0)
        ret = curr->a[0].info; // this will be the returned value，
    else
        ret = (uint32_t)mt.info;
        // ，
#define CHECK_ADD_MEM(pos, intv, mem)                                                                           \
    if (((int)((intv).info) - (pos) >= min_seed_len) && (mem->n == 0 || (pos) < mem->a[mem->n - 1].info >> 32)) \
    {                                                                                                           \
        (intv).info |= (uint64_t)(pos) << 32;                                                                   \
        kv_push(bwtintv_t, *mem, (intv));                                                                       \
    }

#define CHECK_INTV_ADD_MEM(ok, pos, intv, mem) \
    if (ok.x[2] < min_intv)                    \
    {                                          \
        CHECK_ADD_MEM(pos, intv, mem);         \
        break;                                 \
    }

    int last_kmer_start = 0;
    for (j = 0; j < curr->n; ++j)
    {
        bwtintv_t *p = &curr->a[j]; // 
        if (p->info - x < HASH_KMER_LEN)
        {
            if (last_kmer_start && kmer_len == HASH_KMER_LEN && p->info == last_kmer_start && p->info - kmer_len > 0 && q[p->info - kmer_len] < 4)
                qbit = ((qbit << 2) | (3 - q[p->info - kmer_len])) & ((1L << (kmer_len << 1)) - 1); // kmer
            else
                qbit = build_backward_kmer(q, p->info - 1, HASH_KMER_LEN, &kmer_len); // kmer
            last_kmer_start = p->info - 1;
            i = 1;
            do
            {
                bwt_kmer_get(&fmt->kmer_hash, &ik, qbit, kmer_len - i++);
            } while (ik.x[2] < min_intv);
            if (i > 2)
                continue;
            p->x[0] = ik.x[1];
            p->x[1] = ik.x[0];
            p->x[2] = ik.x[2];
            i = p->info - (kmer_len - i + 3);
        }
        else
        {
            i = x - 1;
        }
        for (; i > 0; i -= 2)
        {
            if (q[i] < 4 && q[i - 1] < 4) // 
            {
                fmt_direct2_extend2(fmt, p, &ok1, &ok2, 1, q[i], q[i - 1]);
                __builtin_prefetch(fmt_occ_intv(fmt, ok2.x[0] - 1), 0, 2);
                __builtin_prefetch(fmt_occ_intv(fmt, ok2.x[0] - 1 + ok2.x[2]), 0, 2);
                CHECK_INTV_ADD_MEM(ok1, i + 1, *p, mem);
                ok1.info = p->info;
                CHECK_INTV_ADD_MEM(ok2, i, ok1, mem);
                ok2.info = p->info;
                *p = ok2;
            }
            else if (q[i] < 4) // 
            {
                fmt_direct_extend1(fmt, p, &ok1, 1, q[i]);
                CHECK_INTV_ADD_MEM(ok1, i + 1, *p, mem);
                ok1.info = p->info;
                CHECK_ADD_MEM(i, ok1, mem);
                goto fmt_smem_end;
            }
            else
            { // 
                CHECK_ADD_MEM(i + 1, *p, mem);
                goto fmt_smem_end;
            }
        }
        for (; i == 0; --i)
        { // 
            if (q[i] < 4)
            {
                fmt_direct_extend1(fmt, p, &ok1, 1, q[i]);
                CHECK_INTV_ADD_MEM(ok1, i + 1, *p, mem);
                ok1.info = p->info;
                *p = ok1;
            }
            else
            {
                CHECK_ADD_MEM(i + 1, *p, mem);
                goto fmt_smem_end;
            }
        }
        if (i == -1)
        {
            CHECK_ADD_MEM(i + 1, *p, mem);
            goto fmt_smem_end;
        }
    }

fmt_smem_end:
    fmt_reverse_intvs(mem); // s.t. sorted by the start coordinate
    return ret;
}

int fmt_seed_strategy1(const FMTIndex *fmt, int len, const uint8_t *q, int x, int min_len, int max_intv, bwtintv_t *mem)
{
    int i, kmer_len, first_extend_len;
    bwtintv_t ik = {0}, ok1 = {0}, ok2 = {0};
    uint64_t qbit;
    memset(mem, 0, sizeof(bwtintv_t));
    if (q[x] > 3)
        return x + 1;
    if (len - x <= min_len)
        return len;

    qbit = build_forward_kmer(&q[x], len - x, HASH_KMER_LEN, &kmer_len);
    bwt_kmer_get(&fmt->kmer_hash, &ik, qbit, kmer_len - 1); // 
    ik.info = x + kmer_len;

#define COND_SET_RETURN(iv, ov, start_pos, end_pos, max_intv, min_len) \
    if (iv.x[2] < max_intv && end_pos - start_pos >= min_len)          \
    {                                                                  \
        (ov) = (iv);                                                   \
        (ov).info = (uint64_t)start_pos << 32 | (end_pos + 1);         \
        return end_pos + 1;                                            \
    }
#if 1

    first_extend_len = x + min_len + 1;
    first_extend_len = MIN(len, first_extend_len);
    for (i = (int)ik.info; i + 1 < first_extend_len; i += 2)
    {
        if (q[i] < 4 && q[i + 1] < 4)
        {
            fmt_extend2(fmt, &ik, &ok1, &ok2, 0, 3 - q[i], 3 - q[i + 1]);
            __builtin_prefetch(fmt_occ_intv(fmt, ok2.x[1] - 1), 0, 2);
            __builtin_prefetch(fmt_occ_intv(fmt, ok2.x[1] - 1 + ok2.x[2]), 0, 2);
            ik = ok2;
        }
        else if (q[i] < 4)
            return i + 2;
        else
            return i + 1;
    }
    COND_SET_RETURN(ik, *mem, x, i - 1, max_intv, min_len);
    for (; i + 1 < len; i += 2)
#else
    for (i = (int)ik.info; i + 1 < len; i += 2)
#endif
    { // forward search
        if (q[i] < 4 && q[i + 1] < 4)
        {
            fmt_extend2(fmt, &ik, &ok1, &ok2, 0, 3 - q[i], 3 - q[i + 1]);
            __builtin_prefetch(fmt_occ_intv(fmt, ok2.x[1] - 1), 0, 2);
            __builtin_prefetch(fmt_occ_intv(fmt, ok2.x[1] - 1 + ok2.x[2]), 0, 2);
            COND_SET_RETURN(ok1, *mem, x, i, max_intv, min_len);
            COND_SET_RETURN(ok2, *mem, x, i + 1, max_intv, min_len);
            ik = ok2;
        }
        else if (q[i] < 4) // q[i+1] >= 4
        {
            fmt_extend1(fmt, &ik, &ok1, 0, 3 - q[i]);
            COND_SET_RETURN(ok1, *mem, x, i, max_intv, min_len);
            return i + 2;
        }
        else // q[i] >= 4
        {
            return i + 1;
        }
    }
    if (i == len - 1 && q[i] < 4)
    {
        fmt_extend1(fmt, &ik, &ok1, 0, 3 - q[i]);
        COND_SET_RETURN(ok1, *mem, x, i, max_intv, min_len);
    }
    return len;
}

// kbwt str
inline static void fmt_get_previous_base(const FMTIndex *fmt, bwtint_t k, uint8_t *b1, uint8_t *b2)
{
    uint32_t *p;
    uint8_t base2;
    // ，check point
    p = fmt_occ_intv(fmt, k); // check point
    p += 20;                  // bwt
    // ，mid check point
    int mk = k & FMT_OCC_INTV_MASK;
    int n_mintv = mk >> FMT_MID_INTV_SHIFT;
    p += n_mintv * (4 + (FMT_MID_INTERVAL >> 3)); // midbwt
    // ，uint32_t
    p += (k & FMT_MID_INTV_MASK) >> 3; // uint32_t8（8bwt）
    // ，
    base2 = *p >> ((~(k) & 0x7) << 2) & 0xf;
    *b2 = base2 >> 2 & 3;
    *b1 = base2 & 3;
}

inline static void fmt_previous_line(const FMTIndex *fmt, bwtint_t k, bwtint_t *k1, bwtint_t *k2)
{
    uint8_t b1, b2;
    uint32_t x = 0;
    bwtint_t cnt[4];
    k = k - (k >= fmt->primary); // kbwtbwt（$，$，1）
    uint32_t *p, *pocc, *ptocc, tmp;
    uint8_t base2;
    bwtint_t str_line = k, cp_line = k & (~FMT_OCC_INTV_MASK);
    // ，check point
    pocc = fmt_occ_intv(fmt, k); // check point
    p = pocc + 20;               // bwt
    // ，mid check point
    int mk = k & FMT_OCC_INTV_MASK;
    int n_mintv = mk >> FMT_MID_INTV_SHIFT;
    p += n_mintv * (4 + (FMT_MID_INTERVAL >> 3)); // midbwt
    ptocc = p;
    // ，uint32_t
    p += (k & FMT_MID_INTV_MASK) >> 3; // uint32_t8（8bwt）
    // ，
    base2 = *p >> ((~(k) & 0x7) << 2) & 0xf;
    b2 = base2 >> 2 & 3;
    b1 = base2 & 3;

    cnt[1] = pocc[b1];
    cnt[3] = (pocc + 4 + b1 * 4)[b2];
    if (n_mintv > 0) {
        ptocc -= 4;
        x = *(ptocc + b1);
        cnt[1] += __fmt_mid_sum(x);
        cnt[3] += x >> (b2 << 3) & 0xff;
        x = 0;
        ptocc += 4;
    }

    uint32_t *end = ptocc + ((k >> 3) - ((k & ~FMT_MID_INTV_MASK) >> 3));
    int ti = b1 << 2 | b2;
    for (; ptocc < end; ++ptocc)
    {
        x += __fmt_occ_e2_aux2(fmt, ti, *ptocc);
    }

    tmp = *ptocc & ~((1U << ((~k & 7) << 2)) - 1);
    x += __fmt_occ_e2_aux2(fmt, ti, tmp);

    if (b1 == 0)
    {
        x -= (~k & 7) << 8;
        if (b2 == 0)
            x -= (~k & 7) << 24;
    }
    if (b1 == fmt->first_base && b2 == fmt->last_base && cp_line < fmt->sec_primary && str_line >= fmt->sec_primary)
    {
        cnt[3] -= 1;
    }
    cnt[1] += x >> 8 & 0xff;
    cnt[3] += x >> 24 & 0xff;

    *k1 = fmt->L2[b1] + cnt[1];
    *k2 = fmt->L2[b2] + cnt[3];
}

bwtint_t fmt_sa(const FMTIndex *fmt, bwtint_t k)
{
    bwtint_t sa = 0, mask = fmt->sa_intv - 1;
    bwtint_t k1, k2;
    while (k & mask)
    {
        ++sa;
        fmt_previous_line(fmt, k, &k1, &k2);
        __builtin_prefetch(fmt_occ_intv(fmt, k2), 0, 2);
        if (!(k1 & mask))
        {
            k = k1;
            break;
        }
        ++sa;
        if (!(k2 & mask))
        {
            k = k2;
            break;
        }
        k = k2;
    }
    sa += bwt_get_sa(fmt->sa, k / fmt->sa_intv);
    return sa;
}

bwtint_t fmt_sa_offset(const FMTIndex *fmt, bwtint_t k)
{
    bwtint_t sa = 0, mask = fmt->sa_intv - 1;
    bwtint_t k1, k2;
    while (k & mask)
    {
        ++sa;
        fmt_previous_line(fmt, k, &k1, &k2);
        __builtin_prefetch(fmt_occ_intv(fmt, k2), 0, 2);
        if (!(k1 & mask))
        {
            k = k1;
            break;
        }
        ++sa;
        if (!(k2 & mask))
        {
            k = k2;
            break;
        }
        k = k2;
    }
    sa = (sa << 48) | (k / fmt->sa_intv);
    return sa;
}