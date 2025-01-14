#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <immintrin.h>
#include <emmintrin.h>
#include "utils.h"
#include "debug.h"
#include "profiling.h"

#define ELIMINATE_DIFF_1
// #define ELIMINATE_DIFF_3

#define NO_VAL -1

#define SIMD_WIDTH 16

extern int ksw_extend2_avx2_u8(int qlen, const uint8_t *query,  int tlen, const uint8_t *target, int is_left, int m,  const int8_t *mat,  int o_del,  int e_del,
		int o_ins, int e_ins, int a, int b, int w, int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore, int *_max_off, buf_t *buf);

int ksw_extend2_origin(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int is_left, int m, const int8_t *mat, int o_del, int e_del,
					   int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore, int *_max_off);

static const uint16_t h_vec_int_mask[SIMD_WIDTH][SIMD_WIDTH] = {
	{0xffff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xffff, 0xffff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xffff, 0xffff, 0xffff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0, 0, 0, 0, 0},
	{0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0, 0, 0, 0},
	{0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0, 0, 0},
	{0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0, 0},
	{0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0},
	{0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0},
	{0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0},
	{0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff}
};

#define permute_mask _MM_SHUFFLE(0, 1, 2, 3)

// 
#define SIMD_INIT \
	int oe_del = o_del + e_del, oe_ins = o_ins + e_ins; \
	__m256i zero_vec; \
	__m256i max_vec; \
	__m256i oe_del_vec; \
	__m256i oe_ins_vec; \
	__m256i e_del_vec; \
	__m256i e_ins_vec; \
	__m256i h_vec_mask[SIMD_WIDTH]; \
	zero_vec = _mm256_setzero_si256(); \
	oe_del_vec = _mm256_set1_epi16(-oe_del); \
	oe_ins_vec = _mm256_set1_epi16(-oe_ins); \
	e_del_vec = _mm256_set1_epi16(-e_del); \
	e_ins_vec = _mm256_set1_epi16(-e_ins); \
	__m256i match_sc_vec = _mm256_set1_epi16(a); \
	__m256i mis_sc_vec = _mm256_set1_epi16(-b); \
	__m256i amb_sc_vec = _mm256_set1_epi16(-1); \
	__m256i amb_vec = _mm256_set1_epi16(4); \
	for (i=0; i<SIMD_WIDTH; ++i) h_vec_mask[i] =  _mm256_loadu_si256((__m256i*) (&h_vec_int_mask[i]));

/*
 * e ref
 * f seq
 * m （，）
 * h 
 */
// load
#define SIMD_LOAD \
	__m256i m1 = _mm256_loadu_si256((__m256i*) (&mA1[j])); \
	__m256i e1 = _mm256_loadu_si256((__m256i*) (&eA1[j])); \
	__m256i m1j1 = _mm256_loadu_si256((__m256i*) (&mA1[j-1])); \
	__m256i f1j1 = _mm256_loadu_si256((__m256i*) (&fA1[j-1])); \
	__m256i h0j1 = _mm256_loadu_si256((__m256i*) (&hA0[j-1])); \
	__m256i qs_vec = _mm256_loadu_si256((__m256i*) (&seq[j-1])); \
	__m256i ts_vec = _mm256_loadu_si256((__m256i*) (&ref[tlen - i]));

// refseq，
#define SIMD_CMP_SEQ \
	__m256i match_mask_vec = _mm256_cmpeq_epi16(qs_vec, ts_vec); \
	__m256i mis_score_vec  = _mm256_andnot_si256(match_mask_vec, mis_sc_vec); \
	__m256i score_vec      = _mm256_and_si256(match_sc_vec, match_mask_vec); \
	score_vec = _mm256_or_si256(score_vec, mis_score_vec); \
	__m256i q_amb_mask_vec = _mm256_cmpeq_epi16(qs_vec, amb_vec); \
	__m256i t_amb_mask_vec = _mm256_cmpeq_epi16(ts_vec, amb_vec); \
	__m256i amb_mask_vec   = _mm256_or_si256(q_amb_mask_vec, t_amb_mask_vec); \
	score_vec = _mm256_andnot_si256(amb_mask_vec, score_vec); \
	__m256i amb_score_vec  = _mm256_and_si256(amb_mask_vec, amb_sc_vec); \
	score_vec = _mm256_or_si256(score_vec, amb_score_vec);

// h, e, f, m
#define SIMD_COMPUTE \
	__m256i en_vec0 = _mm256_add_epi16(m1, oe_del_vec); \
	__m256i en_vec1 = _mm256_add_epi16(e1, e_del_vec); \
	__m256i en_vec  = _mm256_max_epi16(en_vec0, en_vec1); \
	__m256i fn_vec0 = _mm256_add_epi16(m1j1, oe_ins_vec); \
	__m256i fn_vec1 = _mm256_add_epi16(f1j1, e_ins_vec); \
	__m256i fn_vec  = _mm256_max_epi16(fn_vec0, fn_vec1); \
	__m256i mn_vec0 = _mm256_add_epi16(h0j1, score_vec); \
	__m256i mn_mask = _mm256_cmpgt_epi16(h0j1, zero_vec); \
	__m256i mn_vec  = _mm256_and_si256(mn_vec0, mn_mask); \
	__m256i hn_vec0 = _mm256_max_epi16(en_vec, fn_vec); \
	__m256i hn_vec  = _mm256_max_epi16(hn_vec0, mn_vec); \
	en_vec = _mm256_max_epi16(en_vec, zero_vec); \
	fn_vec = _mm256_max_epi16(fn_vec, zero_vec); \
	mn_vec = _mm256_max_epi16(mn_vec, zero_vec); \
	hn_vec = _mm256_max_epi16(hn_vec, zero_vec); 

// 
#define SIMD_STORE \
	max_vec = _mm256_max_epi16(max_vec, hn_vec); \
	_mm256_storeu_si256((__m256i*)&eA2[j], en_vec); \
	_mm256_storeu_si256((__m256i*)&fA2[j], fn_vec); \
	_mm256_storeu_si256((__m256i*)&mA2[j], mn_vec); \
	_mm256_storeu_si256((__m256i*)&hA2[j], hn_vec);

// 
#define SIMD_REMOVE_EXTRA \
    en_vec = _mm256_and_si256(en_vec, h_vec_mask[end-j]); \
    fn_vec = _mm256_and_si256(fn_vec, h_vec_mask[end-j]); \
    mn_vec = _mm256_and_si256(mn_vec, h_vec_mask[end-j]); \
    hn_vec = _mm256_and_si256(hn_vec, h_vec_mask[end-j]);

// 
#define SIMD_FIND_MAX \
	max_vec = _mm256_max_epu16(max_vec, _mm256_alignr_epi8(max_vec, max_vec, 2)); \
	max_vec = _mm256_max_epu16(max_vec, _mm256_alignr_epi8(max_vec, max_vec, 4)); \
	max_vec = _mm256_max_epu16(max_vec, _mm256_alignr_epi8(max_vec, max_vec, 6)); \
	max_vec = _mm256_max_epu16(max_vec, _mm256_alignr_epi8(max_vec, max_vec, 8)); \
	max_vec = _mm256_max_epu16(max_vec, _mm256_permute2x128_si256(max_vec, max_vec, 0x01)); \
	int16_t *maxVal = (int16_t*)&max_vec; \
    m = MAX(m, maxVal[0]); /*BSW()*/ \
    if (maxVal[0] > 0 && m >= max) { \
		more_num++;	\
        for(j=beg, i=iend; j<=end; j+=SIMD_WIDTH, i-=SIMD_WIDTH) { \
			/*calc_num += 16;*/ \
            __m256i h2_vec = _mm256_loadu_si256((__m256i*) (&hA2[j])); \
            __m256i vcmp = _mm256_cmpeq_epi16(h2_vec, max_vec); \
            uint32_t mask = _mm256_movemask_epi8(vcmp); \
            if (mask > 0) { \
                int pos = SIMD_WIDTH - 1 - (( __builtin_clz(mask)) >> 1); \
                mj = j - 1 + pos; \
                mi = i - 1 - pos; \
				/*if (m >= max) fprintf(stderr, "%d %d %d %d %d %d %d\n", iend, beg, mi, mj, mask, pos, m);*/  \
            } \
        } \
    } else not_more_num++;

// ，
#define SWAP_DATA_POINTER \
	int16_t * tmp=hA0; \
	hA0 = hA1; hA1 = hA2; hA2 = tmp; \
	tmp = eA1; eA1 = eA2; eA2 = tmp; \
	tmp = fA1; fA1 = fA2; fA2 = tmp; \
	tmp = mA1; mA1 = mA2; mA2 = tmp;

static void write_query_target_sequence(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int h0, int fnum)
{
#ifdef DEBUG_FILE_OUTPUT
	// ，query.fa，target.fa，，info.txt，h0，qlen，tlen
	FILE *query_f = gfq[fnum], 
	*target_f = gft[fnum],
	*info_f = gfi[fnum];
	const char seq_map[5] = {'A', 'C', 'G', 'T', 'N'};
	int i;
	// query
	for (i = 0; i < qlen; ++i)
		fprintf(query_f, "%c", seq_map[query[i]]);
	fprintf(query_f, "\n");
	// target
	for (i = 0; i < tlen; ++i)
		fprintf(target_f, "%c", seq_map[target[i]]);
	fprintf(target_f, "\n");
	// 
	fprintf(info_f, "%-8d%-8d%-8d\n", qlen, tlen, h0);
#endif
}

int ksw_extend2_avx2(int qlen, // query length  query
		const uint8_t *query, // read
		int tlen, // target length reference
		const uint8_t *target, // reference
		int is_left, // 
		int m, //  (5)
		const int8_t *mat, // querytarget m*m
		int o_del, // deletion 
		int e_del, // deletion extension
		int o_ins, // insertion 
		int e_ins, // insertion extensionSIMD_BTYES
		int a, // match
		int b, // mismatch（）
		int w, // ，w =100   beg 
		int end_bonus,
		int zdrop,
		int h0, // seed（query）
		int *_qle, // query
		int *_tle, // reference
		int *_gtle, // querytarget
		int *_gscore, // query
		int *_max_off, // queryreference 
		buf_t *buf) // 
{
//	return ksw_extend2_origin(qlen, query, tlen, target, is_left, m, mat, o_del, e_del, o_ins, e_ins, w, end_bonus, zdrop, h0, _qle, _tle, _gtle, _gscore, _max_off);

#ifdef DEBUG_FILE_OUTPUT
	//fprintf(gf[0], "%d\n", qlen);
#ifdef GET_DIFFERENT_EXTENSION_LENGTH
	if (qlen <= 30) {
		write_query_target_sequence(qlen, query, tlen, target, h0, 0);
	} else if (qlen < 60) {
		write_query_target_sequence(qlen, query, tlen, target, h0, 1);
	} else if (qlen < 90) {
		write_query_target_sequence(qlen, query, tlen, target, h0, 2);
	} else {
		write_query_target_sequence(qlen, query, tlen, target, h0, 3);
	}
#endif
#endif

	if (qlen * a + h0 < 255) return ksw_extend2_avx2_u8(qlen, query, tlen, target, is_left, m, mat, o_del, e_del, o_ins, e_ins, a, b, w, end_bonus, zdrop, h0, _qle, _tle, _gtle, _gscore, _max_off, buf);
	
	int16_t *mA,*hA, *eA, *fA, *mA1, *mA2, *hA0, *hA1, *eA1, *fA1, *hA2, *eA2, *fA2; // hA0colH，H E F M
	int16_t *seq, *ref;
	uint8_t *mem;
	int16_t *qtmem, *vmem;
	int seq_size = qlen + SIMD_WIDTH, ref_size = tlen + SIMD_WIDTH;
	int i, ibeg, D, j, k, beg, end, max, max_i, max_j, max_ins, max_del, max_ie, gscore, max_off;
	int Dloop = tlen + qlen; // 
	int span, beg1, end1; // 
	int col_size = qlen + 2 + SIMD_WIDTH;
	int val_mem_size = (col_size * 9 * 2 + 31) >> 5 << 5; // 32
	int mem_size = (seq_size + ref_size) * 2 + val_mem_size;

	SIMD_INIT; // simd

	assert(h0 > 0);

	// allocate memory
	//mem = malloc(mem_size);

	if (buf->m < mem_size) {
		buf->m = mem_size;
		buf->addr = realloc(buf->addr, mem_size);
	}
	mem = buf->addr;

	qtmem = (int16_t *)&mem[0];
	seq=&qtmem[0]; ref=&qtmem[seq_size];
	if (is_left) {
		for (i=0; i<qlen; ++i) seq[i] = query[qlen - 1 - i];
		for (i=0; i<tlen; ++i) ref[i] = target[i];
	} else {
		for (i=0; i<qlen; ++i) seq[i] = query[i];
		for (i=0; i<tlen; ++i) ref[i] = target[tlen - 1 - i];
	}

	vmem = &ref[ref_size];
	for (i=0; i<(val_mem_size>>1); i+=SIMD_WIDTH) {
		_mm256_storeu_si256((__m256i*)&vmem[i], zero_vec);
	}
	hA = &vmem[0];
	mA = &vmem[col_size * 3];
	eA = &vmem[col_size * 5];
	fA = &vmem[col_size * 7];

	hA0 = &hA[0]; hA1 = &hA[col_size]; hA2 = &hA1[col_size];
	mA1 = &mA[0]; mA2 = &mA[col_size];
	eA1 = &eA[0]; eA2 = &eA[col_size];
	fA1 = &fA[0]; fA2 = &fA[col_size];

	// adjust $w if it is too large
	k = m * m;
	// get the max score
	for (i = 0, max = 0; i < k; ++i) max = max > mat[i]? max : mat[i];
	max_ins = (int)((double)(qlen * max + end_bonus - o_ins) / e_ins + 1.);
	max_ins = max_ins > 1? max_ins : 1;
	w = w < max_ins? w : max_ins;
	max_del = (int)((double)(qlen * max + end_bonus - o_del) / e_del + 1.);
	max_del = max_del > 1? max_del : 1;
	w = w < max_del? w : max_del; // TODO: is this necessary?
	if (tlen < qlen) w = MIN(tlen - 1, w);

	// DP loop
	max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1;;
	max_off = 0;
	beg = 1; end = qlen;
	// init h0
	hA0[0] = h0; // 

	if (qlen == 0 || tlen == 0) Dloop = 0; // 
    if (w >= qlen) { max_ie = 0; gscore = 0; }

	int m_last=0;
    int iend;

#ifdef ELIMINATE_DIFF_1
	int midx = 1, icheck = 0, checkspecial = 1;
	int m3 = 0, m2 = 0, m1 = 0;
	//int marr[10] = {0};
	//int marr[b]; memset(marr, 0, 4 * b);
#endif

	//int print_flag = 0; //(qlen == 64 && tlen == 123);
#ifdef DEBUG_SW_EXTEND
	int dii, djj;
	int16_t ins[tlen + 1][qlen + 2];
	int16_t del[tlen + 1][qlen + 2];
	int16_t score[tlen + 1][qlen + 2];
	for (dii = 0; dii <= tlen; ++dii)
	{
		for (djj = 0; djj <= qlen; ++djj)
		{
			ins[dii][djj] = del[dii][djj] = score[dii][djj] = NO_VAL;
		}
	}
	for (dii = 1; dii <= tlen; ++dii)
	{
		del[dii][0] = MAX(0, h0 - o_del - e_del * dii);
		score[dii][0] = del[dii][0];
	}
	for (djj = 1; djj <= qlen; ++djj)
	{
		ins[0][djj] = MAX(0, h0 - o_ins - e_ins * djj);
		score[0][djj] = ins[0][djj];
	}
	ins[0][0] = del[0][0] = score[0][0] = h0;
#endif

	for (D = 1; LIKELY(D < Dloop); ++D) {
		// ！ tlen ，， qlen
		if (D > tlen) {
			span = MIN(Dloop-D, w);
			beg1 = MAX(D-tlen+1, ((D-w) / 2) + 1);
		} else {
			span = MIN(D-1, w);
			beg1 = MAX(1, ((D-w) / 2) + 1);
		}
		end1 = MIN(qlen, beg1+span);

		if (beg < beg1) beg = beg1;
		if (end > end1) end = end1;
		if (beg > end) break; // ，，hA2，hA0，bug

        iend = D - (beg - 1); // ref，
		span = end - beg;
		ibeg = iend - span - 1; // 0ref

		// 
		int m = 0, mj = -1, mi = -1;
	    max_vec = zero_vec;
		//if (print_flag)
		//{
			//fprintf(stderr, "D: %d, iend: %d, jbeg: %d\n", D, iend, beg);
		//}
		// 
		//  f (insert)
		if (ibeg == 0) { hA1[end] = MAX(0, h0 - (o_ins + e_ins * end)); m = hA1[end];}
		// 
		if (beg == 1) { hA1[0] = MAX(0, h0 - (o_del + e_del * iend)); }
		else if (D & 1) {
			hA1[beg - 1] = 0;
			hA2[beg - 1] = 0;
		}

		for (j=beg, i=iend; j<=end+1-SIMD_WIDTH; j+=SIMD_WIDTH, i-=SIMD_WIDTH) {
			// 
			SIMD_LOAD;
			// seq，
			SIMD_CMP_SEQ;
			// 
			SIMD_COMPUTE;
            //calc_num += 16;
            // 
            SIMD_STORE;
		}
		// 
		if (j <= end) {
			// 
			SIMD_LOAD;
			// seq，
			SIMD_CMP_SEQ;
			// 
			SIMD_COMPUTE;
            //calc_num += 16;
            // 
            SIMD_REMOVE_EXTRA;
			// 
			SIMD_STORE;
		}

		SIMD_FIND_MAX;

#ifdef ELIMINATE_DIFF_1
// BSW()
#if 0
		if (hA1[0] < b && checkspecial) {
			int mi;
			if (hA1[0] == b - 1) {
				icheck = iend + 1;
			}
			for (mi = 0; mi < b - 1; ++mi) {
				if (midx - mi > 0)
					marr[mi] = MAX(marr[mi], hA2[midx - mi]);
			}
			midx += 1;
			if (ibeg > icheck)
			{
				int stopCalc = 0;
				for (mi = 0; mi < b - 1; ++mi)
				{
					stopCalc |= !marr[mi];
				}
				if (stopCalc)
					break;
				else
					checkspecial = 0;
			}
		}
#else
		if (hA1[0] < 4 && checkspecial) { // b == 4
			if (hA1[0] == 3) {
				icheck = iend + 1;
			} else if (midx == 2) {
				m2 = MAX(m2, hA2[midx - 1]);
			} else {
				m2 = MAX(m2, hA2[midx - 1]);
				m1 = MAX(m1, hA2[midx - 2]);
			}
			m3 = MAX(m3, hA2[midx]);
			midx += 1;
			if (ibeg > icheck)
			{
				if (!m1 || !m2 || !m3)
					break;
				else
					checkspecial = 0;
			}

			//if (print_flag) {
				//fprintf(stderr, "jbeg: %d, ibeg: %d, iend: %d, icheck: %d, score: %d %d %d, j: %d\n", beg, ibeg, iend, icheck, hA2[midx + 1], hA2[midx + 2], hA2[midx + 3], midx);
				//if (midx > 2) fprintf(stderr, "%d, %d, %d\n", hA2[midx-1], hA2[midx-2], hA2[midx-3]);
				//fprintf(stderr, "jbeg: %d, ibeg: %d, iend: %d, icheck: %d, hA1: %d, score: %d %d %d, j: %d\n", beg, ibeg, iend, icheck, hA1[0], m1, m2, m3, midx);
			//}
		}
#endif
#endif

#ifdef DEBUG_SW_EXTEND
		for (djj = beg; djj <= end; ++djj)
		{
			dii = D - djj + 1;
			ins[dii][djj] = fA2[djj];
			del[dii][djj] = eA2[djj];
			score[dii][djj] = hA2[djj];
		}
		//if (print_flag)
		//{
			//fprintf(stderr, "score: %d %d %d\n", hA2[beg], hA2[beg+1], hA2[beg+2]);
		//}
#endif

		// j
		j = end + 1;

		if (j == qlen + 1) {
			max_ie = gscore > hA2[qlen] ? max_ie : ibeg;
			gscore = gscore > hA2[qlen] ? gscore : hA2[qlen];
		}
		if (m == 0 && m_last==0) break; // ，
		//if (m == 0 && m_last < 2) break;
		if (m > max) {
			max = m, max_i = mi, max_j = mj;
			max_off = max_off > abs(mj - mi) ? max_off : abs(mj - mi);
		} else if (m == max && max_i >= mi && mj > max_j) {
			max_i = mi, max_j = mj;
			max_off = max_off > abs(mj - mi) ? max_off : abs(mj - mi);
		}
		else if (zdrop > 0 && mi > -1) {
			if (mi - max_i > mj - max_j) {
				if (max - m - ((mi - max_i) - (mj - max_j)) * e_del > zdrop) break;
			} else {
				if (max - m - ((mj - max_j) - (mi - max_i)) * e_ins > zdrop) break;
			}
		}

		// 
		for (j = beg; LIKELY(j <= end); ++j) { int has_val = hA1[j-1] | hA2[j]; if (has_val) break; }
        beg = j;
		for (j = end+1; LIKELY(j >= beg); --j) { int has_val = hA1[j-1] | hA2[j]; if (has_val) break; else hA0[j-1]=0; }
		end = j + 1 <= qlen? j + 1 : qlen;

        m_last = m;
		// swap m, h, e, f
		SWAP_DATA_POINTER;
	}

#ifdef DEBUG_FILE_OUTPUT
#ifdef DEBUG_SW_EXTEND
	fprintf(gf[0], "qlen: %d, tlen: %d, h0: %d, w: %d, mi: %d, mj: %d, mie: %d, max_off: %d, score: %d, max: %d\n", qlen, tlen, h0, w, max_i + 1, max_j + 1, max_ie + 1, max_off, gscore, max);
	fprintf(gf[1], "qlen: %d, tlen: %d, h0: %d, w: %d, mi: %d, mj: %d, mie: %d, max_off: %d, score: %d, max: %d\n", qlen, tlen, h0, w, max_i + 1, max_j + 1, max_ie + 1, max_off, gscore, max);
	fprintf(gf[2], "qlen: %d, tlen: %d, h0: %d, w: %d, mi: %d, mj: %d, mie: %d, max_off: %d, score: %d, max: %d\n", qlen, tlen, h0, w, max_i + 1, max_j + 1, max_ie + 1, max_off, gscore, max);

	fprintf(gf[0], "%-4d", -1);
	fprintf(gf[1], "%-4d", -1);
	fprintf(gf[2], "%-4d", -1);
	fprintf(gf[0], "%-4d", -1);
	fprintf(gf[1], "%-4d", -1);
	fprintf(gf[2], "%-4d", -1);
	for (djj = 0; djj < qlen; ++djj) {
		fprintf(gf[0], "%-4c", "ACGTN"[query[djj]]);
		fprintf(gf[1], "%-4c", "ACGTN"[query[djj]]);
		fprintf(gf[2], "%-4c", "ACGTN"[query[djj]]);
	}
	fprintf(gf[0], "\n");
	fprintf(gf[1], "\n");
	fprintf(gf[2], "\n");
	for (dii = 0; dii <= tlen; ++dii)
	{
		if (dii > 0) {
			fprintf(gf[0], "%-4c", "ACGTN"[target[dii - 1]]);
			fprintf(gf[1], "%-4c", "ACGTN"[target[dii - 1]]);
			fprintf(gf[2], "%-4c", "ACGTN"[target[dii - 1]]);
		} else {
			fprintf(gf[0], "%-4d", -1);
			fprintf(gf[1], "%-4d", -1);
			fprintf(gf[2], "%-4d", -1);
		}
		for (djj = 0; djj <= qlen; ++djj)
		{
			fprintf(gf[0], "%-4d", score[dii][djj]);
			fprintf(gf[1], "%-4d", ins[dii][djj]);
			fprintf(gf[2], "%-4d", del[dii][djj]);
		}
		fprintf(gf[0], "\n");
		fprintf(gf[1], "\n");
		fprintf(gf[2], "\n");
	}
#endif
#endif
	// if (max > h0)
    //     more_num++;
	// else
    //     not_more_num++;
    //	free(mem);
    if (_qle) *_qle = max_j + 1;
    if (_tle) *_tle = max_i + 1;
    if (_gtle) *_gtle = max_ie + 1;
    if (_gscore) *_gscore = gscore;
    if (_max_off) *_max_off = max_off;
	return max;
}

typedef struct {
	int32_t h, e;
} eh_t;

int ksw_extend2_origin(int qlen, // query length  query
                const uint8_t *query, // read
                int tlen, // target length reference
                const uint8_t *target, // reference
				int is_left, // 
                int m, //  (5)
                const int8_t *mat, // querytarget m*m
                int o_del, // deletion 
                int e_del, // deletion extension
                int o_ins, // insertion 
                int e_ins, // insertion extension
                int w, // ，w =100   beg 
                int end_bonus, 
                int zdrop, 
                int h0, // seed（query）
                int *_qle, // query
                int *_tle, // reference
                int *_gtle, // querytarget
                int *_gscore, // query
                int *_max_off) // queryreference 
{
	eh_t *eh; // score array
	int8_t *qp; // query profile
	int i, j, k, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, beg, end, max, max_i, max_j, max_ins, max_del, max_ie, gscore, max_off;
	uint8_t *qmem, *ref, *seq;
	assert(h0 > 0);
	// allocate memory
	qp = malloc(qlen * m);
	eh = calloc(qlen + 1, 8);
	qmem = malloc(qlen + tlen);
	seq=(uint8_t*)&qmem[0]; ref=(uint8_t*)&qmem[qlen];
	if (is_left) {
		for (i=0; i<qlen; ++i) seq[i] = query[qlen - 1 - i];
		for (i=0; i<tlen; ++i) ref[i] = target[tlen - 1 - i];
	} else {
		for (i=0; i<qlen; ++i) seq[i] = query[i];
		for (i=0; i<tlen; ++i) ref[i] = target[i];
	}
	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[seq[j]];
	}
	// fill the first row
	eh[0].h = h0; eh[1].h = h0 > oe_ins? h0 - oe_ins : 0;
	for (j = 2; j <= qlen && eh[j-1].h > e_ins; ++j)
		eh[j].h = eh[j-1].h - e_ins;
	// adjust $w if it is too large
	k = m * m;
	for (i = 0, max = 0; i < k; ++i) // get the max score
		max = max > mat[i]? max : mat[i];
	max_ins = (int)((double)(qlen * max + end_bonus - o_ins) / e_ins + 1.);
	max_ins = max_ins > 1? max_ins : 1;
	w = w < max_ins? w : max_ins;
	max_del = (int)((double)(qlen * max + end_bonus - o_del) / e_del + 1.);
	max_del = max_del > 1? max_del : 1;
	w = w < max_del? w : max_del; // TODO: is this necessary?
    //fprintf(stderr, "%d\n", w);
	// DP loop
	max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1;
	max_off = 0;
	beg = 0, end = qlen;

	//int print_flag = 0; //(qlen == 116 && tlen == 241);
	//fprintf(stderr, "%d %d %d\n", print_flag, qlen, tlen);
#ifdef DEBUG_SW_EXTEND
	int dii, djj;
	int16_t ins[tlen + 1][qlen + 2];
	int16_t del[tlen + 1][qlen + 2];
	int16_t score[tlen + 1][qlen + 2];
	for (dii = 0; dii <= tlen; ++dii)
	{
		for (djj = 0; djj <= qlen; ++djj)
		{
			ins[dii][djj] = del[dii][djj] = score[dii][djj] = NO_VAL;
		}
	}
	for (dii = 1; dii <= tlen; ++dii)
	{
		del[dii][0] = MAX(0, h0 - o_del - e_del * dii);
		score[dii][0] = del[dii][0];
	}
	for (djj = 1; djj <= qlen; ++djj)
	{
		ins[0][djj] = MAX(0, h0 - o_ins - e_ins * djj);
		score[0][djj] = ins[0][djj];
	}
	ins[0][0] = del[0][0] = score[0][0] = h0;
#endif

#ifdef DEBUG_FILE_OUTPUT
#ifdef COUNT_CALC_NUM
	int bsw_cal_num = 0;
	int real_cal_num = 0;
	for (i = 0; i < tlen; ++i)
	{
		int beg = MAX(0, i - w);
		int end = MIN(qlen, i + w + 1);
		if (beg >= end) break;
		bsw_cal_num += end - beg;
	}
	fprintf(gf[0], "start\n%d\n", bsw_cal_num);
#endif
#endif

#ifdef ELIMINATE_DIFF_3
	int prun_end = qlen; // for test diff_3
#endif

	for (i = 0; LIKELY(i < tlen); ++i) {
		int t, f = 0, h1, m = 0, mj = -1;
		int8_t *q = &qp[ref[i] * qlen];
		// apply the band and the constraint (if provided)
		if (beg < i - w) beg = i - w;
		if (end > i + w + 1) end = i + w + 1;
		if (end > qlen) end = qlen; // 
		// compute the first column
		if (beg == 0) {
			h1 = h0 - (o_del + e_del * (i + 1));
			if (h1 < 0) h1 = 0;
		} else h1 = 0;
		//m = h1; // VP-BSW()
		for (j = beg; LIKELY(j < end); ++j) {
            //calc_num++;

#ifdef DEBUG_FILE_OUTPUT
#ifdef COUNT_CALC_NUM
			real_cal_num++;
#endif
#endif

#ifdef DEBUG_SW_EXTEND
			ins[i+1][j+1] = f;
#endif
			// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
			// Similar to SSE2-SW, cells are computed in the following order:
			//   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			//   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
			//   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
			eh_t *p = &eh[j];
			int h, M = p->h, e = p->e; // get H(i-1,j-1) and E(i-1,j)
			p->h = h1;          // set H(i,j-1) for the next row
			M = M? M + q[j] : 0;// separating H and M to disallow a cigar like "100M3I3D20M",0，swnw
			h = M > e? M : e;   // e and f are guaranteed to be non-negative, so h>=0 even if M<0
			h = h > f? h : f;
#ifdef ELIMINATE_DIFF_3
			if (j >= prun_end && h==0) break; // for test diff_3
#endif
			h1 = h;             // save H(i,j) to h1 for the next column

#ifdef DEBUG_SW_EXTEND
			score[i+1][j+1] = h;
#endif
			mj = m > h? mj : j; // record the position where max score is achieved
			m = m > h? m : h;   // m is stored at eh[mj+1]
			t = M - oe_del;
			t = t > 0? t : 0;
			e -= e_del;
#ifdef DEBUG_SW_EXTEND
			del[i + 1][j + 1] = e;
#endif
			e = e > t? e : t;   // computed E(i+1,j)

#ifdef DEBUG_SW_EXTEND
//			del[i+1][j+1] = e;
#endif
			p->e = e;           // save E(i+1,j) for the next row
			t = M - oe_ins;
			t = t > 0? t : 0;
			f -= e_ins;
			f = f > t? f : t;   // computed F(i,j+1)
		}
		eh[end].h = h1; eh[end].e = 0;
		if (j == qlen) {
			max_ie = gscore > h1? max_ie : i;
			gscore = gscore > h1? gscore : h1;
		}
		if (m == 0) break;
		if (m > max) {
			max = m, max_i = i, max_j = mj;
			max_off = max_off > abs(mj - i)? max_off : abs(mj - i);
			//fprintf(stderr, "%d %d %d %d\n", i, mj, max_off, m);
		} else if (zdrop > 0) {
			if (i - max_i > mj - max_j) {
				if (max - m - ((i - max_i) - (mj - max_j)) * e_del > zdrop) break;
			} else {
				if (max - m - ((mj - max_j) - (i - max_i)) * e_ins > zdrop) break;
			}
		}
		// update beg and end for the next round
		for (j = beg; LIKELY(j < end) && eh[j].h == 0 && eh[j].e == 0; ++j); // f（insert score）
		beg = j;
		for (j = end; LIKELY(j >= beg) && eh[j].h == 0 && eh[j].e == 0; --j);
#ifdef ELIMINATE_DIFF_3
		prun_end = j + 2 < qlen ? j + 2 : qlen;  end = qlen; // for test diff_3
#else
		end = j + 2 < qlen? j + 2 : qlen;
#endif
		// beg = 0; end = qlen; // uncomment this line for debugging
		// if (print_flag) {
		//	fprintf(stderr, "beg: %d; end: %d\n", beg, end);
		// }
	}
#ifdef DEBUG_FILE_OUTPUT
#ifdef DEBUG_SW_EXTEND
	fprintf(gf[0], "qlen: %d, tlen: %d, h0: %d, w: %d, mi: %d, mj: %d, mie: %d, max_off: %d, score: %d, max: %d\n", qlen, tlen, h0, w, max_i + 1, max_j + 1, max_ie + 1, max_off, gscore, max);
	fprintf(gf[1], "qlen: %d, tlen: %d, h0: %d, w: %d, mi: %d, mj: %d, mie: %d, max_off: %d, score: %d, max: %d\n", qlen, tlen, h0, w, max_i + 1, max_j + 1, max_ie + 1, max_off, gscore, max);
	fprintf(gf[2], "qlen: %d, tlen: %d, h0: %d, w: %d, mi: %d, mj: %d, mie: %d, max_off: %d, score: %d, max: %d\n", qlen, tlen, h0, w, max_i + 1, max_j + 1, max_ie + 1, max_off, gscore, max);

	fprintf(gf[0], "%-4d", -1);
	fprintf(gf[1], "%-4d", -1);
	fprintf(gf[2], "%-4d", -1);
	fprintf(gf[0], "%-4d", -1);
	fprintf(gf[1], "%-4d", -1);
	fprintf(gf[2], "%-4d", -1);
	for (djj = 0; djj < qlen; ++djj)
	{
		fprintf(gf[0], "%-4c", "ACGTN"[query[djj]]);
		fprintf(gf[1], "%-4c", "ACGTN"[query[djj]]);
		fprintf(gf[2], "%-4c", "ACGTN"[query[djj]]);
	}
	fprintf(gf[0], "\n");
	fprintf(gf[1], "\n");
	fprintf(gf[2], "\n");
	for (dii = 0; dii <= tlen; ++dii)
	{
		if (dii > 0)
		{
			fprintf(gf[0], "%-4c", "ACGTN"[target[dii - 1]]);
			fprintf(gf[1], "%-4c", "ACGTN"[target[dii - 1]]);
			fprintf(gf[2], "%-4c", "ACGTN"[target[dii - 1]]);
		}
		else
		{
			fprintf(gf[0], "%-4d", -1);
			fprintf(gf[1], "%-4d", -1);
			fprintf(gf[2], "%-4d", -1);
		}
		for (djj = 0; djj <= qlen; ++djj)
		{
			fprintf(gf[0], "%-4d", score[dii][djj]);
			fprintf(gf[1], "%-4d", ins[dii][djj]);
			fprintf(gf[2], "%-4d", del[dii][djj]);
		}
		fprintf(gf[0], "\n");
		fprintf(gf[1], "\n");
		fprintf(gf[2], "\n");
	}
#endif
#endif

#ifdef DEBUG_FILE_OUTPUT
#ifdef COUNT_CALC_NUM
	fprintf(gf[0], "%d\nend\n", real_cal_num);
#endif
#endif

	free(eh); free(qp); free(qmem);
	if (_qle) *_qle = max_j + 1;
	if (_tle) *_tle = max_i + 1;
	if (_gtle) *_gtle = max_ie + 1;
	if (_gscore) *_gscore = gscore;
	if (_max_off) *_max_off = max_off;
	return max;
}
