#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <stdio.h>
#include <immintrin.h>
#include <emmintrin.h>
#include "utils.h"

#define ELIMINATE_DIFF_1

#define SIMD_WIDTH 32

static const uint8_t h_vec_int_mask[SIMD_WIDTH][SIMD_WIDTH] = {
	{0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0},
	{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff}
};

//static const uint8_t reverse_mask[SIMD_WIDTH] = {7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8};
#define permute_mask _MM_SHUFFLE(0, 1, 2, 3)
//const int permute_mask = _MM_SHUFFLE(0, 1, 2, 3);
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
	oe_del_vec = _mm256_set1_epi8(oe_del); \
	oe_ins_vec = _mm256_set1_epi8(oe_ins); \
	e_del_vec = _mm256_set1_epi8(e_del); \
	e_ins_vec = _mm256_set1_epi8(e_ins); \
	__m256i match_sc_vec = _mm256_set1_epi8(a); \
	__m256i mis_sc_vec = _mm256_set1_epi8(b); \
	__m256i amb_sc_vec = _mm256_set1_epi8(1); \
	__m256i amb_vec = _mm256_set1_epi8(4); \
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
	__m256i match_mask_vec = _mm256_cmpeq_epi8(qs_vec, ts_vec); \
	__m256i mis_score_vec  = _mm256_andnot_si256(match_mask_vec, mis_sc_vec); \
	__m256i match_score_vec   = _mm256_and_si256(match_sc_vec, match_mask_vec); \
	__m256i q_amb_mask_vec = _mm256_cmpeq_epi8(qs_vec, amb_vec); \
	__m256i t_amb_mask_vec = _mm256_cmpeq_epi8(ts_vec, amb_vec); \
	__m256i amb_mask_vec   = _mm256_or_si256(q_amb_mask_vec, t_amb_mask_vec); \
	__m256i amb_score_vec  = _mm256_and_si256(amb_mask_vec, amb_sc_vec); \
    mis_score_vec = _mm256_andnot_si256(amb_mask_vec, mis_score_vec); \
    mis_score_vec = _mm256_or_si256(amb_score_vec, mis_score_vec); \
    match_score_vec = _mm256_andnot_si256(amb_mask_vec, match_score_vec);

// h, e, f, m
#define SIMD_COMPUTE \
	__m256i en_vec0 = _mm256_max_epu8(m1, oe_del_vec); \
    en_vec0 = _mm256_subs_epu8(en_vec0, oe_del_vec); \
	__m256i en_vec1 = _mm256_max_epu8(e1, e_del_vec); \
    en_vec1 = _mm256_subs_epu8(en_vec1, e_del_vec); \
	__m256i en_vec  = _mm256_max_epu8(en_vec0, en_vec1); \
	__m256i fn_vec0 = _mm256_max_epu8(m1j1, oe_ins_vec); \
    fn_vec0 = _mm256_subs_epu8(fn_vec0, oe_ins_vec); \
	__m256i fn_vec1 = _mm256_max_epu8(f1j1, e_ins_vec); \
    fn_vec1 = _mm256_subs_epu8(fn_vec1, e_ins_vec); \
	__m256i fn_vec  = _mm256_max_epu8(fn_vec0, fn_vec1); \
	__m256i mn_vec0 = _mm256_adds_epu8(h0j1, match_score_vec); \
    mn_vec0 = _mm256_max_epu8(mn_vec0, mis_score_vec); \
    mn_vec0 = _mm256_subs_epu8(mn_vec0, mis_score_vec); \
	__m256i mn_mask = _mm256_cmpeq_epi8(h0j1, zero_vec); \
	__m256i mn_vec  = _mm256_andnot_si256(mn_mask, mn_vec0); \
	__m256i hn_vec0 = _mm256_max_epu8(en_vec, fn_vec); \
	__m256i hn_vec  = _mm256_max_epu8(hn_vec0, mn_vec); 

// 
#define SIMD_STORE \
	max_vec = _mm256_max_epu8(max_vec, hn_vec); \
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
	uint8_t *maxVal = (uint8_t*)&max_vec; \
	max_vec = _mm256_max_epu8(max_vec, _mm256_alignr_epi8(max_vec, max_vec, 1)); \
	max_vec = _mm256_max_epu8(max_vec, _mm256_alignr_epi8(max_vec, max_vec, 2)); \
	max_vec = _mm256_max_epu8(max_vec, _mm256_alignr_epi8(max_vec, max_vec, 3)); \
	max_vec = _mm256_max_epu8(max_vec, _mm256_alignr_epi8(max_vec, max_vec, 4)); \
	max_vec = _mm256_max_epu8(max_vec, _mm256_alignr_epi8(max_vec, max_vec, 5)); \
	max_vec = _mm256_max_epu8(max_vec, _mm256_alignr_epi8(max_vec, max_vec, 6)); \
	max_vec = _mm256_max_epu8(max_vec, _mm256_alignr_epi8(max_vec, max_vec, 7)); \
	max_vec = _mm256_max_epu8(max_vec, _mm256_alignr_epi8(max_vec, max_vec, 8)); \
	max_vec = _mm256_max_epu8(max_vec, _mm256_permute2x128_si256(max_vec, max_vec, 0x01)); \
    m = MAX(m, maxVal[0]); \
    if (maxVal[0] > 0 && m >= max) { \
        for(j=beg, i=iend; j<=end; j+=SIMD_WIDTH, i-=SIMD_WIDTH) { \
            __m256i h2_vec = _mm256_loadu_si256((__m256i*) (&hA2[j])); \
            __m256i vcmp = _mm256_cmpeq_epi8(h2_vec, max_vec); \
            uint32_t mask = _mm256_movemask_epi8(vcmp); \
            if (mask > 0) { \
                int pos = SIMD_WIDTH - 1 - __builtin_clz(mask); \
                mj = j - 1 + pos; \
                mi = i - 1 - pos; \
            } \
        } \
    }

// ，
#define SWAP_DATA_POINTER \
	uint8_t * tmp=hA0; \
	hA0 = hA1; hA1 = hA2; hA2 = tmp; \
	tmp = eA1; eA1 = eA2; eA2 = tmp; \
	tmp = fA1; fA1 = fA2; fA2 = tmp; \
	tmp = mA1; mA1 = mA2; mA2 = tmp;


int ksw_extend2_avx2_u8(int qlen, // query length  query
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
	uint8_t *mA,*hA, *eA, *fA, *mA1, *mA2, *hA0, *hA1, *eA1, *fA1, *hA2, *eA2, *fA2; // hA0colH，H E F M
	uint8_t *seq, *ref;
	uint8_t *mem, *qtmem, *vmem;
	int seq_size = qlen + SIMD_WIDTH, ref_size = tlen + SIMD_WIDTH;
	int i, ibeg, D, j, k, beg, end, max, max_i, max_j, max_ins, max_del, max_ie, gscore, max_off;
	int Dloop = tlen + qlen; // 
	int span, beg1, end1; // 
	int col_size = qlen + 2 + SIMD_WIDTH;
	int val_mem_size = (col_size * 9 + 31) >> 5 << 5; // 32
	int mem_size = seq_size + ref_size + val_mem_size;

	SIMD_INIT; // simd

	assert(h0 > 0);

	// allocate memory
	//mem = malloc(mem_size);
	if (buf->m < mem_size) {
		buf->m = mem_size;
		buf->addr = realloc(buf->addr, mem_size);
	}
	mem = buf->addr;

	qtmem = &mem[0];
	seq=(uint8_t*)&qtmem[0]; ref=(uint8_t*)&qtmem[seq_size];
	if (is_left) {
		for (i=0; i<qlen; ++i) seq[i] = query[qlen - 1 - i];
		for (i=0; i<tlen; ++i) ref[i] = target[i];
	} else {
		for (i=0; i<qlen; ++i) seq[i] = query[i];
		for (i=0; i<tlen; ++i) ref[i] = target[tlen - 1 - i];
	}

	vmem = &ref[ref_size];
	for (i=0; i<val_mem_size; i+=SIMD_WIDTH) {
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
	// int marr[10] = {0};
	// int marr[b]; memset(marr, 0, 4 * b);
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

		// 
		//  f (insert)
		if (ibeg == 0) { hA1[end] = MAX(0, h0 - (o_ins + e_ins * end)); m = hA1[end]; }
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
            //calc_num += 32;
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
            //calc_num += 32;
            // 
            SIMD_REMOVE_EXTRA;
			// 
			SIMD_STORE;
		}

		SIMD_FIND_MAX;

#ifdef ELIMINATE_DIFF_1
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
		if (hA1[0] < 4 && checkspecial)
		{ // b == 4
			if (hA1[0] == 3)
			{
				icheck = iend + 1;
			}
			else if (midx == 2)
			{
				m2 = MAX(m2, hA2[midx - 1]);
			}
			else
			{
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
		}
#endif
#endif

		// j
		j = end + 1;

		if (j == qlen + 1) {
			max_ie = gscore > hA2[qlen] ? max_ie : ibeg;
			gscore = gscore > hA2[qlen] ? gscore : hA2[qlen];
		}
		if (m == 0 && m_last==0) break; // ，
		if (m > max) {
			max = m, max_i = mi, max_j = mj;
			max_off = max_off > abs(mj - mi)? max_off : abs(mj - mi);
		}
		else if (m == max && max_i >= mi && mj > max_j) {
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

	//free(mem);
    if (_qle) *_qle = max_j + 1;
    if (_tle) *_tle = max_i + 1;
    if (_gtle) *_gtle = max_ie + 1;
    if (_gscore) *_gscore = gscore;
    if (_max_off) *_max_off = max_off;
	return max;
}
