#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <immintrin.h>
#include <emmintrin.h>
#include "utils.h"
#include "ksw.h"

// static const kswr_t g_defr = { 0, -1, -1, -1, -1, -1, -1 };

typedef struct _kswq_t {
	int qlen, slen;
	uint8_t shift, mdiff, max, size;
	__m256i *qp, *H0, *H1, *E, *Hmax;
} kswq_t;

/**
 * Initialize the query data structure
 *
 * @param size   Number of bytes used to store a score; valid valures are 1 or 2
 * @param qlen   Length of the query sequence
 * @param query  Query sequence
 * @param m      Size of the alphabet (ACGTN)
 * @param mat    Scoring matrix in a one-dimension array
 *
 * @return       Query data structure
 */
kswq_t *ksw_qinit_avx2(int size, int qlen, const uint8_t *query, int m, const int8_t *mat)
{
    kswq_t *q;
    int slen, a, tmp, p;
    size = size > 1 ? 2 : 1;
    p = 16 * (3 - size); // # values per __m256i
    slen = (qlen + p - 1) / p; // segmented length
    q = (kswq_t *)malloc(sizeof(kswq_t) + 512 + 32 * slen * (m + 4)); // a single block of memory
    q->qp = (__m256i *)(((size_t)q + sizeof(kswq_t) + 31) >> 5 << 5); // align memory
	q->H0 = q->qp + slen * m;
	q->H1 = q->H0 + slen;
	q->E  = q->H1 + slen;
	q->Hmax = q->E + slen;
	q->slen = slen; q->qlen = qlen; q->size = size;
    // compute shift
    tmp = m * m;
    	for (a = 0, q->shift = 127, q->mdiff = 0; a < tmp; ++a) { // find the minimum and maximum score
		if (mat[a] < (int8_t)q->shift) q->shift = mat[a];
		if (mat[a] > (int8_t)q->mdiff) q->mdiff = mat[a];
	}
    q->max = q->mdiff;
    q->shift = 256 - q->shift; // NB: q->shift is uint8_t
    q->mdiff += q->shift;      // this is the difference between the min and max scores
    // An example: p=8, qlen=19, slen=3 and segmentation:
	//  {{0,3,6,9,12,15,18,-1},{1,4,7,10,13,16,-1,-1},{2,5,8,11,14,17,-1,-1}}
	if (size == 1) {
		int8_t *t = (int8_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? 0 : ma[query[k]]) + q->shift;
		}
	} else {
		int16_t *t = (int16_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? 0 : ma[query[k]]);
		}
	}
    return q;
}