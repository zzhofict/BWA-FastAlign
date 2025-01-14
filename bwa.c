/* The MIT License

   Copyright (c) 2018-     Dana-Farber Cancer Institute
                 2009-2018 Broad Institute, Inc.
                 2008-2009 Genome Research Ltd. (GRL)

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <assert.h>
#include <unistd.h>
#include <pthread.h>
#include <limits.h>
#include "bntseq.h"
#include "bwa.h"
#include "ksw.h"
#include "utils.h"
#include "kstring.h"
#include "kvec.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

int bwa_verbose = 3;
int bwa_dbg = 0;
char bwa_rg_id[256];
char *bwa_pg;

/************************
 * Batch FASTA/Q reader *
 ************************/

#include "kseq.h"
KSEQ_DECLARE(gzFile)

static inline void trim_readno(kstring_t *s)
{
	if (s->l > 2 && s->s[s->l-2] == '/' && isdigit(s->s[s->l-1]))
		s->l -= 2, s->s[s->l] = 0;
}

static inline void dupkstring(const kstring_t *str, int dupempty, char **dstp, int *sm)
{
	if (!dupempty && str->l == 0) {
		if (*dstp) free(*dstp);
		*dstp = 0; *sm = 0;
	} else if (*dstp == 0 || *sm < str->l) {
		*sm = str->l;
		*dstp = realloc(*dstp, str->l + 1);
	}

	char *s = *dstp;
	if (!s) return;

	memcpy(s, str->s, str->l);
	s[str->l] = '\0';
}

static inline void kseq2bseq1(const kseq_t *ks, bseq1_t *s, int copy_comment)
{ // TODO: it would be better to allocate one chunk of memory, but probably it does not matter in practice
	dupkstring(&ks->name, 1, &s->name, &s->m_name);
	if (copy_comment) dupkstring(&ks->comment, 0, &s->comment, &s->m_comment);
	dupkstring(&ks->seq, 1, &s->seq, &s->m_seq);
	dupkstring(&ks->qual, 0, &s->qual, &s->m_qual);
	s->l_seq = ks->seq.l;
}

typedef struct {
    kseq_t *ks;
    bseq1_t *seq;
    int start_pos;
    int n_bound;
	int copy_comment;
    int ret_n;
    int ret_size;
	int ret_status;
	int chunk_size;
} read_data_t;

static void *thread_bseq_read(void *data) {
    read_data_t *d = (read_data_t*) data;
    kseq_t *ks = d->ks;
    bseq1_t *seqs = d->seq;
	int copy_comment = d->copy_comment;
	int chunk_size = d->chunk_size;
	int cur_n = 0, cur_pos = d->start_pos, size = 0;
	int ret_status = 1;

    //pthread_t thread_id = pthread_self();
    //fprintf(stderr, "Thread ID: %lu\n", thread_id);
    PROF_START(parse);
    while (cur_n < d->n_bound && (ret_status = kseq_read(ks)) >= 0) {

        trim_readno(&ks->name);
        kseq2bseq1(ks, seqs + cur_pos, copy_comment);
		seqs[cur_pos].id = cur_pos;
		size += seqs[cur_pos].l_seq;
        cur_pos += 2; cur_n += 1;

        if (size >= chunk_size) break;
	}
    PROF_END(gprof[G_parse_seq], parse);
    d->ret_n = cur_n; d->ret_size = size; d->ret_status = ret_status;
    return 0;
}

#define READ_ONE_SEQ(ksin)                    \
	trim_readno(&(ksin)->name);               \
	kseq2bseq1(ksin, &seqs[n], copy_comment); \
	seqs[n].id = n;                           \
	size += seqs[n++].l_seq;

// multi thread reading input seqs
void bseq_read_pe_mt(int chunk_size, int *n_, void *ks1_, void *ks2_, int copy_comment, int64_t *size_, int *m_, bseq1_t **seqs_ptr)
{
	kseq_t *ks = (kseq_t *)ks1_, *ks2 = (kseq_t *)ks2_;
	int size = 0, m = *m_, n = 0;
	bseq1_t *seqs = *seqs_ptr;
	read_data_t d[2];
	pthread_t tid[2];
	const int chunk_size_narrow = 4 * 1024 * 1024;
	const int init_n_reads = 20;
	if (m == 0) { // ，
		seqs = calloc(init_n_reads, sizeof(bseq1_t)); // 20reads，readschunk sizereads
#if 1
		int ks1_ret = 0, ks2_ret = 0;
		int i = init_n_reads >> 1;
		while (i-- > 0) {
			ks1_ret = kseq_read(ks);
			if (ks1_ret < 0) break;
			ks2_ret = kseq_read(ks2);
			if (ks2_ret < 0) {
				fprintf(stderr, "[W::%s] the 2nd file has fewer sequences.\n", __func__);
				break;
			}
			READ_ONE_SEQ(ks);
			READ_ONE_SEQ(ks2);
		}
		if (ks1_ret < 0 || ks2_ret < 0) {
			if (size == 0 && kseq_read(ks2) >= 0) { // test if the 2nd file is finished
				fprintf(stderr, "[W::%s] the 1st file has fewer sequences.\n", __func__);
			}
			*n_ = n; *seqs_ptr = seqs; *size_ = size; *m_ = n;
			return;
		}
		m = (chunk_size + size / init_n_reads - 1) / (size / init_n_reads);
#else
		m = 50000;
#endif
		seqs = realloc(seqs, m * sizeof(bseq1_t));
		memset(seqs + n, 0, sizeof(bseq1_t) * (m - n));
	}

	d[0].copy_comment = copy_comment; d[1].copy_comment = copy_comment;
    d[0].ks = ks ; d[0].seq = &seqs[0]; d[0].n_bound = (m >> 1) - (n>>1); d[0].start_pos = n;
    d[1].ks = ks2; d[1].seq = &seqs[0]; d[1].n_bound = (m >> 1) - (n>>1); d[1].start_pos = n+1;
	d[0].chunk_size = d[1].chunk_size = (chunk_size - chunk_size_narrow - size) >> 1;

	pthread_create(&tid[0], 0, thread_bseq_read, &d[0]);
    pthread_create(&tid[1], 0, thread_bseq_read, &d[1]);
    pthread_join(tid[0], 0); pthread_join(tid[1], 0);

	size += d[0].ret_size + d[1].ret_size;

	// reads
	if (d[0].ret_n < d[1].ret_n)
	{
		int num_to_read = d[1].ret_n - d[0].ret_n;
		int offset = n + d[0].ret_n * 2;
		while (num_to_read-- > 0 && kseq_read(ks) >= 0) {
			trim_readno(&ks->name);
			kseq2bseq1(ks, &seqs[offset], copy_comment);
			seqs[offset].id = offset;
			size += seqs[offset].l_seq;
			offset += 2;
		}
		d[0].ret_n = d[1].ret_n;
	}
	else if (d[1].ret_n < d[0].ret_n)
	{
		int num_to_read = d[0].ret_n - d[1].ret_n;
		int offset = n + 1 + d[1].ret_n * 2;
		while (num_to_read-- > 0 && kseq_read(ks2) >= 0) {
			trim_readno(&ks2->name);
			kseq2bseq1(ks2, &seqs[offset], copy_comment);
			seqs[offset].id = offset;
			size += seqs[offset].l_seq;
			offset += 2;
		}
		d[1].ret_n = d[0].ret_n;
	}

	n += d[0].ret_n + d[1].ret_n;

	if (size == 0 && kseq_read(ks2) >= 0) { // test if the 2nd file is finished
		fprintf(stderr, "[W::%s] the 1st file has fewer sequences.\n", __func__);
	} else if (size < chunk_size && d[0].ret_status > 0 && d[1].ret_status > 0) {
		while (kseq_read(ks) >= 0) {
			if (kseq_read(ks2) < 0) { // the 2nd file has fewer reads
				fprintf(stderr, "[W::%s] the 2nd file has fewer sequences.\n", __func__);
				break;
			}
			if (n >= m) {
				m = m? m<<1 : 256;
				seqs = realloc(seqs, m * sizeof(bseq1_t));
				memset(seqs + n, 0, (m-n) * sizeof(bseq1_t));
			}
			READ_ONE_SEQ(ks);
			READ_ONE_SEQ(ks2);
			if (size >= chunk_size && (n&1) == 0) break;
		}
		if (size == 0) { // test if the 2nd file is finished
			if (kseq_read(ks2) >= 0) fprintf(stderr, "[W::%s] the 1st file has fewer sequences.\n", __func__);
		}
	}
	*n_ = n; *size_ = size;
	if (m > *m_) *m_ = m;
	*seqs_ptr = seqs;
}

void bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_, int copy_comment, int64_t *size_, int *m_, bseq1_t **seqs_ptr)
{
#ifdef USE_MT_READ
	if (ks2_) return bseq_read_pe_mt(chunk_size, n_, ks1_, ks2_, copy_comment, size_, m_, seqs_ptr);
#endif
	kseq_t *ks = (kseq_t*)ks1_, *ks2 = (kseq_t*)ks2_;
	int size = 0, m, n;
	bseq1_t *seqs = *seqs_ptr;
	n = 0; m = *m_;
    while (kseq_read(ks) >= 0) {
		if (ks2 && kseq_read(ks2) < 0) { // the 2nd file has fewer reads
			fprintf(stderr, "[W::%s] the 2nd file has fewer sequences.\n", __func__);
			break;
		}
        PROF_START(parse);
        if (n >= m) {
			m = m? m<<1 : 256;
			seqs = realloc(seqs, m * sizeof(bseq1_t));
			memset(seqs + n, 0, (m-n) * sizeof(bseq1_t));
		}
		READ_ONE_SEQ(ks);
		if (ks2) {
			READ_ONE_SEQ(ks2);
		}
        PROF_END(gprof[G_parse_seq], parse);
        if (size >= chunk_size && (n&1) == 0) break;
	}

    if (size == 0) { // test if the 2nd file is finished
		if (ks2 && kseq_read(ks2) >= 0)
			fprintf(stderr, "[W::%s] the 1st file has fewer sequences.\n", __func__);
	}
	*n_ = n; *size_ = size;
	if (m > *m_) *m_ = m;
	*seqs_ptr = seqs;
}

void bseq_classify(int n, bseq1_t *seqs, int m[2], bseq1_t *sep[2])
{
	int i, has_last;
	kvec_t(bseq1_t) a[2] = {{0,0,0}, {0,0,0}};
	for (i = 1, has_last = 1; i < n; ++i) {
		if (has_last) {
			if (strcmp(seqs[i].name, seqs[i-1].name) == 0) {
				kv_push(bseq1_t, a[1], seqs[i-1]);
				kv_push(bseq1_t, a[1], seqs[i]);
				has_last = 0;
			} else kv_push(bseq1_t, a[0], seqs[i-1]);
		} else has_last = 1;
	}
	if (has_last) kv_push(bseq1_t, a[0], seqs[i-1]);
	sep[0] = a[0].a, m[0] = a[0].n;
	sep[1] = a[1].a, m[1] = a[1].n;
}

/*****************
 * CIGAR related *
 *****************/

void bwa_fill_scmat(int a, int b, int8_t mat[25])
{
	int i, j, k;
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? a : -b;
		mat[k++] = -1; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = -1;
}

// Generate CIGAR when the alignment end points are known
uint32_t *bwa_gen_cigar2(const int8_t mat[25], int o_del, int e_del, int o_ins, int e_ins, int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM)
{
	uint32_t *cigar = 0;
	uint8_t tmp, *rseq;
	int i;
	int64_t rlen;
	kstring_t str;
	const char *int2base;

	if (n_cigar) *n_cigar = 0;
	if (NM) *NM = -1;
	if (l_query <= 0 || rb >= re || (rb < l_pac && re > l_pac)) return 0; // reject if negative length or bridging the forward and reverse strand
	rseq = bns_get_seq(l_pac, pac, rb, re, &rlen);
	if (re - rb != rlen) goto ret_gen_cigar; // possible if out of range
	if (rb >= l_pac) { // then reverse both query and rseq; this is to ensure indels to be placed at the leftmost position
		for (i = 0; i < l_query>>1; ++i)
			tmp = query[i], query[i] = query[l_query - 1 - i], query[l_query - 1 - i] = tmp;
		for (i = 0; i < rlen>>1; ++i)
			tmp = rseq[i], rseq[i] = rseq[rlen - 1 - i], rseq[rlen - 1 - i] = tmp;
	}
	if (l_query == re - rb && w_ == 0) { // no gap; no need to do DP
		// UPDATE: we come to this block now... FIXME: due to an issue in mem_reg2aln(), we never come to this block. This does not affect accuracy, but it hurts performance.
		if (n_cigar) {
			cigar = malloc(4);
			cigar[0] = l_query<<4 | 0;
			*n_cigar = 1;
		}
		for (i = 0, *score = 0; i < l_query; ++i)
			*score += mat[rseq[i]*5 + query[i]];
	} else {
		int w, max_gap, max_ins, max_del, min_w;
		// set the band-width
		max_ins = (int)((double)(((l_query+1)>>1) * mat[0] - o_ins) / e_ins + 1.);
		max_del = (int)((double)(((l_query+1)>>1) * mat[0] - o_del) / e_del + 1.);
		max_gap = max_ins > max_del? max_ins : max_del;
		max_gap = max_gap > 1? max_gap : 1;
		w = (max_gap + abs((int)rlen - l_query) + 1) >> 1;
		w = w < w_? w : w_;
		min_w = abs((int)rlen - l_query) + 3;
		w = w > min_w? w : min_w;
		// NW alignment
		if (bwa_verbose >= 4) {
			printf("* Global bandwidth: %d\n", w);
			printf("* Global ref:   "); for (i = 0; i < rlen; ++i) putchar("ACGTN"[(int)rseq[i]]); putchar('\n');
			printf("* Global query: "); for (i = 0; i < l_query; ++i) putchar("ACGTN"[(int)query[i]]); putchar('\n');
		}
		*score = ksw_global2(l_query, query, rlen, rseq, 5, mat, o_del, e_del, o_ins, e_ins, w, n_cigar, &cigar);
	}
	if (NM && n_cigar) {// compute NM and MD
		int k, x, y, u, n_mm = 0, n_gap = 0;
		str.l = str.m = *n_cigar * 4; str.s = (char*)cigar; // append MD to CIGAR
		int2base = rb < l_pac? "ACGTN" : "TGCAN";
		for (k = 0, x = y = u = 0; k < *n_cigar; ++k) {
			int op, len;
			cigar = (uint32_t*)str.s;
			op  = cigar[k]&0xf, len = cigar[k]>>4;
			if (op == 0) { // match
				for (i = 0; i < len; ++i) {
					if (query[x + i] != rseq[y + i]) {
						kputw(u, &str);
						kputc(int2base[rseq[y+i]], &str);
						++n_mm; u = 0;
					} else ++u;
				}
				x += len; y += len;
			} else if (op == 2) { // deletion
				if (k > 0 && k < *n_cigar - 1) { // don't do the following if D is the first or the last CIGAR
					kputw(u, &str); kputc('^', &str);
					for (i = 0; i < len; ++i)
						kputc(int2base[rseq[y+i]], &str);
					u = 0; n_gap += len;
				}
				y += len;
			} else if (op == 1) x += len, n_gap += len; // insertion
		}
		kputw(u, &str); kputc(0, &str);
		*NM = n_mm + n_gap;
		cigar = (uint32_t*)str.s;
	}
	if (rb >= l_pac) // reverse back query
		for (i = 0; i < l_query>>1; ++i)
			tmp = query[i], query[i] = query[l_query - 1 - i], query[l_query - 1 - i] = tmp;

ret_gen_cigar:
	free(rseq);
	return cigar;
}

uint32_t *bwa_gen_cigar(const int8_t mat[25], int q, int r, int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM)
{
	return bwa_gen_cigar2(mat, q, r, q, r, w_, l_pac, pac, l_query, query, rb, re, score, n_cigar, NM);
}

/*********************
 * Full index reader *
 *********************/

char *bwa_idx_infer_prefix(const char *hint)
{
	char *prefix;
	int l_hint;
	FILE *fp;
	l_hint = strlen(hint);
	prefix = malloc(l_hint + 3 + 4 + 1 + 10);
	strcpy(prefix, hint);
	//strcpy(prefix + l_hint, ".64.bwt");
	strcpy(prefix + l_hint, ".ne.bwt");
	if ((fp = fopen(prefix, "rb")) != 0) {
		fclose(fp);
		prefix[l_hint + 3] = 0;
		return prefix;
	} else {
		strcpy(prefix + l_hint, ".bwt");
		if ((fp = fopen(prefix, "rb")) == 0) {
			free(prefix);
			return 0;
		} else {
			fclose(fp);
			prefix[l_hint] = 0;
			return prefix;
		}
	}
}

bwt_t *bwa_idx_load_bwt(const char *hint)
{
	char *tmp, *prefix;
	bwt_t *bwt;
	prefix = bwa_idx_infer_prefix(hint);
	if (prefix == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to locate the index files\n", __func__);
		return 0;
	}
	tmp = calloc(strlen(prefix) + 5, 1);
	strcat(strcpy(tmp, prefix), ".bwt"); // FM-index
	bwt = bwt_restore_bwt(tmp);
	//strcat(strcpy(tmp, prefix), ".33.4.sa");  // partial suffix array (SA)
	strcat(strcpy(tmp, prefix), ".sa");  // partial suffix array (SA)
	bwt_restore_sa(tmp, bwt);
	free(tmp); free(prefix);
	return bwt;
}

FMTIndex *bwa_idx_load_fmt(const char *hint)
{
	char *fmt_idx_fn, *kmer_idx_fn, *sa_fn;
	// char *kmer_bit_fn;
	FMTIndex *fmt;
	char suffix[32];
	int l_hint = strlen(hint);
	fmt_idx_fn = malloc(l_hint + 32);
	kmer_idx_fn = malloc(l_hint + 32);
	//kmer_bit_fn = malloc(l_hint + 32);
	sa_fn = malloc(l_hint + 32);
//	sprintf(suffix, ".256.%d.fmt", FMT_MID_INTERVAL);
	sprintf(suffix, ".fmt");
	strcpy(fmt_idx_fn, hint);
	strcpy(fmt_idx_fn + l_hint, suffix);
//	sprintf(suffix, ".%d.kmer", HASH_KMER_LEN);
	sprintf(suffix, ".kmer");
	strcpy(kmer_idx_fn, hint);
	strcpy(kmer_idx_fn + l_hint, suffix);

	if (access(fmt_idx_fn, F_OK) != 0 || access(kmer_idx_fn, F_OK) != 0)
	{
		if (bwa_verbose >= 1)
			fprintf(stderr, "[E::%s] fail to locate the index files\n", __func__);
		return 0;
	}
	fmt = fmt_restore_fmt(fmt_idx_fn);
	// fprintf(stderr, "%s\n", kmer_idx_fn);
	fmt->kmer_hash = fmt_restore_kmer_idx(kmer_idx_fn);

	strcpy(sa_fn, hint);
//	sprintf(suffix, ".33.%d.sa", SA_INTV);
	sprintf(suffix, ".bytesa");
	strcpy(sa_fn + l_hint, suffix); // partial suffix array (SA)
	fmt_restore_sa(sa_fn, fmt);

	free(fmt_idx_fn);
	free(kmer_idx_fn);
	free(sa_fn);

	return fmt;
}

bwaidx_t *bwa_idx_load_from_disk(const char *hint, int which)
{
	bwaidx_t *idx;
	char *prefix;
	prefix = bwa_idx_infer_prefix(hint);
	if (prefix == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to locate the index files\n", __func__);
		return 0;
	}
	idx = calloc(1, sizeof(bwaidx_t));


	if (which & BWA_IDX_FMT)  idx->fmt = bwa_idx_load_fmt(hint);
	if (which & BWA_IDX_BWT)
	{
		idx->bwt = bwa_idx_load_bwt(hint);
		if (which & BWA_IDX_FMT)
			idx->bwt->kmer_hash = idx->fmt->kmer_hash;
	}
	
	if (which & BWA_IDX_BNS)
	{
		int i, c;
		idx->bns = bns_restore(prefix);
		for (i = c = 0; i < idx->bns->n_seqs; ++i)
			if (idx->bns->anns[i].is_alt) ++c;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] read %d ALT contigs\n", __func__, c);
		if (which & BWA_IDX_PAC) {
			idx->pac = calloc(idx->bns->l_pac/4+1, 1);
			err_fread_noeof(idx->pac, 1, idx->bns->l_pac/4+1, idx->bns->fp_pac); // concatenated 2-bit encoded sequence
			err_fclose(idx->bns->fp_pac);
			idx->bns->fp_pac = 0;
			// fmtpac
			if (which & BWA_IDX_FMT) {
				idx->fmt->l_pac = idx->bns->l_pac;
				idx->fmt->pac = idx->pac;
			}
		}
	}
	free(prefix);
	return idx;
}

bwaidx_t *bwa_idx_load(const char *hint, int which)
{
	return bwa_idx_load_from_disk(hint, which);
}

void bwa_idx_destroy(bwaidx_t *idx)
{
	if (idx == 0) return;
	if (idx->mem == 0) {
		if (idx->bwt) bwt_destroy(idx->bwt);
		if (idx->bns) bns_destroy(idx->bns);
		if (idx->pac) free(idx->pac);
	} else {
		free(idx->bwt); free(idx->bns->anns); free(idx->bns);
		if (!idx->is_shm) free(idx->mem);
	}
	free(idx);
}

int bwa_mem2idx(int64_t l_mem, uint8_t *mem, bwaidx_t *idx)
{
	int64_t k = 0, x;
	int i;

	// generate idx->bwt
	x = sizeof(bwt_t); idx->bwt = malloc(x); memcpy(idx->bwt, mem + k, x); k += x;
	x = idx->bwt->bwt_size * 4; idx->bwt->bwt = (uint32_t*)(mem + k); k += x;
	x = idx->bwt->n_sa * sizeof(bwtint_t); idx->bwt->sa = (bwtint_t*)(mem + k); k += x;
	
	// generate idx->bns and idx->pac
	x = sizeof(bntseq_t); idx->bns = malloc(x); memcpy(idx->bns, mem + k, x); k += x;
	x = idx->bns->n_holes * sizeof(bntamb1_t); idx->bns->ambs = (bntamb1_t*)(mem + k); k += x;
	x = idx->bns->n_seqs  * sizeof(bntann1_t); idx->bns->anns = malloc(x); memcpy(idx->bns->anns, mem + k, x); k += x;
	for (i = 0; i < idx->bns->n_seqs; ++i) {
		idx->bns->anns[i].name = (char*)(mem + k); k += strlen(idx->bns->anns[i].name) + 1;
		idx->bns->anns[i].anno = (char*)(mem + k); k += strlen(idx->bns->anns[i].anno) + 1;
	}
	idx->pac = (uint8_t*)(mem + k); k += idx->bns->l_pac/4+1;
	assert(k == l_mem);

	idx->l_mem = k; idx->mem = mem;
	return 0;
}

static void mem_to_bnspac(bwaidx_t *idx, uint8_t **mem_, int64_t *k_)
{
	int i;
	int64_t x, k = *k_;
	uint8_t *mem = *mem_;
	// generate idx->bns and idx->pac
	x = sizeof(bntseq_t); idx->bns = malloc(x); memcpy(idx->bns, mem + k, x); k += x;
	x = idx->bns->n_holes * sizeof(bntamb1_t); idx->bns->ambs = (bntamb1_t*)(mem + k); k += x;
	x = idx->bns->n_seqs  * sizeof(bntann1_t); idx->bns->anns = malloc(x); memcpy(idx->bns->anns, mem + k, x); k += x;
	for (i = 0; i < idx->bns->n_seqs; ++i) {
		idx->bns->anns[i].name = (char*)(mem + k); k += strlen(idx->bns->anns[i].name) + 1;
		idx->bns->anns[i].anno = (char*)(mem + k); k += strlen(idx->bns->anns[i].anno) + 1;
	}
	idx->pac = (uint8_t*)(mem + k); k += idx->bns->l_pac/4+1;
	*k_ = k;
	*mem_ = mem;
}

int bwa_mem2fmtidx(int64_t l_mem, uint8_t *mem, bwaidx_t *idx)
{
	int64_t k = 0, x;
	// generate idx->fmt
	x = sizeof(FMTIndex); idx->fmt = malloc(x); memcpy(idx->fmt, mem + k, x); k += x;
	x = idx->fmt->bwt_size * 4; idx->fmt->bwt = (uint32_t*)(mem + k); k += x;
	x = SA_BYTES(idx->fmt->n_sa); idx->fmt->sa = (uint8_t*)(mem + k); k += x;
	// kmer hash
	x = (1 << (10 << 1)) * sizeof(KmerEntryArr);
	idx->fmt->kmer_hash.ke10 = (KmerEntryArr*)(mem + k); k += x;
	x = (1 << (11 << 1)) * sizeof(KmerEntry);
	idx->fmt->kmer_hash.ke11 = (KmerEntry*)(mem + k); k += x;
	x = (1 << (12 << 1)) * sizeof(KmerEntry);
	idx->fmt->kmer_hash.ke12 = (KmerEntry*)(mem + k); k += x;
#if HASH_KMER_LEN > 12
	x = (1 << (13 << 1)) * sizeof(KmerEntry);
	idx->fmt->kmer_hash.ke13 = (KmerEntry*)(mem + k); k += x;
#endif
#if HASH_KMER_LEN > 13
	x = (1 << (14 << 1)) * sizeof(KmerEntry);
	idx->fmt->kmer_hash.ke14 = (KmerEntry*)(mem + k); k += x;
#endif

	// generate idx->bns and idx->pac
	mem_to_bnspac(idx, &mem, &k);
	assert(k == l_mem);
	idx->l_mem = k; idx->mem = mem;
	idx->fmt->pac = idx->pac;
	return 0;
}

static void move_bns_to_mem(bwaidx_t *idx, uint8_t **mem_, int64_t *k_)
{
	int i;
	int64_t x, tmp, k=*k_;
	uint8_t *mem = *mem_;
	// copy idx->bns
	tmp = idx->bns->n_seqs * sizeof(bntann1_t) + idx->bns->n_holes * sizeof(bntamb1_t);
	for (i = 0; i < idx->bns->n_seqs; ++i) // compute the size of heap-allocated memory
		tmp += strlen(idx->bns->anns[i].name) + strlen(idx->bns->anns[i].anno) + 2;
	mem = realloc(mem, k + sizeof(bntseq_t) + tmp);
	x = sizeof(bntseq_t); memcpy(mem + k, idx->bns, x); k += x;
	x = idx->bns->n_holes * sizeof(bntamb1_t); memcpy(mem + k, idx->bns->ambs, x); k += x;
	free(idx->bns->ambs);
	x = idx->bns->n_seqs * sizeof(bntann1_t); memcpy(mem + k, idx->bns->anns, x); k += x;
	for (i = 0; i < idx->bns->n_seqs; ++i) {
		x = strlen(idx->bns->anns[i].name) + 1; memcpy(mem + k, idx->bns->anns[i].name, x); k += x;
		x = strlen(idx->bns->anns[i].anno) + 1; memcpy(mem + k, idx->bns->anns[i].anno, x); k += x;
		free(idx->bns->anns[i].name); free(idx->bns->anns[i].anno);
	}
	free(idx->bns->anns);
	*k_ = k;
	*mem_ = mem;
}

static void move_pac_to_mem(bwaidx_t *idx, uint8_t **mem_, int64_t *k_)
{
	int64_t x, k = *k_;
	uint8_t *mem = *mem_;
	// copy idx->pac
	x = idx->bns->l_pac/4+1;
	mem = realloc(mem, k + x);
	memcpy(mem + k, idx->pac, x); k += x;
	free(idx->bns); idx->bns = 0;
	free(idx->pac); idx->pac = 0;
	*k_ = k;
	*mem_ = mem;
}

int bwa_fmtidx2mem(bwaidx_t *idx)
{
	int64_t k, x;
	uint8_t *mem;

	// copy idx->fmt
	x = idx->fmt->bwt_size * 4;
	mem = realloc(idx->fmt->bwt, sizeof(FMTIndex) + x); idx->fmt->bwt = 0;
	memmove(mem + sizeof(FMTIndex), mem, x);
	memcpy(mem, idx->fmt, sizeof(FMTIndex)); k = sizeof(FMTIndex) + x;
	x = SA_BYTES(idx->fmt->n_sa); mem = realloc(mem, k + x); memcpy(mem + k, idx->fmt->sa, x); k += x;
	// kmer hash
	x = (1 << (10 << 1)) * sizeof(KmerEntryArr);
	mem = realloc(mem, k + x); memcpy(mem + k, idx->fmt->kmer_hash.ke10, x); k += x;
	x = (1 << (11 << 1)) * sizeof(KmerEntry);
	mem = realloc(mem, k + x); memcpy(mem + k, idx->fmt->kmer_hash.ke11, x); k += x;
	x = (1 << (12 << 1)) * sizeof(KmerEntry);
	mem = realloc(mem, k + x); memcpy(mem + k, idx->fmt->kmer_hash.ke12, x); k += x;
#if HASH_KMER_LEN > 12
	x = (1 << (13 << 1)) * sizeof(KmerEntry);
	mem = realloc(mem, k + x); memcpy(mem + k, idx->fmt->kmer_hash.ke13, x); k += x;
#endif
#if HASH_KMER_LEN > 13
	x = (1 << (14 << 1)) * sizeof(KmerEntry);
	mem = realloc(mem, k + x); memcpy(mem + k, idx->fmt->kmer_hash.ke14, x); k += x;
#endif
	free(idx->fmt->kmer_hash.ke10);
	free(idx->fmt->kmer_hash.ke11);
	free(idx->fmt->kmer_hash.ke12);
	free(idx->fmt->kmer_hash.ke13);
	free(idx->fmt->kmer_hash.ke14);
	free(idx->fmt->sa);
	free(idx->fmt); idx->fmt = 0;

	// copy idx->bns
	move_bns_to_mem(idx, &mem, &k);
	// copy idx->pac
	move_pac_to_mem(idx, &mem, &k);

	return bwa_mem2fmtidx(k, mem, idx);
}

int bwa_idx2mem(bwaidx_t *idx)
{
	int i;
	int64_t k, x, tmp;
	uint8_t *mem;

	// copy idx->bwt

	x = idx->bwt->bwt_size * 4;
	mem = realloc(idx->bwt->bwt, sizeof(bwt_t) + x); idx->bwt->bwt = 0;
	memmove(mem + sizeof(bwt_t), mem, x);
	memcpy(mem, idx->bwt, sizeof(bwt_t)); k = sizeof(bwt_t) + x;
	x = idx->bwt->n_sa * sizeof(bwtint_t); mem = realloc(mem, k + x); memcpy(mem + k, idx->bwt->sa, x); k += x;
	free(idx->bwt->sa);
	free(idx->bwt); idx->bwt = 0;


	// copy idx->bns
	tmp = idx->bns->n_seqs * sizeof(bntann1_t) + idx->bns->n_holes * sizeof(bntamb1_t);
	for (i = 0; i < idx->bns->n_seqs; ++i) // compute the size of heap-allocated memory
		tmp += strlen(idx->bns->anns[i].name) + strlen(idx->bns->anns[i].anno) + 2;
	mem = realloc(mem, k + sizeof(bntseq_t) + tmp);
	x = sizeof(bntseq_t); memcpy(mem + k, idx->bns, x); k += x;
	x = idx->bns->n_holes * sizeof(bntamb1_t); memcpy(mem + k, idx->bns->ambs, x); k += x;
	free(idx->bns->ambs);
	x = idx->bns->n_seqs * sizeof(bntann1_t); memcpy(mem + k, idx->bns->anns, x); k += x;
	for (i = 0; i < idx->bns->n_seqs; ++i) {
		x = strlen(idx->bns->anns[i].name) + 1; memcpy(mem + k, idx->bns->anns[i].name, x); k += x;
		x = strlen(idx->bns->anns[i].anno) + 1; memcpy(mem + k, idx->bns->anns[i].anno, x); k += x;
		free(idx->bns->anns[i].name); free(idx->bns->anns[i].anno);
	}
	free(idx->bns->anns);

	// copy idx->pac
	x = idx->bns->l_pac/4+1;
	mem = realloc(mem, k + x);
	memcpy(mem + k, idx->pac, x); k += x;
	free(idx->bns); idx->bns = 0;
	free(idx->pac); idx->pac = 0;

	return bwa_mem2idx(k, mem, idx);
}

/***********************
 * SAM header routines *
 ***********************/

void bwa_print_sam_hdr(const bntseq_t *bns, const char *hdr_line)
{
	int i, n_HD = 0, n_SQ = 0;
	extern char *bwa_pg;
	
	if (hdr_line) {
		// check for HD line
		const char *p = hdr_line;
		if ((p = strstr(p, "@HD")) != 0) {
			++n_HD;
		}	
		// check for SQ lines
		p = hdr_line;
		while ((p = strstr(p, "@SQ\t")) != 0) {
			if (p == hdr_line || *(p-1) == '\n') ++n_SQ;
			p += 4;
		}
	}
	if (n_SQ == 0) {
		for (i = 0; i < bns->n_seqs; ++i) {
			err_printf("@SQ\tSN:%s\tLN:%d", bns->anns[i].name, bns->anns[i].len);
			if (bns->anns[i].is_alt) err_printf("\tAH:*\n");
			else err_fputc('\n', stdout);
		}
	} else if (n_SQ != bns->n_seqs && bwa_verbose >= 2)
		fprintf(stderr, "[W::%s] %d @SQ lines provided with -H; %d sequences in the index. Continue anyway.\n", __func__, n_SQ, bns->n_seqs);
	if (n_HD == 0) {
		err_printf("@HD\tVN:1.5\tSO:unsorted\tGO:query\n");
	}
	if (hdr_line) err_printf("%s\n", hdr_line);
	if (bwa_pg) err_printf("%s\n", bwa_pg);
}

static char *bwa_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == 'n') *q++ = '\n';
			else if (*p == 'r') *q++ = '\r';
			else if (*p == '\\') *q++ = '\\';
		} else *q++ = *p;
	}
	*q = '\0';
	return s;
}

char *bwa_set_rg(const char *s)
{
	char *p, *q, *r, *rg_line = 0;
	memset(bwa_rg_id, 0, 256);
	if (strstr(s, "@RG") != s) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] the read group line is not started with @RG\n", __func__);
		goto err_set_rg;
	}
	if (strstr(s, "\t") != NULL) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] the read group line contained literal <tab> characters -- replace with escaped tabs: \\t\n", __func__);
		goto err_set_rg;
	}
	rg_line = strdup(s);
	bwa_escape(rg_line);
	if ((p = strstr(rg_line, "\tID:")) == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] no ID within the read group line\n", __func__);
		goto err_set_rg;
	}
	p += 4;
	for (q = p; *q && *q != '\t' && *q != '\n'; ++q);
	if (q - p + 1 > 256) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] @RG:ID is longer than 255 characters\n", __func__);
		goto err_set_rg;
	}
	for (q = p, r = bwa_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
		*r++ = *q;
	return rg_line;

err_set_rg:
	free(rg_line);
	return 0;
}

char *bwa_insert_header(const char *s, char *hdr)
{
	int len = 0;
	if (s == 0 || s[0] != '@') return hdr;
	if (hdr) {
		len = strlen(hdr);
		hdr = realloc(hdr, len + strlen(s) + 2);
		hdr[len++] = '\n';
		strcpy(hdr + len, s);
	} else hdr = strdup(s);
	bwa_escape(hdr + len);
	return hdr;
}
