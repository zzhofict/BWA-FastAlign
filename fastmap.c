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
#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <pthread.h>
#include <immintrin.h>
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "bntseq.h"
#include "kseq.h"
#include "fmt_idx.h"
#include "yarn.h"
#include "profiling.h"
#include "debug.h"
KSEQ_DECLARE(gzFile)

extern unsigned char nst_nt4_table[256];

void *kopen(const char *fn, int *_fd);
int kclose(void *a);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
	int n_seqs;
	int n_sams;
	int m_seqs;
	int m_sams;
	bseq1_t *seqs;
	seq_sam_t *sams;
} ktp_data_t;

typedef struct {
	kseq_t *ks, *ks2;
	mem_opt_t *opt;
	mem_pestat_t *pes0;
	int64_t n_processed;
	int copy_comment, actual_chunk_size;
	bwaidx_t *idx;
	mem_worker_t *w;
	int data_idx; // pingpong buffer index
	ktp_data_t *data;
	int wbuf_size;
	char *wbuf;
	volatile int read_complete;
	volatile int calc_complete;
	long read_idx;
	long calc_idx;
	long write_idx;
} ktp_aux_t;

// read
static inline void *read_data(ktp_aux_t *aux, ktp_data_t *data)
{

	PROF_START(read);
	ktp_data_t *ret = aux->data + aux->data_idx;
	aux->data_idx = !aux->data_idx;
	int64_t size = 0;
	bseq_read(aux->actual_chunk_size, &ret->n_seqs, aux->ks, aux->ks2, aux->copy_comment, &size, &ret->m_seqs, &ret->seqs);
	PROF_END(gprof[G_READ], read);

#if 0
	// 
	int i;
	int64_t ms = 0;
	for (i = 0; i < ret->n_seqs; ++i)
	{
		ms += ret->seqs[i].m_name + ret->seqs[i].m_seq + ret->seqs[i].m_comment + ret->seqs[i].m_qual;
	}
	fprintf(stderr, "seq size: %ld M\n", ms / 1024 / 1024);
#endif

	if (ret->n_seqs == 0) {
		return 0;
	}
	if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] read %d sequences (%ld bp)...\n", __func__, ret->n_seqs, (long)size);
	return ret;
}

// calculate
static inline void *calc_data(ktp_aux_t *aux, ktp_data_t *data)
{
	PROF_START(compute);
	const mem_opt_t *opt = aux->opt;
	if (data->n_sams != data->n_seqs) {
		if (data->m_sams < data->m_seqs) {
			data->m_sams = data->m_seqs;
			data->sams = realloc(data->sams, data->m_sams * sizeof(seq_sam_t));
			memset(data->sams + data->n_sams, 0, (data->m_sams - data->n_sams) * sizeof(seq_sam_t));
		}
		data->n_sams = data->n_seqs;
	}
	if (opt->flag & MEM_F_SMARTPE) {
		int i;
		fprintf(stderr, "here MEM_F_SMARTPE\n");
		bseq1_t *sep[2];
		seq_sam_t *ss[2];
		int n_sep[2];
		mem_opt_t tmp_opt = *opt;
		bseq_classify(data->n_seqs, data->seqs, n_sep, sep);
		ss[0] = calloc(0, n_sep[0] * sizeof(seq_sam_t));
		ss[1] = calloc(0, n_sep[1] * sizeof(seq_sam_t));
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] %d single-end sequences; %d paired-end sequences\n", __func__, n_sep[0], n_sep[1]);
		if (n_sep[0]) {
			tmp_opt.flag &= ~MEM_F_PE;
			mem_process_seqs(&tmp_opt, aux->w, aux->n_processed, n_sep[0], sep[0], 0, ss[0]);
			for (i = 0; i < n_sep[0]; ++i)
				data->sams[sep[0][i].id].sam = ss[0][i].sam;
		}
		if (n_sep[1]) {
			tmp_opt.flag |= MEM_F_PE;
			mem_process_seqs(&tmp_opt, aux->w, aux->n_processed + n_sep[0], n_sep[1], sep[1], aux->pes0, ss[1]);
			for (i = 0; i < n_sep[1]; ++i)
				data->sams[sep[1][i].id].sam = ss[1][i].sam;
		}
		free(sep[0]); free(sep[1]);
		free(ss[0]); free(ss[1]);
	} else mem_process_seqs(opt, aux->w, aux->n_processed, data->n_seqs, data->seqs, aux->pes0, data->sams);
	aux->n_processed += data->n_seqs;
	PROF_END(gprof[G_COMPUTE], compute);

	return data;
}

// write
static inline void *write_data(ktp_aux_t *aux, ktp_data_t *data)
{
	int i;
	PROF_START(write);
	// int64_t ms = 0;
	int buf_written = 0;
	for (i = 0; i < data->n_sams; ++i)
	{
		const int slen = data->sams[i].sam.l;
		if (slen && (buf_written + slen) < aux->wbuf_size)
		{
			memcpy(&aux->wbuf[buf_written], data->sams[i].sam.s, slen);
			buf_written += slen;
		}
		else if (buf_written > 0)
		{
			err_fwrite(aux->wbuf, 1, buf_written, stdout);
			if ((buf_written + slen) >= aux->wbuf_size)
			{
				memcpy(&aux->wbuf[0], data->sams[i].sam.s, slen);
				buf_written = slen;
			}
			else
			{
				buf_written = 0;
			}
		}
	}
	if (buf_written > 0) {
		err_fwrite(aux->wbuf, 1, buf_written, stdout);
	}

	//for (i = 0; i < data->n_sams; ++i)
	//{
	//	//ms += data->sams[i].sam.m;
	//	if (data->sams[i].sam.l)
	//		err_fputs(data->sams[i].sam.s, stdout);
	//}

	//fprintf(stderr, "sam size: %ld M\n", ms / 1024 / 1024);
	PROF_END(gprof[G_WRITE], write);
	return 0;
}

// io ，
static void *process(void *shared, int step, void *_data)
{
	ktp_aux_t *aux = (ktp_aux_t*)shared;
	ktp_data_t *data = (ktp_data_t*)_data;
	if (step == 0) {
		return read_data(aux, data);
	} else if (step == 1) {
		return calc_data(aux, data);
	} else if (step == 2) {
		return write_data(aux, data);
	}
	return 0;
}

////////////// pipeline
static lock_t *input_have = NULL;
static lock_t *output_have = NULL;

static void *thread_read(void *data)
{
	ktp_aux_t *aux = (ktp_aux_t *)data;
	while (1)
	{
        PROF_START(w1);
        POSSESS(input_have);
        WAIT_FOR(input_have, NOT_TO_BE, 0);
		RELEASE(input_have);
        PROF_END(gprof[G_read_wait_1], w1);
        // fprintf(stderr, "start read: %ld\n", aux->read_idx);
		if (read_data(aux, aux->data) == 0) {
			POSSESS(input_have);
			aux->read_complete = 1;
			TWIST(input_have, BY, -1);
			break;
		}
        PROF_START(w2);
        POSSESS(input_have);
		aux->read_idx++;
		TWIST(input_have, BY, -1);
        PROF_END(gprof[G_read_wait_2], w2);
        // fprintf(stderr, "next read: %ld\n", aux->read_idx);
	}
	return 0;
}

static void *thread_calc(void *data)
{
	ktp_aux_t *aux = (ktp_aux_t *)data;
	int d_idx = 0;
	int add_idx = 0;
	while (1)
	{
		// fprintf(stderr, "start calc: %ld\n", aux->calc_idx);
        PROF_START(w1);
        POSSESS(input_have);
		WAIT_FOR(input_have, NOT_TO_BE, 2);
		RELEASE(input_have); // ？

		POSSESS(output_have);
		WAIT_FOR(output_have, NOT_TO_BE, 2);
		RELEASE(output_have);
        PROF_END(gprof[G_comp_wait_1], w1);

        if (aux->calc_idx < aux->read_idx) {
			calc_data(aux, aux->data + d_idx);
			d_idx = !d_idx;
			add_idx = 1;
		}
		if (aux->read_complete)
		{
			POSSESS(output_have);
			if (add_idx) aux->calc_idx ++;
			aux->calc_complete = 1;
			TWIST(output_have, BY, 1); // 
			break;					   // 
		}
        PROF_START(w2);
        POSSESS(output_have);
		if (add_idx) aux->calc_idx ++;
		TWIST(output_have, BY, 1);

		POSSESS(input_have);
		TWIST(input_have, BY, 1);
        PROF_END(gprof[G_comp_wait_2], w2);
        // fprintf(stderr, "next calc: %ld\n", aux->calc_idx);
	}
	return 0;
}

static void *thread_write(void *data)
{
	ktp_aux_t *aux = (ktp_aux_t *)data;
	int d_idx = 0;
	while (1)
	{
		// fprintf(stderr, "start write: %ld\n", aux->write_idx);
        PROF_START(w1);
        POSSESS(output_have);
		WAIT_FOR(output_have, NOT_TO_BE, 0);
		RELEASE(output_have);
        PROF_END(gprof[G_write_wait_1], w1);
        if (aux->write_idx < aux->calc_idx) {
			write_data(aux, aux->data + d_idx);
			d_idx = !d_idx;
			aux->write_idx++;
		}
		if (aux->calc_complete) {
			if (aux->write_idx < aux->calc_idx)
				write_data(aux, aux->data + d_idx);
			break;
		}
        PROF_START(w2);
        POSSESS(output_have);
		TWIST(output_have, BY, -1);
        PROF_END(gprof[G_write_wait_2], w2);
        // fprintf(stderr, "next write: %ld\n", aux->write_idx);
	}
	return 0;
}

static void new_pipeline(ktp_aux_t *aux)
{
	input_have = NEW_LOCK(2);
	output_have = NEW_LOCK(0);
	pthread_t tid[3];
	int i;

	pthread_create(&tid[0], 0, thread_read, aux);
	pthread_create(&tid[1], 0, thread_calc, aux);
	pthread_create(&tid[2], 0, thread_write, aux);

	for (i = 0; i < 3; ++i) pthread_join(tid[i], 0);
}

//////////////

static void update_a(mem_opt_t *opt, const mem_opt_t *opt0)
{
	if (opt0->a) { // matching score is changed
		if (!opt0->b) opt->b *= opt->a;
		if (!opt0->T) opt->T *= opt->a;
		if (!opt0->o_del) opt->o_del *= opt->a;
		if (!opt0->e_del) opt->e_del *= opt->a;
		if (!opt0->o_ins) opt->o_ins *= opt->a;
		if (!opt0->e_ins) opt->e_ins *= opt->a;
		if (!opt0->zdrop) opt->zdrop *= opt->a;
		if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
		if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
		if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
	}
}

int main_mem(int argc, char *argv[])
{
#ifdef DEBUG_FILE_OUTPUT
	open_debug_files();
	// open_debug_files
#endif

#ifdef SHOW_PERF
#if USE_RDTSC
	uint64_t tmp_time = __rdtsc();
	sleep(1);
	proc_freq = __rdtsc() - tmp_time;
#else
	proc_freq = 1000;
#endif
#endif

	PROF_START(all);
	PROF_START(prepare);
	mem_opt_t *opt, opt0;
	int fd, fd2, i, c, ignore_alt = 0, no_mt_io = 0;
	int fixed_chunk_size = -1;
	gzFile fp, fp2 = 0;
	char *p, *rg_line = 0, *hdr_line = 0;
	const char *mode = 0;
	void *ko = 0, *ko2 = 0;
	mem_pestat_t pes[4];
	ktp_aux_t aux;
	int useERT = 0;

	memset(&aux, 0, sizeof(ktp_aux_t));
	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	for (i = 0; i < 4; ++i) pes[i].failed = 1;

	aux.opt = opt = mem_opt_init();
	memset(&opt0, 0, sizeof(mem_opt_t));
	while ((c = getopt(argc, argv, "512qpaMCSPVYjuk:c:v:s:r:t:b:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:o:f:W:x:G:h:y:K:X:H:F:z:Z")) >= 0) {
		if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
		else if (c == '1') no_mt_io = 1;
		else if (c == '2') no_mt_io = 2;
		else if (c == 'x') mode = optarg;
		else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1;
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1;
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1;
		else if (c == 'U') opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1;
		else if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'b') opt->batch_size = atoi(optarg) >> 1 << 1, opt->batch_size = opt->batch_size > 1? opt->batch_size : 512;
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else if (c == 'p') opt->flag |= MEM_F_PE | MEM_F_SMARTPE;
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
		else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
		else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
		else if (c == '5') opt->flag |= MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ; // always apply MEM_F_KEEP_SUPP_MAPQ with -5
		else if (c == 'q') opt->flag |= MEM_F_KEEP_SUPP_MAPQ;
		else if (c == 'u') opt->flag |= MEM_F_XB;
		else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
		else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'j') ignore_alt = 1;
		else if (c == 'r') opt->split_factor = atof(optarg), opt0.split_factor = 1.;
		else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
		else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
		else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
		else if (c == 'G') opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
		else if (c == 'N') opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
		else if (c == 'o' || c == 'f') xreopen(optarg, "wb", stdout);
		else if (c == 'W') opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
		else if (c == 'y') opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;
		else if (c == 'C') aux.copy_comment = 1;
		else if (c == 'K') fixed_chunk_size = atoi(optarg);
		else if (c == 'X') opt->mask_level = atof(optarg);
		else if (c == 'F') bwa_dbg = atoi(optarg);
		else if (c == 'h') {
			opt0.max_XA_hits = opt0.max_XA_hits_alt = 1;
			opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->max_XA_hits_alt = strtol(p+1, &p, 10);
		}
		else if (c == 'z') opt->XA_drop_ratio = atof(optarg);
		else if (c == 'Q') {
			opt0.mapQ_coef_len = 1;
			opt->mapQ_coef_len = atoi(optarg);
			opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
		} else if (c == 'O') {
			opt0.o_del = opt0.o_ins = 1;
			opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->o_ins = strtol(p+1, &p, 10);
		} else if (c == 'E') {
			opt0.e_del = opt0.e_ins = 1;
			opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->e_ins = strtol(p+1, &p, 10);
		} else if (c == 'L') {
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
			opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->pen_clip3 = strtol(p+1, &p, 10);
		} else if (c == 'R') {
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1; // FIXME: memory leak
		} else if (c == 'H') {
			if (optarg[0] != '@') {
				FILE *fp;
				if ((fp = fopen(optarg, "r")) != 0) {
					char *buf;
					buf = calloc(1, 0x10000);
					while (fgets(buf, 0xffff, fp)) {
						i = strlen(buf);
						assert(buf[i-1] == '\n'); // a long line
						buf[i-1] = 0;
						hdr_line = bwa_insert_header(buf, hdr_line);
					}
					free(buf);
					fclose(fp);
				}
			} else hdr_line = bwa_insert_header(optarg, hdr_line);
		} else if (c == 'I') { // specify the insert size distribution
			aux.pes0 = pes;
			pes[1].failed = 0;
			pes[1].avg = strtod(optarg, &p);
			pes[1].std = pes[1].avg * .1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].std = strtod(p+1, &p);
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
			if (pes[1].low < 1) pes[1].low = 1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].high = (int)(strtod(p+1, &p) + .499);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].low  = (int)(strtod(p+1, &p) + .499);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] mean insert size: %.3f, stddev: %.3f, max: %d, min: %d\n",
						__func__, pes[1].avg, pes[1].std, pes[1].high, pes[1].low);
		}
		else return 1;
	}

	if (rg_line) {
		hdr_line = bwa_insert_header(rg_line, hdr_line);
		free(rg_line);
	}

	if (opt->n_threads < 1) opt->n_threads = 1;
	if (opt->batch_size < 1) opt->batch_size = 256;
	// fprintf(stderr, "batch_size: %d\n", opt->batch_size);

	if (optind + 1 >= argc || optind + 3 < argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: fastbwa mem [options] <idxbase> <in1.fq> [in2.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT        number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -b INT        batch size of reads to process at one time [%d]\n", opt->batch_size);
		fprintf(stderr, "       -k INT        minimum seed length [%d]\n", opt->min_seed_len);
		fprintf(stderr, "       -w INT        band width for banded alignment [%d]\n", opt->w);
		fprintf(stderr, "       -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
		fprintf(stderr, "       -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
		fprintf(stderr, "       -y INT        seed occurrence for the 3rd round seeding [%ld]\n", (long)opt->max_mem_intv);
//		fprintf(stderr, "       -s INT        look for internal seeds inside a seed with less than INT occ [%d]\n", opt->split_width);
		fprintf(stderr, "       -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
		fprintf(stderr, "       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->drop_ratio);
		fprintf(stderr, "       -W INT        discard a chain if seeded bases shorter than INT [0]\n");
		fprintf(stderr, "       -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
		fprintf(stderr, "       -S            skip mate rescue\n");
		fprintf(stderr, "       -2            read and write in different disks \n");
		fprintf(stderr, "       -P            skip pairing; mate rescue performed unless -S also in use\n");
		fprintf(stderr, "\nScoring options:\n\n");
		fprintf(stderr, "       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [%d]\n", opt->a);
		fprintf(stderr, "       -B INT        penalty for a mismatch [%d]\n", opt->b);
		fprintf(stderr, "       -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
		fprintf(stderr, "       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", opt->e_del, opt->e_ins);
		fprintf(stderr, "       -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
		fprintf(stderr, "       -U INT        penalty for an unpaired read pair [%d]\n\n", opt->pen_unpaired);
		fprintf(stderr, "       -x STR        read type. Setting -x changes multiple parameters unless overridden [null]\n");
		fprintf(stderr, "                     pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)\n");
		fprintf(stderr, "                     ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)\n");
		fprintf(stderr, "                     intractg: -B9 -O16 -L5  (intra-species contigs to ref)\n");
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p            smart pairing (ignoring in2.fq)\n");
		fprintf(stderr, "       -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
		fprintf(stderr, "       -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]\n");
		fprintf(stderr, "       -o FILE       sam file to output results to [stdout]\n");
		fprintf(stderr, "       -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n");
		fprintf(stderr, "       -5            for split alignment, take the alignment with the smallest query (not genomic) coordinate as primary\n");
		fprintf(stderr, "       -q            don't modify mapQ of supplementary alignments\n");
		fprintf(stderr, "       -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility) []\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT        verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "       -T INT        minimum score to output [%d]\n", opt->T);
		fprintf(stderr, "       -h INT[,INT]  if there are <INT hits with score >%.2f%% of the max score, output all in XA [%d,%d]\n", 
				opt->XA_drop_ratio * 100.0,
				opt->max_XA_hits, opt->max_XA_hits_alt);
		fprintf(stderr, "                     A second value may be given for alternate sequences.\n");
		fprintf(stderr, "       -z FLOAT      The fraction of the max score to use with -h [%f].\n", opt->XA_drop_ratio);
		fprintf(stderr, "                     specify the mean, standard deviation (10%% of the mean if absent), max\n");
		fprintf(stderr, "       -a            output all alignments for SE or unpaired PE\n");
		fprintf(stderr, "       -C            append FASTA/FASTQ comment to SAM output\n");
		fprintf(stderr, "       -V            output the reference FASTA header in the XR tag\n");
		fprintf(stderr, "       -Y            use soft clipping for supplementary alignments\n");
		fprintf(stderr, "       -M            mark shorter split hits as secondary\n\n");
		fprintf(stderr, "       -I FLOAT[,FLOAT[,INT[,INT]]]\n");
		fprintf(stderr, "                     specify the mean, standard deviation (10%% of the mean if absent), max\n");
		fprintf(stderr, "                     (4 sigma from the mean if absent) and min of the insert size distribution.\n");
		fprintf(stderr, "                     FR orientation only. [inferred]\n");
		fprintf(stderr, "       -u            output XB instead of XA; XB is XA with the alignment score and mapping quality added.\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		free(opt);
		return 1;
	}

	if (mode) {
		if (strcmp(mode, "intractg") == 0) {
			if (!opt0.o_del) opt->o_del = 16;
			if (!opt0.o_ins) opt->o_ins = 16;
			if (!opt0.b) opt->b = 9;
			if (!opt0.pen_clip5) opt->pen_clip5 = 5;
			if (!opt0.pen_clip3) opt->pen_clip3 = 5;
		} else if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "ont2d") == 0) {
			if (!opt0.o_del) opt->o_del = 1;
			if (!opt0.e_del) opt->e_del = 1;
			if (!opt0.o_ins) opt->o_ins = 1;
			if (!opt0.e_ins) opt->e_ins = 1;
			if (!opt0.b) opt->b = 1;
			if (opt0.split_factor == 0.) opt->split_factor = 10.;
			if (strcmp(mode, "ont2d") == 0) {
				if (!opt0.min_chain_weight) opt->min_chain_weight = 20;
				if (!opt0.min_seed_len) opt->min_seed_len = 14;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			} else {
				if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
				if (!opt0.min_seed_len) opt->min_seed_len = 17;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			}
		} else {
			fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);
			return 1; // FIXME memory leak
		}
	} else update_a(opt, &opt0);
	bwa_fill_scmat(opt->a, opt->b, opt->mat);
	PROF_END(gprof[G_PREPARE], prepare);

	PROF_START(idx);
	if (useERT) aux.idx = bwa_ertidx_load_from_shm(argv[optind]);
	else aux.idx = bwa_fmtidx_load_from_shm(argv[optind]);
	if (aux.idx == 0) {
		if (!useERT) {
			if ((aux.idx = bwa_idx_load(argv[optind], BWA_IDX_BNS | BWA_IDX_PAC | BWA_IDX_FMT)) == 0) return 1; // FIXME: memory leak
		}
		else {
			if ((aux.idx = bwa_ertidx_load_from_disk(argv[optind])) == 0) return 1; // FIXME: memory leak
		}

	} else if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] load the bwa index from shared memory\n", __func__);
	if (ignore_alt)
		for (i = 0; i < aux.idx->bns->n_seqs; ++i)
			aux.idx->bns->anns[i].is_alt = 0;
	PROF_END(gprof[G_LOAD_IDX], idx);

	ko = kopen(argv[optind + 1], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
		return 1;
	}
	fp = gzdopen(fd, "r");
    start_async_read(fp);
    // stop_async_read(fp);
    aux.ks = kseq_init(fp);
    if (optind + 2 < argc) {
		if (opt->flag&MEM_F_PE) {
			if (bwa_verbose >= 2)
				fprintf(stderr, "[W::%s] when '-p' is in use, the second query file is ignored.\n", __func__);
		} else {
			ko2 = kopen(argv[optind + 2], &fd2);
			if (ko2 == 0) {
				if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 2]);
				return 1;
			}
			fp2 = gzdopen(fd2, "r");
            start_async_read(fp2);
            aux.ks2 = kseq_init(fp2);
			opt->flag |= MEM_F_PE;
		}
	}

	//fprintf(stderr, "%ld %ld %ld %ld %ld\n", aux.idx->fmt->L2[0], aux.idx->fmt->L2[1], aux.idx->fmt->L2[2], aux.idx->fmt->L2[3], aux.idx->fmt->L2[4]);
	aux.w = init_mem_worker(opt, aux.idx->fmt, aux.idx->ert, aux.idx->bns, aux.idx->pac, useERT);
	aux.data = calloc(2, sizeof(ktp_data_t));

	bwa_print_sam_hdr(aux.idx->bns, hdr_line);
	aux.actual_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads;

	// allocate write buffer
	aux.wbuf_size = 16777216;
	aux.wbuf = malloc(aux.wbuf_size);

	PROF_START(pipeline);

	if (no_mt_io == 2) new_pipeline(&aux);
	else kt_pipeline(no_mt_io? 1 : 2, process, &aux, 3);

	PROF_END(gprof[G_PIPELINE], pipeline);

//	free(hdr_line);
//	free(opt);
//	bwa_idx_destroy(aux.idx);
//	kseq_destroy(aux.ks);
    stop_async_read(fp);
    err_gzclose(fp); kclose(ko);
    if (aux.ks2) {
//		kseq_destroy(aux.ks2);
		stop_async_read(fp2);
		err_gzclose(fp2); kclose(ko2);
	}
	PROF_END(gprof[G_ALL], all);

#ifdef SHOW_PERF
	display_stats(opt->n_threads);
#endif

#ifdef SHOW_DATA_PERF
	fprintf(stderr, "\n");
	// fprintf(stderr, "full_match: %ld, all: %ld\n", gdat[1], gdat[0]);
	fprintf(stderr, "\n");
#endif

#ifdef DEBUG_FILE_OUTPUT
	close_files();
#endif

	return 0;
}

int main_fastmap(int argc, char *argv[])
{
	int c, i, min_iwidth = 20, min_len = 17, print_seq = 0, min_intv = 1, max_len = INT_MAX;
	uint64_t max_intv = 0;
	kseq_t *seq;
	bwtint_t k;
	gzFile fp;
	smem_i *itr;
	const bwtintv_v *a;
	bwaidx_t *idx;

	while ((c = getopt(argc, argv, "w:l:pi:I:L:")) >= 0) {
		switch (c) {
			case 'p': print_seq = 1; break;
			case 'w': min_iwidth = atoi(optarg); break;
			case 'l': min_len = atoi(optarg); break;
			case 'i': min_intv = atoi(optarg); break;
			case 'I': max_intv = atol(optarg); break;
			case 'L': max_len  = atoi(optarg); break;
		    default: return 1;
		}
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fastbwa fastmap [options] <idxbase> <in.fq>\n\n");
		fprintf(stderr, "Options: -l INT    min SMEM length to output [%d]\n", min_len);
		fprintf(stderr, "         -w INT    max interval size to find coordiantes [%d]\n", min_iwidth);
		fprintf(stderr, "         -i INT    min SMEM interval size [%d]\n", min_intv);
		fprintf(stderr, "         -L INT    max MEM length [%d]\n", max_len);
		fprintf(stderr, "         -I INT    stop if MEM is longer than -l with a size less than INT [%ld]\n", (long)max_intv);
		fprintf(stderr, "\n");
		return 1;
	}

	fp = xzopen(argv[optind + 1], "r");
	seq = kseq_init(fp);
	if ((idx = bwa_idx_load(argv[optind], BWA_IDX_BWT|BWA_IDX_BNS)) == 0) return 1;
	itr = smem_itr_init(idx->bwt);
	smem_config(itr, min_intv, max_len, max_intv);
	while (kseq_read(seq) >= 0) {
		err_printf("SQ\t%s\t%ld", seq->name.s, seq->seq.l);
		if (print_seq) {
			err_putchar('\t');
			err_puts(seq->seq.s);
		} else err_putchar('\n');
		for (i = 0; i < seq->seq.l; ++i)
			seq->seq.s[i] = nst_nt4_table[(int)seq->seq.s[i]];
		smem_set_query(itr, seq->seq.l, (uint8_t*)seq->seq.s);
		while ((a = smem_next(itr)) != 0) {
			for (i = 0; i < a->n; ++i) {
				bwtintv_t *p = &a->a[i];
				if ((uint32_t)p->info - (p->info>>32) < min_len) continue;
				err_printf("EM\t%d\t%d\t%ld", (uint32_t)(p->info>>32), (uint32_t)p->info, (long)p->x[2]);
				if (p->x[2] <= min_iwidth) {
					for (k = 0; k < p->x[2]; ++k) {
						bwtint_t pos;
						int len, is_rev, ref_id;
						len  = (uint32_t)p->info - (p->info>>32);
						pos = bns_depos(idx->bns, bwt_sa(idx->bwt, p->x[0] + k), &is_rev);
						if (is_rev) pos -= len - 1;
						bns_cnt_ambi(idx->bns, pos, len, &ref_id);
						err_printf("\t%s:%c%ld", idx->bns->anns[ref_id].name, "+-"[is_rev], (long)(pos - idx->bns->anns[ref_id].offset) + 1);
					}
				} else err_puts("\t*");
				err_putchar('\n');
			}
		}
		err_puts("//");
	}

	smem_itr_destroy(itr);
	bwa_idx_destroy(idx);
	kseq_destroy(seq);
	err_gzclose(fp);
	return 0;
}
