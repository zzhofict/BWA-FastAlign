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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>
#include "bntseq.h"
#include "bwa.h"
#include "bwt.h"
#include "utils.h"
#include "rle.h"
#include "rope.h"

#ifdef _DIVBWT
#include "divsufsort.h"
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif


int is_bwt(ubyte_t *T, int n);

int64_t bwa_seq_len(const char *fn_pac)
{
	FILE *fp;
	int64_t pac_len;
	ubyte_t c;
	fp = xopen(fn_pac, "rb");
	err_fseek(fp, -1, SEEK_END);
	pac_len = err_ftell(fp);
	err_fread_noeof(&c, 1, 1, fp);
	err_fclose(fp);
	return (pac_len - 1) * 4 + (int)c;
}

bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is)
{
	bwt_t *bwt;
	ubyte_t *buf, *buf2;
	int64_t i, pac_size;
	FILE *fp;

	// initialization
	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	bwt->seq_len = bwa_seq_len(fn_pac);
	bwt->bwt_size = (bwt->seq_len + 15) >> 4;
	fp = xopen(fn_pac, "rb");

	// prepare sequence
	pac_size = (bwt->seq_len>>2) + ((bwt->seq_len&3) == 0? 0 : 1);
	buf2 = (ubyte_t*)calloc(pac_size, 1);
	err_fread_noeof(buf2, 1, pac_size, fp);
	err_fclose(fp);
	memset(bwt->L2, 0, 5 * 4);
	buf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
	for (i = 0; i < bwt->seq_len; ++i) {
		buf[i] = buf2[i>>2] >> ((3 - (i&3)) << 1) & 3;
		++bwt->L2[1+buf[i]];
	}
	for (i = 2; i <= 4; ++i) bwt->L2[i] += bwt->L2[i-1];
	free(buf2);

	// Burrows-Wheeler Transform
	if (use_is) {
		bwt->primary = is_bwt(buf, bwt->seq_len);
	} else {
		rope_t *r;
		int64_t x;
		rpitr_t itr;
		const uint8_t *blk;

		r = rope_init(ROPE_DEF_MAX_NODES, ROPE_DEF_BLOCK_LEN);
		for (i = bwt->seq_len - 1, x = 0; i >= 0; --i) {
			int c = buf[i] + 1;
			x = rope_insert_run(r, x, c, 1, 0) + 1;
			while (--c >= 0) x += r->c[c];
		}
		bwt->primary = x;
		rope_itr_first(r, &itr);
		x = 0;
		while ((blk = rope_itr_next_block(&itr)) != 0) {
			const uint8_t *q = blk + 2, *end = blk + 2 + *rle_nptr(blk);
			while (q < end) {
				int c = 0;
				int64_t l;
				rle_dec1(q, c, l);
				for (i = 0; i < l; ++i)
					buf[x++] = c - 1;
			}
		}
		rope_destroy(r);
	}
	bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
	for (i = 0; i < bwt->seq_len; ++i)
		bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
	free(buf);
	return bwt;
}

int bwa_pac2bwt(int argc, char *argv[]) // the "pac2bwt" command; IMPORTANT: bwt generated at this step CANNOT be used with BWA. bwtupdate is required!
{
	bwt_t *bwt;
	int c, use_is = 1;
	while ((c = getopt(argc, argv, "d")) >= 0) {
		switch (c) {
		case 'd': use_is = 0; break;
		default: return 1;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: fastbwa pac2bwt [-d] <in.pac> <out.bwt>\n");
		return 1;
	}
	bwt = bwt_pac2bwt(argv[optind], use_is);
	bwt_dump_bwt(argv[optind+1], bwt);
	bwt_destroy(bwt);
	return 0;
}

int bwa_pac2bref(int argc, char *argv[])
{
	FILE *fp_pac, *fp_bref;
	uint64_t pac_size, pac_len, seq_len, bref_size;
	int64_t i;
	uint8_t *pac_buf, *bref_buf;
	uint8_t c;
	char pac_file_name[PATH_MAX];
	char binary_ref_name[PATH_MAX];
	if (optind + 1 > argc)
	{
		fprintf(stderr, "Usage: fastbwa pac2ref <index_prefix>\n");
		return 1;
	}
	sprintf(pac_file_name, "%s.pac", argv[optind]);
	sprintf(binary_ref_name, "%s.ert.0123", argv[optind]);

	fp_pac = xopen(pac_file_name, "rb");
	err_fseek(fp_pac, -1, SEEK_END);
	pac_len = err_ftell(fp_pac);
	err_fread_noeof(&c, 1, 1, fp_pac);
	seq_len = (pac_len - 1) * 4 + (int)c;

	err_fseek(fp_pac, 0, SEEK_SET);
	pac_size = (seq_len >> 2) + ((seq_len & 3) == 0 ? 0 : 1);
	pac_buf = (uint8_t *)calloc(pac_size, 1);
	err_fread_noeof(pac_buf, 1, pac_size, fp_pac);
	err_fclose(fp_pac);

	bref_size = seq_len * 2;
	bref_buf = (uint8_t *)malloc(bref_size);
	for (i = 0; i < seq_len; ++i) {
		int nt = pac_buf[i >> 2] >> ((3 - (i & 3)) << 1) & 3;
		bref_buf[i] = nt;
	}

	for (i = seq_len - 1; i >= 0; --i) {
		int64_t j = bref_size - i - 1;
		bref_buf[j] = 3 - bref_buf[i];
	}

	fp_bref = xopen(binary_ref_name, "wb");
	err_fwrite(bref_buf, sizeof(uint8_t), bref_size, fp_bref);

	free(bref_buf);
	free(pac_buf);
	err_fclose(fp_bref);

	return 0;
}

#define bwt_B00(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

void bwt_bwtupdate_core(bwt_t *bwt)
{
	bwtint_t i, k, c[4], n_occ;
	uint32_t *buf;

	n_occ = (bwt->seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
	bwt->bwt_size += n_occ * sizeof(bwtint_t); // the new size
	buf = (uint32_t*)calloc(bwt->bwt_size, 4); // will be the new bwt
	c[0] = c[1] = c[2] = c[3] = 0;
	for (i = k = 0; i < bwt->seq_len; ++i) {
		if (i % OCC_INTERVAL == 0) {
			memcpy(buf + k, c, sizeof(bwtint_t) * 4);
			k += sizeof(bwtint_t); // in fact: sizeof(bwtint_t)=4*(sizeof(bwtint_t)/4)
		}
		if (i % 16 == 0) buf[k++] = bwt->bwt[i/16]; // 16 == sizeof(uint32_t)/2
		++c[bwt_B00(bwt, i)];
	}
	// the last element
	memcpy(buf + k, c, sizeof(bwtint_t) * 4);
	xassert(k + sizeof(bwtint_t) == bwt->bwt_size, "inconsistent bwt_size");
	// update bwt
	free(bwt->bwt); bwt->bwt = buf;
}

int bwa_bwtupdate(int argc, char *argv[]) // the "bwtupdate" command
{
	bwt_t *bwt;
	if (argc != 2) {
		fprintf(stderr, "Usage: fastbwa bwtupdate <the.bwt>\n");
		return 1;
	}
	bwt = bwt_restore_bwt(argv[1]);
	bwt_bwtupdate_core(bwt);
	bwt_dump_bwt(argv[1], bwt);
	bwt_destroy(bwt);
	return 0;
}

int bwa_bwt2sa(int argc, char *argv[]) // the "bwt2sa" command
{
	bwt_t *bwt;
	int c, sa_intv = 32;
	while ((c = getopt(argc, argv, "i:")) >= 0) {
		switch (c) {
		case 'i': sa_intv = atoi(optarg); break;
		default: return 1;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: fastbwa bwt2sa [-i %d] <in.bwt> <out.sa>\n", sa_intv);
		return 1;
	}
	bwt = bwt_restore_bwt(argv[optind]);
	bwt_cal_sa(bwt, sa_intv);
	bwt_dump_sa(argv[optind+1], bwt);
	bwt_destroy(bwt);
	return 0;
}

int bwa_bwt2bytesa(int argc, char *argv[]) // the "bwt2bytesa" command
{
	bwt_t *bwt;
	int c, sa_intv = 32;
	while ((c = getopt(argc, argv, "i:")) >= 0) {
		switch (c) {
		case 'i': sa_intv = atoi(optarg); break;
		default: return 1;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: fastbwa bwt2bytesa [-i %d] <in.bwt> <out.sa>\n", sa_intv);
		return 1;
	}
	bwt = bwt_restore_bwt(argv[optind]);
	bwt_cal_byte_sa(bwt, sa_intv);
	bwt_dump_byte_sa(argv[optind+1], bwt);
	bwt_destroy(bwt);
	return 0;
}

int bwa_bwt2fmt(int argc, char *argv[]) // create fmt index
{
	bwt_t *bwt;
	//char buf[1024];
	if (optind + 1 > argc) {
		fprintf(stderr, "Usage: fastbwa bwt2fmt <in.bwt> <out.fmt>\n");
		return 1;
	}
	bwt = bwt_restore_bwt(argv[optind]);
	FMTIndex *fmt;
	fmt = create_fmt_from_bwt(bwt);
	dump_fmt(argv[optind + 1], fmt);
	// sprintf(buf, "%s.kmer", argv[optind + 1]);
	return 0;
}

int bwa_build_kmer(int argc, char *argv[])
{
	char buf[1024];
	if (optind + 1 > argc)
	{
		fprintf(stderr, "Usage: fastbwa build_kmerhash <in.fmt> <out.kmerhash>\n");
		return 1;
	}
	FMTIndex *fmt = fmt_restore_fmt(argv[optind]);
	fmt_create_kmer_index(fmt);
	sprintf(buf, "%s", argv[optind + 1]);
	fmt_dump_kmer_idx(buf, &fmt->kmer_hash);
	return 0;
}

int bwa_index(int argc, char *argv[]) // the "index" command
{
	int num_threads = 1;
	int c, algo_type = BWTALGO_AUTO, is_64 = 0, block_size = 10000000;
	char *prefix = 0, *str;
	while ((c = getopt(argc, argv, "6a:p:b:")) >= 0) {
		switch (c) {
		case 'a': // if -a is not set, algo_type will be determined later
			if (strcmp(optarg, "rb2") == 0) algo_type = BWTALGO_RB2;
			else if (strcmp(optarg, "bwtsw") == 0) algo_type = BWTALGO_BWTSW;
			else if (strcmp(optarg, "is") == 0) algo_type = BWTALGO_IS;
			else err_fatal(__func__, "unknown algorithm: '%s'.", optarg);
			break;
		case 'p': prefix = strdup(optarg); break;
		case '6': is_64 = 1; break;
		case 't':
			num_threads = atoi(optarg);
			assert(num_threads > 0 && num_threads < 256);
			break;
		case 'b':
			block_size = strtol(optarg, &str, 10);
			if (*str == 'G' || *str == 'g') block_size *= 1024 * 1024 * 1024;
			else if (*str == 'M' || *str == 'm') block_size *= 1024 * 1024;
			else if (*str == 'K' || *str == 'k') block_size *= 1024;
			break;
		default: return 1;
		}
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fastbwa index [options] <in.fasta>\n\n");
		fprintf(stderr, "Options: -a STR    BWT construction algorithm: bwtsw, is or rb2 [auto]\n");
		fprintf(stderr, "         -p STR    prefix of the index [same as fasta name]\n");
		fprintf(stderr, "         -t INT    number of threads for ERT index building [%d]\n", num_threads);
		fprintf(stderr, "         -b INT    block size for the bwtsw algorithm (effective with -a bwtsw) [%d]\n", block_size);
		fprintf(stderr, "         -6        index files named as <in.fasta>.64.* instead of <in.fasta>.* \n");
		fprintf(stderr, "\n");
		fprintf(stderr,	"Warning: `-a bwtsw' does not work for short genomes, while `-a is' and\n");
		fprintf(stderr, "         `-a div' do not work not for long genomes.\n\n");
		return 1;
	}
	if (prefix == 0) {
		prefix = malloc(strlen(argv[optind]) + 4);
		strcpy(prefix, argv[optind]);
		if (is_64) strcat(prefix, ".64");
	}
	bwa_idx_build(argv[optind], prefix, algo_type, block_size);
	free(prefix);
	return 0;
}

int bwa_idx_build(const char *fa, const char *prefix, int algo_type, int block_size)
{
	extern void bwa_pac_rev_core(const char *fn, const char *fn_rev);

	char *str, *str2, *str3;
	clock_t t;
	int64_t l_pac;

	str  = (char*)calloc(strlen(prefix) + 10, 1);
	str2 = (char*)calloc(strlen(prefix) + 10, 1);
	str3 = (char*)calloc(strlen(prefix) + 10, 1);

	{ // nucleotide indexing
		gzFile fp = xzopen(fa, "r");
		t = clock();
		if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Pack FASTA... ");
		l_pac = bns_fasta2bntseq(fp, prefix, 0);
		if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		err_gzclose(fp);
	}
	if (algo_type == 0) algo_type = l_pac > 50000000? 2 : 3; // set the algorithm for generating BWT
	{
		strcpy(str, prefix); strcat(str, ".pac");
		strcpy(str2, prefix); strcat(str2, ".bwt");
		t = clock();
		if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Construct BWT for the packed sequence...\n");
		if (algo_type == 2) bwt_bwtgen2(str, str2, block_size);
		else if (algo_type == 1 || algo_type == 3) {
			bwt_t *bwt;
			bwt = bwt_pac2bwt(str, algo_type == 3);
			bwt_dump_bwt(str2, bwt);
			bwt_destroy(bwt);
		}
		if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] %.2f seconds elapse.\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	{
		bwt_t *bwt;
		strcpy(str, prefix); strcat(str, ".bwt");
		t = clock();
		if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Update BWT... ");
		bwt = bwt_restore_bwt(str);
		bwt_bwtupdate_core(bwt);
		bwt_dump_bwt(str, bwt);
		bwt_destroy(bwt);
		if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	{
		gzFile fp = xzopen(fa, "r");
		t = clock();
		if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Pack forward-only FASTA... ");
		l_pac = bns_fasta2bntseq(fp, prefix, 1);
		if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		err_gzclose(fp);
	}
	{
		bwt_t *bwt;
		strcpy(str, prefix); strcat(str, ".bwt");
		strcpy(str3, prefix); strcat(str3, ".sa");
		t = clock();
		if (bwa_verbose >= 3) fprintf(stderr, "[bwa_index] Construct SA from BWT and Occ... ");
		bwt = bwt_restore_bwt(str);
		bwt_cal_sa(bwt, 32);
		bwt_dump_sa(str3, bwt);
//		bwt_destroy(bwt);
		if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		{
			// build FMT-Index
			t = clock();
			FMTIndex *fmt;
			strcpy(str, prefix); strcat(str, ".fmt");
			fmt = create_fmt_from_bwt(bwt);
			dump_fmt(str, fmt);
			if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
			// create Kmer-Hash
			t = clock();
			fmt_create_kmer_index(fmt);
			strcpy(str, prefix); strcat(str, ".kmer");
			fmt_dump_kmer_idx(str, &fmt->kmer_hash);
			if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
			// create byte sa
			t = clock();
			bwt_cal_byte_sa(bwt, 4);
			strcpy(str, prefix); strcat(str, ".bytesa");
			bwt_dump_byte_sa(str, bwt);
			if (bwa_verbose >= 3) fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		}
	}
	free(str3); free(str2); free(str);
	return 0;
}
