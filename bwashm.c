#include <sys/types.h>
#include <sys/mman.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <stdio.h>
#include "bwa.h"

int bwa_shm_stage(bwaidx_t *idx, const char *hint, int useERT)
{

	const char *name;
	uint8_t *shm, *shm_idx;
	uint16_t *cnt;
	int shmid, to_init = 0, l;
	char path[PATH_MAX + 1], *tmpfn = 0;

	if (hint == 0 || hint[0] == 0) return -1;
	for (name = hint + strlen(hint) - 1; name >= hint && *name != '/'; --name);
	++name;

	if ((shmid = shm_open("/bwactl", O_RDWR, 0)) < 0) {
		shmid = shm_open("/bwactl", O_CREAT|O_RDWR|O_EXCL, 0644);
		to_init = 1;
	}
	if (shmid < 0) return -1;
	if (ftruncate(shmid, BWA_CTL_SIZE) < 0) return -1;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ|PROT_WRITE, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	if (to_init) {
		memset(shm, 0, BWA_CTL_SIZE);
		cnt[1] = 4;
	}

	// 
	if (idx->mem == 0) {
		bwa_fmtidx2mem(idx);
	}

	if (tmpfn) {
		FILE *fp;
		if ((fp = fopen(tmpfn, "wb")) != 0) {
			int64_t rest = idx->l_mem;
			while (rest > 0) {
				int64_t l = rest < 0x1000000? rest : 0x1000000;
				rest -= fwrite(&idx->mem[idx->l_mem - rest], 1, l, fp);
			}
			fclose(fp);
			free(idx->mem); idx->mem = 0;
		} else {
			fprintf(stderr, "[W::%s] fail to create the temporary file. Option '-f' is ignored.\n", __func__);
			tmpfn = 0;
		}
	}

	strcat(strcpy(path, "/bwaidx-"), name);
	if ((shmid = shm_open(path, O_CREAT|O_RDWR|O_EXCL, 0644)) < 0) {
		shm_unlink(path);
		perror("shm_open()");
		return -1;
	}
	l = 8 + strlen(name) + 1;
	if (cnt[1] + l > BWA_CTL_SIZE) return -1;
	memcpy(shm + cnt[1], &idx->l_mem, 8);
	memcpy(shm + cnt[1] + 8, name, l - 8);
	cnt[1] += l; ++cnt[0];
	if (ftruncate(shmid, idx->l_mem) < 0) return -1;
	shm_idx = mmap(0, idx->l_mem, PROT_READ|PROT_WRITE, MAP_SHARED, shmid, 0);
	if (tmpfn) {
		FILE *fp;
		fp = fopen(tmpfn, "rb");
		int64_t rest = idx->l_mem;
		while (rest > 0) {
			int64_t l = rest < 0x1000000? rest : 0x1000000;
			rest -= fread(&shm_idx[idx->l_mem - rest], 1, l, fp);
		}
		fclose(fp);
		unlink(tmpfn);
	} else {
		memcpy(shm_idx, idx->mem, idx->l_mem);
		free(idx->mem);
	}
	bwa_mem2fmtidx(idx->l_mem, shm_idx, idx);
	idx->is_shm = 1;
	return 0;
}

#define INIT_SHM_LOAD \
	const char *name; \
	uint8_t *shm, *shm_idx; \
	uint16_t *cnt, i; \
	char *p, path[PATH_MAX + 1]; \
	int shmid; \
	int64_t l_mem; \
	bwaidx_t *idx; \
	if (hint == 0 || hint[0] == 0) return 0; \
	for (name = hint + strlen(hint) - 1; name >= hint && *name != '/'; --name); \
	++name; \
	if ((shmid = shm_open("/bwactl", O_RDONLY, 0)) < 0) return 0; \
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ, MAP_SHARED, shmid, 0); \
	cnt = (uint16_t*)shm; \
	if (cnt[0] == 0) return 0; \
	for (i = 0, p = (char*)(shm + 4); i < cnt[0]; ++i) { \
		memcpy(&l_mem, p, 8); p += 8; \
		if (strcmp(p, name) == 0) break; \
		p += strlen(p) + 1; \
	} \
	if (i == cnt[0]) return 0; \
	strcat(strcpy(path, "/bwaidx-"), name); \
	if ((shmid = shm_open(path, O_RDONLY, 0)) < 0) return 0; \
	shm_idx = mmap(0, l_mem, PROT_READ, MAP_SHARED, shmid, 0); \
	idx = calloc(1, sizeof(bwaidx_t));


bwaidx_t *bwa_fmtidx_load_from_shm(const char *hint)
{
	INIT_SHM_LOAD;
	bwa_mem2fmtidx(l_mem, shm_idx, idx);
	idx->is_shm = 1;
	return idx;
}

bwaidx_t *bwa_idx_load_from_shm(const char *hint)
{
	const char *name;
	uint8_t *shm, *shm_idx;
	uint16_t *cnt, i;
	char *p, path[PATH_MAX + 1];
	int shmid;
	int64_t l_mem;
	bwaidx_t *idx;

	if (hint == 0 || hint[0] == 0) return 0;
	for (name = hint + strlen(hint) - 1; name >= hint && *name != '/'; --name);
	++name;
	if ((shmid = shm_open("/bwactl", O_RDONLY, 0)) < 0) return 0;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	if (cnt[0] == 0) return 0;
	for (i = 0, p = (char*)(shm + 4); i < cnt[0]; ++i) {
		memcpy(&l_mem, p, 8); p += 8;
		if (strcmp(p, name) == 0) break;
		p += strlen(p) + 1;
	}
	if (i == cnt[0]) return 0;

	strcat(strcpy(path, "/bwaidx-"), name);
	if ((shmid = shm_open(path, O_RDONLY, 0)) < 0) return 0;
	shm_idx = mmap(0, l_mem, PROT_READ, MAP_SHARED, shmid, 0);
	idx = calloc(1, sizeof(bwaidx_t));
	bwa_mem2idx(l_mem, shm_idx, idx);
	idx->is_shm = 1;
	return idx;
}

int bwa_shm_test(const char *hint)
{
	int shmid;
	uint16_t *cnt, i;
	char *p, *shm;
	const char *name;

	if (hint == 0 || hint[0] == 0) return 0;
	for (name = hint + strlen(hint) - 1; name >= hint && *name != '/'; --name);
	++name;
	if ((shmid = shm_open("/bwactl", O_RDONLY, 0)) < 0) return 0;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	for (i = 0, p = shm + 4; i < cnt[0]; ++i) {
		if (strcmp(p + 8, name) == 0) return 1;
		p += strlen(p) + 9;
	}
	return 0;
}

int bwa_shm_list(void)
{
	int shmid;
	uint16_t *cnt, i;
	char *p, *shm;
	if ((shmid = shm_open("/bwactl", O_RDONLY, 0)) < 0) return -1;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	for (i = 0, p = shm + 4; i < cnt[0]; ++i) {
		int64_t l_mem;
		memcpy(&l_mem, p, 8); p += 8;
		printf("%s\t%ld\n", p, (long)l_mem);
		p += strlen(p) + 1;
	}
	return 0;
}

int bwa_shm_destroy(void)
{
	int shmid;
	uint16_t *cnt, i;
	char *p, *shm;
	char path[PATH_MAX + 1];

	if ((shmid = shm_open("/bwactl", O_RDONLY, 0)) < 0) return -1;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	for (i = 0, p = shm + 4; i < cnt[0]; ++i) {
		int64_t l_mem;
		memcpy(&l_mem, p, 8); p += 8;
		strcat(strcpy(path, "/bwaidx-"), p);
		shm_unlink(path);
		p += strlen(p) + 1;
	}
	munmap(shm, BWA_CTL_SIZE);
	shm_unlink("/bwactl");
	return 0;
}

int bwa_shm_stage_fmt(const char *idx_prefix)
{
	const char *name;
	uint8_t *shm, *mem;
	uint16_t *cnt;
	int shmid, to_init = 0, l;
	char path[PATH_MAX];
	int64_t l_mem = 0, x = 0, k = 0;
	char fn[PATH_MAX];
	FILE *fp;
	FMTIndex fmt;
	bntseq_t *bns;
	int i;
	// clac l_mem
	// fmt
	x = sizeof(FMTIndex); l_mem += x;
	sprintf(fn, "%s.fmt", idx_prefix);
	fp = xopen(fn, "rb");
	err_fseek(fp, 0, SEEK_END);
	x = ftell(fp) - sizeof(bwtint_t) * 6 - 3; l_mem += x;
	fmt.bwt_size = x >> 2;
	fseek(fp, 0, SEEK_SET);
	err_fread_noeof(&fmt.primary, sizeof(bwtint_t), 1, fp);
	err_fread_noeof(&fmt.sec_primary, sizeof(bwtint_t), 1, fp);
	err_fread_noeof(&fmt.sec_bcp, sizeof(uint8_t), 1, fp);
	err_fread_noeof(&fmt.first_base, sizeof(uint8_t), 1, fp);
	err_fread_noeof(&fmt.last_base, sizeof(uint8_t), 1, fp);
	err_fread_noeof(fmt.L2 + 1, sizeof(bwtint_t), 4, fp);
	fmt.seq_len = fmt.L2[4];
	fmt_gen_cnt_occ(&fmt);
	err_fclose(fp);
	// sa
	sprintf(fn, "%s.bytesa", idx_prefix);
	fp = xopen(fn, "rb");
	err_fseek(fp, sizeof(bwtint_t) * 5, SEEK_SET);
	err_fread_noeof(&fmt.sa_intv, sizeof(bwtint_t), 1, fp);
	fmt.n_sa = (fmt.seq_len + fmt.sa_intv) / fmt.sa_intv;
	x = SA_BYTES(fmt.n_sa); l_mem += x;
	err_fclose(fp);
	// kmaer hash
	x = (1 << (10 << 1)) * sizeof(KmerEntryArr)
	  + (1 << (11 << 1)) * sizeof(KmerEntry)
	  + (1 << (12 << 1)) * sizeof(KmerEntry);
#if HASH_KMER_LEN > 12
	x += (1 << (13 << 1)) * sizeof(KmerEntry);
#endif
#if HASH_KMER_LEN > 13
	x += (1 << (14 << 1)) * sizeof(KmerEntry);
#endif
	l_mem += x;
	// bns
	x = sizeof(bntseq_t); l_mem += x;
	bns = bns_restore(idx_prefix);
	x = bns->n_holes * sizeof(bntamb1_t); l_mem += x;
	x = bns->n_seqs * sizeof(bntann1_t); l_mem += x;
	for (i = 0; i < bns->n_seqs; ++i) {
		x = strlen(bns->anns[i].name) + 1; l_mem += x;
		x = strlen(bns->anns[i].anno) + 1; l_mem += x;
	}
	// pac
	x = bns->l_pac / 4 + 1; l_mem += x;
	fmt.l_pac = bns->l_pac;

	fprintf(stderr, "l_mem: %ld\n", l_mem);

	for (name = idx_prefix + strlen(idx_prefix) - 1; name >= idx_prefix && *name != '/'; --name) ;
	++name;
	if ((shmid = shm_open("/bwactl", O_RDWR, 0)) < 0) {
		shmid = shm_open("/bwactl", O_CREAT|O_RDWR|O_EXCL, 0644);
		to_init = 1;
	}
	if (shmid < 0) return -1;
	if (ftruncate(shmid, BWA_CTL_SIZE) < 0) return -1;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ|PROT_WRITE, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	if (to_init) {
		memset(shm, 0, BWA_CTL_SIZE);
		cnt[1] = 4;
	}
	strcat(strcpy(path, "/bwaidx-"), name);
	if ((shmid = shm_open(path, O_CREAT|O_RDWR|O_EXCL, 0644)) < 0) {
		shm_unlink(path);
		perror("shm_open()");
		return -1;
	}
	l = 8 + strlen(name) + 1;
	if (cnt[1] + l > BWA_CTL_SIZE) return -1;
	memcpy(shm + cnt[1], &l_mem, 8);
	memcpy(shm + cnt[1] + 8, name, l - 8);
	cnt[1] += l; ++cnt[0];
	if (ftruncate(shmid, l_mem) < 0) return -1;
	mem = mmap(0, l_mem, PROT_READ|PROT_WRITE, MAP_SHARED, shmid, 0);

	// write to share mem
	
	// fmt
	x = sizeof(FMTIndex);
	memcpy(mem, &fmt, x); k = x;
	sprintf(fn, "%s.fmt", idx_prefix);
	fp = xopen(fn, "rb");
	err_fseek(fp, sizeof(bwtint_t) * 6 + 3, SEEK_SET);
	x = fmt.bwt_size * 4;
	fread_fix(fp, x, mem + k); k += x;
	err_fclose(fp);
	// sa
	sprintf(fn, "%s.bytesa", idx_prefix);
	fp = xopen(fn, "rb");
	err_fseek(fp, sizeof(bwtint_t) * 7, SEEK_SET);
	x = SA_BYTES(fmt.n_sa);
	fread_fix(fp, x, mem + k); k += x;
	err_fclose(fp);
	// kmer hash
	sprintf(fn, "%s.kmer", idx_prefix);
	fp = xopen(fn, "rb");
	x = (1 << (10 << 1)) * sizeof(KmerEntryArr);
	fread_fix(fp, x, mem + k); k += x;
	x = (1 << (11 << 1)) * sizeof(KmerEntry);
	fread_fix(fp, x, mem + k); k += x;
	x = (1 << (12 << 1)) * sizeof(KmerEntry);
	fread_fix(fp, x, mem + k); k += x;
#if HASH_KMER_LEN > 12
	x = (1 << (13 << 1)) * sizeof(KmerEntry);
	fread_fix(fp, x, mem + k); k += x;
#endif
#if HASH_KMER_LEN > 13
	x = (1 << (14 << 1)) * sizeof(KmerEntry);
	fread_fix(fp, x, mem + k); k += x;
#endif
	err_fclose(fp);
	// bns
	x = sizeof(bntseq_t); memcpy(mem + k, bns, x); k += x;
	x = bns->n_holes * sizeof(bntamb1_t); memcpy(mem + k, bns->ambs, x); k += x;
	x = bns->n_seqs * sizeof(bntann1_t); memcpy(mem + k, bns->anns, x); k += x;
		for (i = 0; i < bns->n_seqs; ++i) {
		x = strlen(bns->anns[i].name) + 1; memcpy(mem + k, bns->anns[i].name, x); k += x;
		x = strlen(bns->anns[i].anno) + 1; memcpy(mem + k, bns->anns[i].anno, x); k += x;
	}
	// pac
	x = bns->l_pac / 4 + 1;
	err_fread_noeof(mem + k, 1, x, bns->fp_pac); k += x;// concatenated 2-bit encoded sequence
	err_fclose(bns->fp_pac);

	fprintf(stderr, "k: %ld\n", k);

	return 0;
}

int main_shm(int argc, char *argv[])
{
	int c, to_list = 0, to_drop = 0, ret = 0;
	char shm_prefix[PATH_MAX];
	while ((c = getopt(argc, argv, "ldZ")) >= 0)
	{
		if (c == 'l') to_list = 1;
		else if (c == 'd') to_drop = 1;
	}
	if (optind == argc && !to_list && !to_drop) {
		fprintf(stderr, "\nUsage: fastbwa shm [-d|-l] [-f tmpFile] [idxbase]\n\n");
		fprintf(stderr, "Options: -d       destroy all indices in shared memory\n");
		fprintf(stderr, "         -l       list names of indices in shared memory\n");
		return 1;
	}
	if (optind < argc && (to_list || to_drop)) {
		fprintf(stderr, "[E::%s] open -l or -d cannot be used when 'idxbase' is present\n", __func__);
		return 1;
	}
	if (optind < argc) {
		sprintf(shm_prefix, "%s", argv[optind]);
		if (bwa_shm_test(shm_prefix) == 0)
		{
#if 0
			bwaidx_t *idx;
			if (useERT)
				idx = bwa_ertidx_load_from_disk(argv[optind]);
			else
				idx = bwa_idx_load_from_disk(argv[optind], BWA_IDX_BNS | BWA_IDX_PAC | BWA_IDX_FMT);
			if (bwa_shm_stage(idx, shm_prefix, useERT) < 0) {
				fprintf(stderr, "[E::%s] failed to stage the index in shared memory\n", __func__);
				ret = 1;
			}
			bwa_idx_destroy(idx);
#else

			if (bwa_shm_stage_fmt(argv[optind]) < 0) {
				fprintf(stderr, "[E::%s] failed to stage the index in shared memory\n", __func__);
				ret = 1;
			}
			
#endif
		}
		else
			fprintf(stderr, "[M::%s] index '%s' is already in shared memory\n", __func__, argv[optind]);
	}
	if (to_list) bwa_shm_list();
	if (to_drop) bwa_shm_destroy();
	return ret;
}
