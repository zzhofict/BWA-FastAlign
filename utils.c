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
#define FSYNC_ON_FLUSH

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <errno.h>
#ifdef FSYNC_ON_FLUSH
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <pthread.h>
#include <sys/resource.h>
#include <sys/time.h>

#include "ksort.h"
#include "kvec.h"
#include "utils.h"
#include "yarn.h"
#include "khash.h"
#define pair64_lt(a, b) ((a).x < (b).x || ((a).x == (b).x && (a).y < (b).y))
KSORT_INIT(128, pair64_t, pair64_lt)
KSORT_INIT(64,  uint64_t, ks_lt_generic)

#include "kseq.h"
KSEQ_INIT2(, gzFile, err_gzread)

/********************
 * System utilities *
 ********************/

FILE *err_xopen_core(const char *func, const char *fn, const char *mode)
{
	FILE *fp = 0;
	if (strcmp(fn, "-") == 0)
		return (strstr(mode, "r"))? stdin : stdout;
	if ((fp = fopen(fn, mode)) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
	}
	return fp;
}

FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp)
{
	if (freopen(fn, mode, fp) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
	}
	return fp;
}

gzFile err_xzopen_core(const char *func, const char *fn, const char *mode)
{
	gzFile fp;
	if (strcmp(fn, "-") == 0) {
		fp = gzdopen(fileno((strstr(mode, "r"))? stdin : stdout), mode);
		/* According to zlib.h, this is the only reason gzdopen can fail */
		if (!fp) err_fatal(func, "Out of memory");
		return fp;
	}
	if ((fp = gzopen(fn, mode)) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, errno ? strerror(errno) : "Out of memory");
	}
	return fp;
}

void err_fatal(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, "\n");
	va_end(args);
	exit(EXIT_FAILURE);
}

void err_fatal_core(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, " Abort!\n");
	va_end(args);
	abort();
}

void _err_fatal_simple(const char *func, const char *msg)
{
	fprintf(stderr, "[%s] %s\n", func, msg);
	exit(EXIT_FAILURE);
}

void _err_fatal_simple_core(const char *func, const char *msg)
{
	fprintf(stderr, "[%s] %s Abort!\n", func, msg);
	abort();
}

size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fwrite(ptr, size, nmemb, stream);
	if (ret != nmemb) 
		_err_fatal_simple("fwrite", strerror(errno));
	return ret;
}

size_t err_fread_noeof(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fread(ptr, size, nmemb, stream);
	if (ret != nmemb)
	{
		_err_fatal_simple("fread", ferror(stream) ? strerror(errno) : "Unexpected end of file");
	}
	return ret;
}

uint64_t fread_fix(FILE *fp, uint64_t size, void *a)
{ // Mac/Darwin has a bug when reading data longer than 2GB. This function fixes this issue by reading data in small chunks
	// const int bufsize = 0x1000000; // 16M block
	const int bufsize = 0x4000000; // 64M block
	uint64_t offset = 0;
	while (size)
	{
		int x = bufsize < size ? bufsize : size;
		if ((x = err_fread_noeof(a + offset, 1, x, fp)) == 0)
			break;
		size -= x;
		offset += x;
	}
	return offset;
}

typedef struct {
    pthread_t tid;
    void *buf[2];
    volatile int readSize[2];
    uint64_t getIdx;
    uint64_t putIdx;
    volatile int finish;
    lock_t *mtx;
} FileKV;

KHASH_MAP_INIT_INT64(fkv, FileKV);
static khash_t(fkv) *fHash = 0;

#define USE_ASYNC_READ

int err_gzread(gzFile file, void *ptr, unsigned int len)
{
    int ret = 0;
    PROF_START(read);
#ifndef USE_ASYNC_READ
    ret = gzread(file, ptr, len);
#else
    khiter_t k = kh_get(fkv, fHash, (int64_t)file);
    FileKV *val = &kh_value(fHash, k);
    POSSESS(val->mtx);
    WAIT_FOR(val->mtx, NOT_TO_BE, 0);  // 
    RELEASE(val->mtx);

    int curIdx = val->getIdx % 2;
    if (val->finish) {
        if (val->getIdx < val->putIdx) {
            ret = val->readSize[curIdx];
			if (ret > 0) memcpy(ptr, val->buf[curIdx], ret);
            ++val->getIdx;
            return ret;
        }
        return 0;
    }
    ret = val->readSize[curIdx];
    memcpy(ptr, val->buf[curIdx], ret);

    POSSESS(val->mtx);
    ++val->getIdx;
    TWIST(val->mtx, BY, -1);
#endif
    PROF_END(gprof[G_read_seq], read);

    if (ret < 0)
	{
		int errnum = 0;
		const char *msg = gzerror(file, &errnum);
		_err_fatal_simple("gzread", Z_ERRNO == errnum ? strerror(errno) : msg);
	}

	return ret;
}

static int64_t kBufSize = 16777216;

static void *async_gzread(void *data) {
    gzFile file = (gzFile)data;
    khiter_t k = kh_get(fkv, fHash, (int64_t)file);
    FileKV *val = &kh_value(fHash, k);

    int ret = 0;
    while (1) {
        POSSESS(val->mtx);
        WAIT_FOR(val->mtx, NOT_TO_BE, 2);  // 
        RELEASE(val->mtx);

        int curIdx = val->putIdx % 2;
        ret = gzread(file, val->buf[curIdx], kBufSize);
        val->readSize[curIdx] = ret;

        if (ret <= 0) {
            POSSESS(val->mtx);
            val->finish = 1;
            TWIST(val->mtx, BY, 1);
            break;
        }

        POSSESS(val->mtx);
        val->putIdx += 1;
        TWIST(val->mtx, BY, 1);
    }

    return NULL;
}

int start_async_read(gzFile file) {
	int ret = 0;
#ifdef USE_ASYNC_READ
    if (fHash == 0) {
        fHash = kh_init(fkv);
    }
    khiter_t k = kh_put(fkv, fHash, (int64_t)file, &ret);
    kh_key(fHash, k) = (int64_t)file;
    FileKV *fv =  &kh_value(fHash, k);

    fv->mtx = NEW_LOCK(0);
    fv->getIdx = fv->putIdx = fv->finish = 0;
    fv->readSize[0] = fv->readSize[1] = 0;
    fv->buf[0] = malloc(kBufSize);
    fv->buf[1] = malloc(kBufSize);
    ret = pthread_create(&fv->tid, 0, async_gzread, file);
#endif
    return ret;
}

int stop_async_read(gzFile file) {
#ifdef USE_ASYNC_READ
    khiter_t k = kh_get(fkv, fHash, (int64_t)file);
    FileKV *val = &kh_value(fHash, k);
    pthread_join(val->tid, 0);
#endif
    return 0;
}

int err_fseek(FILE *stream, long offset, int whence)
{
	int ret = fseek(stream, offset, whence);
	if (0 != ret)
	{
		_err_fatal_simple("fseek", strerror(errno));
	}
	return ret;
}

long err_ftell(FILE *stream)
{
	long ret = ftell(stream);
	if (-1 == ret)
	{
		_err_fatal_simple("ftell", strerror(errno));
	}
	return ret;
}

int err_printf(const char *format, ...) 
{
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stdout, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf(stdout)", strerror(saveErrno));
	return done;
}

int err_fprintf(FILE *stream, const char *format, ...) 
{
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stream, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf", strerror(saveErrno));
	return done;
}

int err_fputc(int c, FILE *stream)
{
	int ret = putc(c, stream);
	if (EOF == ret)
	{
		_err_fatal_simple("fputc", strerror(errno));
	}

	return ret;
}

int err_fputs(const char *s, FILE *stream)
{
	int ret = fputs(s, stream);
	if (EOF == ret)
	{
		_err_fatal_simple("fputs", strerror(errno));
	}

	return ret;
}

int err_puts(const char *s)
{
	int ret = puts(s);
	if (EOF == ret)
	{
		_err_fatal_simple("puts", strerror(errno));
	}

	return ret;
}

int err_fflush(FILE *stream) 
{
    int ret = fflush(stream);
    if (ret != 0) _err_fatal_simple("fflush", strerror(errno));

#ifdef FSYNC_ON_FLUSH
	/* Calling fflush() ensures that all the data has made it to the
	   kernel buffers, but this may not be sufficient for remote filesystems
	   (e.g. NFS, lustre) as an error may still occur while the kernel
	   is copying the buffered data to the file server.  To be sure of
	   catching these errors, we need to call fsync() on the file
	   descriptor, but only if it is a regular file.  */
	{
		struct stat sbuf;
		if (0 != fstat(fileno(stream), &sbuf))
			_err_fatal_simple("fstat", strerror(errno));
		
		if (S_ISREG(sbuf.st_mode))
		{
			if (0 != fsync(fileno(stream)))
				_err_fatal_simple("fsync", strerror(errno));
		}
	}
#endif
    return ret;
}

int err_fclose(FILE *stream) 
{
	int ret = fclose(stream);
	if (ret != 0) _err_fatal_simple("fclose", strerror(errno));
	return ret;
}

int err_gzclose(gzFile file)
{
	int ret = gzclose(file);
	if (Z_OK != ret)
	{
		_err_fatal_simple("gzclose", Z_ERRNO == ret ? strerror(errno) : zError(ret));
	}

	return ret;
}

/*********
 * Timer *
 *********/

double cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime(void)
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

uint64_t realtime_msec(void)
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (uint64_t)1000 * (tv.tv_sec + ((1e-6) * tv.tv_usec));
}

long peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

int memcpy_bwamem(void *dest, size_t dmax, const void *src, size_t smax, char *file_name, int line_num)
{
#define RSIZE_MAX_MEM (256UL << 20) /* 256MB */
	if (dmax < smax)
	{
		fprintf(stderr, "[%s: %d] src size is lager than dest size.(src: %ld; dest: %ld)\n", file_name, line_num, smax, dmax);
		exit(EXIT_FAILURE);
	}
	int64_t bytes_copied;
	for (bytes_copied = 0; bytes_copied < smax; bytes_copied += RSIZE_MAX_MEM)
	{
		int64_t bytes_remaining = smax - bytes_copied;
		int64_t bytes_to_copy = (bytes_remaining > RSIZE_MAX_MEM) ? RSIZE_MAX_MEM : bytes_remaining;
		memcpy((char *)dest + bytes_copied, (const char *)src + bytes_copied, bytes_to_copy);
	}
	return 0;
}