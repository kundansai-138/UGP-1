/*
 * adios2_wait_shim.c - time every adios2_begin_step / adios2_end_step / adios2_close call of a process,
 * without rebuilding it. Loaded with LD_PRELOAD in front of the ADIOS2 C library; the Fortran bindings
 * (Xcompact3d / 2decomp, streaming SPOD) call these C functions, so they pass through here.
 *
 * Build (no ADIOS2 headers or libraries needed, plain gcc):
 *     gcc -O2 -fPIC -shared -o libadios2_wait.so adios2_wait_shim.c -ldl
 *
 * Use (MPICH):
 *     mpirun -genv LD_PRELOAD /path/libadios2_wait.so -genv ADIOS_WAIT_DIR $PWD/wait \
 *            -genv ADIOS_WAIT_TAG reader ... ./streamingSPOD
 *
 * Environment
 *     ADIOS_WAIT_DIR     folder for the CSV files (must exist). Not set -> pass through, log nothing.
 *     ADIOS_WAIT_TAG     file prefix and first column, e.g. reader / writer (default: proc)
 *     ADIOS_WAIT_ENGINE  only log engines whose name contains this text (default: all engines)
 *
 * Output: <dir>/<tag>_r<rank>.csv, one line per call:
 *     tag,rank,host,engine,call,step,t_enter,dur,status,rc
 *     t_enter = CLOCK_REALTIME (epoch s) at entry, dur = CLOCK_MONOTONIC seconds inside the call,
 *     step = adios2_current_step after a successful begin_step (-1 otherwise),
 *     status = begin_step status (0 ok, 1 not ready, 2 end of stream, -1 error; -9 for other calls),
 *     rc = adios2_error returned by the call.
 * Rank from PMI_RANK (MPICH Hydra), OMPI_COMM_WORLD_RANK, PMIX_RANK or SLURM_PROCID, else the pid.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef struct adios2_engine adios2_engine; /* opaque, same as in adios2_c_types.h */
typedef int adios2_error;                   /* C enums are int-sized on x86-64 */

typedef adios2_error (*begin_fn)(adios2_engine *, int, float, int *);
typedef adios2_error (*plain_fn)(adios2_engine *);
typedef adios2_error (*cstep_fn)(size_t *, const adios2_engine *);
typedef adios2_error (*name_fn)(char *, size_t *, const adios2_engine *);

static begin_fn real_begin;
static plain_fn real_end, real_close;
static cstep_fn real_cstep;
static name_fn real_name;

static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
static FILE *out;
static int disabled;          /* 1 = ADIOS_WAIT_DIR not set or file could not be opened */
static char tag[64] = "proc", host[256] = "?", rank[32] = "?", filt[256] = "";

#define NCACHE 64
static struct { const adios2_engine *e; char name[512]; int keep; } cache[NCACHE];

static double now_real(void) { struct timespec t; clock_gettime(CLOCK_REALTIME, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static double now_mono(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static void *need(const char *sym)
{
    void *p = dlsym(RTLD_NEXT, sym);
    if (!p) { fprintf(stderr, "adios2_wait_shim: %s not found (is ADIOS2 loaded?)\n", sym); abort(); }
    return p;
}

static void init_once(void)
{
    static int done;
    if (done) return;
    done = 1;
    real_begin = (begin_fn)need("adios2_begin_step");
    real_end = (plain_fn)need("adios2_end_step");
    real_close = (plain_fn)need("adios2_close");
    real_cstep = (cstep_fn)dlsym(RTLD_NEXT, "adios2_current_step");
    real_name = (name_fn)dlsym(RTLD_NEXT, "adios2_engine_name");

    const char *dir = getenv("ADIOS_WAIT_DIR");
    if (!dir || !*dir) { disabled = 1; return; }
    const char *s;
    if ((s = getenv("ADIOS_WAIT_TAG")) && *s) snprintf(tag, sizeof tag, "%s", s);
    if ((s = getenv("ADIOS_WAIT_ENGINE")) && *s) snprintf(filt, sizeof filt, "%s", s);
    const char *rk[] = {"PMI_RANK", "OMPI_COMM_WORLD_RANK", "PMIX_RANK", "SLURM_PROCID"};
    snprintf(rank, sizeof rank, "p%d", (int)getpid());
    for (int i = 0; i < 4; i++)
        if ((s = getenv(rk[i])) && *s) { snprintf(rank, sizeof rank, "%s", s); break; }
    gethostname(host, sizeof host - 1);

    char path[4096];
    snprintf(path, sizeof path, "%s/%s_r%s.csv", dir, tag, rank);
    out = fopen(path, "w");
    if (!out) { fprintf(stderr, "adios2_wait_shim: cannot write %s\n", path); disabled = 1; return; }
    fprintf(out, "tag,rank,host,engine,call,step,t_enter,dur,status,rc\n");
    fflush(out);
}

/* engine name (cached per handle); keep = passes the ADIOS_WAIT_ENGINE filter */
static int engine_info(const adios2_engine *e, const char **name)
{
    int free_slot = -1;
    for (int i = 0; i < NCACHE; i++) {
        if (cache[i].e == e) { *name = cache[i].name; return cache[i].keep; }
        if (!cache[i].e && free_slot < 0) free_slot = i;
    }
    static char tmp[512];
    char *dst = free_slot >= 0 ? cache[free_slot].name : tmp;
    snprintf(dst, 512, "?");
    if (real_name) {
        size_t n = 0;
        if (real_name(NULL, &n, e) == 0 && n > 0) {
            char *buf = malloc(n + 1);
            if (buf && real_name(buf, &n, e) == 0) {
                buf[n] = 0;
                for (char *c = buf; *c; c++) if (*c == ',' || *c == '\n') *c = ';';
                snprintf(dst, 512, "%s", buf);
            }
            free(buf);
        }
    }
    int keep = !*filt || strstr(dst, filt) != NULL;
    if (free_slot >= 0) { cache[free_slot].e = e; cache[free_slot].keep = keep; }
    *name = dst;
    return keep;
}

static void log_call(const adios2_engine *e, const char *call, long step, double t0, double dur, int status, int rc)
{
    const char *name;
    if (!engine_info(e, &name)) return;
    fprintf(out, "%s,%s,%s,%s,%s,%ld,%.6f,%.6f,%d,%d\n", tag, rank, host, name, call, step, t0, dur, status, rc);
    fflush(out);
}

adios2_error adios2_begin_step(adios2_engine *engine, int mode, float timeout_seconds, int *status)
{
    pthread_mutex_lock(&lock); init_once(); pthread_mutex_unlock(&lock);
    double t0 = now_real(), m0 = now_mono();
    adios2_error rc = real_begin(engine, mode, timeout_seconds, status);
    double dur = now_mono() - m0;
    if (!disabled) {
        long step = -1;
        if (rc == 0 && status && *status == 0 && real_cstep) {
            size_t s = 0;
            if (real_cstep(&s, engine) == 0) step = (long)s;
        }
        pthread_mutex_lock(&lock);
        log_call(engine, "begin_step", step, t0, dur, status ? *status : -9, rc);
        pthread_mutex_unlock(&lock);
    }
    return rc;
}

adios2_error adios2_end_step(adios2_engine *engine)
{
    pthread_mutex_lock(&lock); init_once(); pthread_mutex_unlock(&lock);
    double t0 = now_real(), m0 = now_mono();
    adios2_error rc = real_end(engine);
    double dur = now_mono() - m0;
    if (!disabled) {
        pthread_mutex_lock(&lock);
        log_call(engine, "end_step", -1, t0, dur, -9, rc);
        pthread_mutex_unlock(&lock);
    }
    return rc;
}

adios2_error adios2_close(adios2_engine *engine)
{
    pthread_mutex_lock(&lock); init_once(); pthread_mutex_unlock(&lock);
    const char *name = NULL;
    int keep = 1;
    if (!disabled) {               /* look the name up before the handle is closed */
        pthread_mutex_lock(&lock);
        keep = engine_info(engine, &name);
        pthread_mutex_unlock(&lock);
    }
    double t0 = now_real(), m0 = now_mono();
    adios2_error rc = real_close(engine);
    double dur = now_mono() - m0;
    if (!disabled && keep) {
        pthread_mutex_lock(&lock);
        fprintf(out, "%s,%s,%s,%s,close,-1,%.6f,%.6f,-9,%d\n", tag, rank, host, name, t0, dur, rc);
        fflush(out);
        for (int i = 0; i < NCACHE; i++) if (cache[i].e == engine) cache[i].e = NULL;  /* handle may be reused */
        pthread_mutex_unlock(&lock);
    }
    return rc;
}
