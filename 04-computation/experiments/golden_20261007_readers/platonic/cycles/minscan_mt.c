/* minscan_mt.c -- multithreaded version of minscan.c (same logic; see that file).
 * usage: minscan_mt XB YB LMAX NTHREADS */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include <stdatomic.h>
typedef unsigned __int128 u128;
static uint64_t XB, YB; static int LMAX;
static _Atomic uint64_t nextchunk = 0; static uint64_t nchunksP, nchunksN;
#define CH 4000000ULL
static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
static uint64_t capped = 0, ovf = 0, maxstopP = 0, maxstopN = 0;
static void report(int sign, uint64_t c, int L) {
    char w[300]; int k = 0; u128 v = c; int t;
    for (t = 0; t < L && t < 299; t++) {
        int b = (int)(v & 1); w[t] = b ? '1' : '0'; k += b;
        v = (sign > 0) ? (b ? (3 * v + 1) >> 1 : v >> 1) : (b ? (3 * v - 1) >> 1 : v >> 1);
    }
    w[t] = 0;
    pthread_mutex_lock(&mu);
    printf("CYCLE sign=%c min|x|=%llu period L=%d odd k=%d word=%s\n", sign > 0 ? '+' : '-',
           (unsigned long long)c, L, k, w); fflush(stdout);
    pthread_mutex_unlock(&mu);
}
static void *worker(void *arg) {
    uint64_t lcap = 0, lovf = 0, lmP = 0, lmN = 0;
    for (;;) {
        uint64_t ch = atomic_fetch_add(&nextchunk, 1);
        if (ch >= nchunksP + nchunksN) break;
        int sign = ch < nchunksP ? +1 : -1;
        uint64_t base = sign > 0 ? ch * CH : (ch - nchunksP) * CH;
        uint64_t lim = sign > 0 ? XB : YB;
        uint64_t lo = base + 1, hi = base + CH; if (hi > lim) hi = lim;
        for (uint64_t c = lo; c <= hi; c++) {
            u128 v = c; int t;
            if (sign > 0) {
                for (t = 1; t <= LMAX; t++) {
                    v = (v & 1) ? (3 * v + 1) >> 1 : v >> 1;
                    if (v == c) { report(+1, c, t); break; }
                    if (v < c) break;
                    if (v >> 120) { lovf++; break; }
                }
                if ((uint64_t)t > lmP && t <= LMAX) lmP = t;
            } else {
                for (t = 1; t <= LMAX; t++) {
                    v = (v & 1) ? (3 * v - 1) >> 1 : v >> 1;
                    if (v == c) { report(-1, c, t); break; }
                    if (v < c) break;
                    if (v >> 120) { lovf++; break; }
                }
                if ((uint64_t)t > lmN && t <= LMAX) lmN = t;
            }
            if (t > LMAX) lcap++;
        }
    }
    pthread_mutex_lock(&mu);
    capped += lcap; ovf += lovf; if (lmP > maxstopP) maxstopP = lmP; if (lmN > maxstopN) maxstopN = lmN;
    pthread_mutex_unlock(&mu);
    return 0;
}
int main(int argc, char **argv) {
    XB = strtoull(argv[1], 0, 10); YB = strtoull(argv[2], 0, 10); LMAX = atoi(argv[3]);
    int nt = atoi(argv[4]);
    nchunksP = (XB + CH - 1) / CH; nchunksN = (YB + CH - 1) / CH;
    pthread_t th[64];
    for (int i = 0; i < nt; i++) pthread_create(&th[i], 0, worker, 0);
    for (int i = 0; i < nt; i++) pthread_join(th[i], 0);
    printf("done XB=%llu YB=%llu LMAX=%d capped=%llu overflow=%llu max_stop_pos=%llu max_stop_neg=%llu\n",
           (unsigned long long)XB, (unsigned long long)YB, LMAX, (unsigned long long)capped,
           (unsigned long long)ovf, (unsigned long long)maxstopP, (unsigned long long)maxstopN);
    return 0;
}
