/*
 * chanpak_f1f2_cover35_bucketed_20261106.c  (opus-2026-11-06-S82)
 *
 * Same families and output as chanpak_f1f2_cover34_20261106.c (HYP-9478 P5-P7), reorganised for memory locality:
 *   F2 values go to bitmap A, F1 values to bitmap B (one write per value), union = A | B at the end;
 *   values are staged in thread-local buckets of 2^BSH consecutive integers and flushed per bucket, so the
 *   atomic ORs of a flush hit one cache-sized region of the bitmap.
 * g: table for v < 2^27, exact suffix law beyond (checked on random v and at the endpoints, as in the 2^34 program).
 * Usage: cover35b LOGX bankfile
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
typedef long long ll;
typedef unsigned long long ull;
#define BSH 23                 /* bucket = 2^23 consecutive integers = 2^22 odd numbers = 512 KB of bitmap */
#define CAP 512                /* staged values per bucket per thread */
static uint32_t *gtab; static ll GT;
static ll Ctab[256], Etab[256], Ktab[256];
static ll g_rec(ll v) { if (v == 0) return 1; if (v & 1) return g_rec(v >> 1); ll u = v >> 1; return g_rec(u) + u + (u & 1); }
static inline ll gval(ll v) { if (v < GT) return gtab[v]; ll k = v >> 8, r = v & 255; return gval(k) + k * Ctab[r] + (k & 1) * Etab[r] + Ktab[r]; }
static ll X; static int NB;    /* number of buckets */
typedef struct { uint32_t *buf; int *cnt; } Stage;
static inline void flushb(ull *bm, Stage *s, int b) {
    uint32_t *p = s->buf + (size_t)b * CAP; int c = s->cnt[b]; ll base = (ll)b << BSH;
    for (int i = 0; i < c; i++) { ll N = base + p[i]; ll idx = N >> 1; __atomic_fetch_or(&bm[idx >> 6], 1ULL << (idx & 63), __ATOMIC_RELAXED); }
    s->cnt[b] = 0;
}
static inline void push(ull *bm, Stage *s, ll N) {
    if (N < 1 || N > X || !(N & 1)) return;
    int b = (int)(N >> BSH); uint32_t off = (uint32_t)(N & ((1LL << BSH) - 1));
    s->buf[(size_t)b * CAP + s->cnt[b]] = off;
    if (++s->cnt[b] == CAP) flushb(bm, s, b);
}
static inline int getb(const ull *b, ll N) { ll i = N >> 1; return (int)((b[i >> 6] >> (i & 63)) & 1); }
int main(int argc, char **argv) {
    int LOGX = atoi(argv[1]); X = 1LL << LOGX; NB = (int)((X >> BSH) + 1);
    FILE *f = fopen(argv[2], "r"); if (!f) { fprintf(stderr, "no bank\n"); return 1; }
    static ll yh[4096], ye[4096], yR[4096]; int nb = 0;
    while (nb < 4096 && fscanf(f, "%lld %lld %lld", &yh[nb], &ye[nb], &yR[nb]) == 3) nb++;
    fclose(f);
    printf("# chanpak_f1f2_cover35_bucketed_20261106 LOGX=%d bank=%s (%d triples) threads=%d buckets=%d\n", LOGX, argv[2], nb, omp_get_max_threads(), NB);
    int GTB = LOGX < 27 ? LOGX + 1 : 27; GT = 1LL << GTB;
    gtab = malloc(sizeof(uint32_t) * (size_t)GT); gtab[0] = 1;
    for (ll v = 1; v < GT; v++) { if (v & 1) gtab[v] = gtab[v >> 1]; else { ll u = v >> 1; gtab[v] = (uint32_t)((ll)gtab[u] + u + (u & 1)); } }
    for (int r = 0; r < 256; r++) {
        ll C = 0, E = 0, K = 0, val = 0; int prev = -1;
        for (int j = 1; j <= 8; j++) { int bit = (r >> (8 - j)) & 1; if (!bit) { C += 1LL << (j - 1); if (j == 1) E = 1; else K += val + prev; } val = 2 * val + bit; prev = bit; }
        Ctab[r] = C; Etab[r] = E; Ktab[r] = K;
    }
    {   ull s = 88172645463325252ULL; ll bad = 0;
        for (int t = 0; t < 2000000; t++) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; ll v = (ll)(s % (ull)(X + 1)); if (g_rec(v) != gval(v)) bad++; }
        ll ends[3] = {X - 2, X - 1, X}; for (int t = 0; t < 3; t++) if (g_rec(ends[t]) != gval(ends[t])) bad++;
        printf("g table + suffix law vs recursion: 2,000,000 random v <= X and v = X-2, X-1, X: failures %lld\n", bad);
        if (bad) return 2;
    }
    size_t words = (size_t)((X >> 1) / 64 + 2);
    ull *A = calloc(words, sizeof(ull)), *B = calloc(words, sizeof(ull));
    if (!A || !B) { fprintf(stderr, "alloc failed\n"); return 3; }
    int nth = omp_get_max_threads(); Stage *st = malloc(sizeof(Stage) * nth);
    for (int t = 0; t < nth; t++) { st[t].buf = malloc(sizeof(uint32_t) * (size_t)NB * CAP); st[t].cnt = calloc(NB, sizeof(int)); if (!st[t].buf || !st[t].cnt) { fprintf(stderr, "stage alloc\n"); return 3; } }
    double t0 = omp_get_wtime(); ll evals = 0, pairs = 0;
    /* F2 -> A */
    #pragma omp parallel reduction(+:evals)
    {
        Stage *s = &st[omp_get_thread_num()];
        #pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < nb; i++)
            for (ll v = 0; (1 + v) * yR[i] <= X; v++) { ll ev = (v == 0 || !(v & 1)); push(A, s, gval(v) * yh[i] - ev * ye[i] + (1 + v) * yR[i]); evals++; }
        for (int b = 0; b < NB; b++) if (s->cnt[b]) flushb(A, s, b);
    }
    printf("F2 done %.1fs (%lld evaluations)\n", omp_get_wtime() - t0, evals); fflush(stdout);
    /* F1 -> B */
    ll umax = 0; while ((umax + 2) * (umax + 2) <= X + 1) umax++;
    #pragma omp parallel reduction(+:pairs)
    {
        Stage *s = &st[omp_get_thread_num()];
        #pragma omp for schedule(dynamic, 1)
        for (ll u = 0; u <= umax; u++) {
            ll gu = gval(u), eu = (u == 0 || !(u & 1));
            for (ll v = u; (1 + u) * (1 + v) <= X; v++) { ll ev = (v == 0 || !(v & 1)); push(B, s, (1 + u) * (1 + v) + gu * gval(v) - eu * ev); pairs++; }
        }
        for (int b = 0; b < NB; b++) if (s->cnt[b]) flushb(B, s, b);
    }
    printf("F1 done %.1fs (%lld pairs)\n", omp_get_wtime() - t0, pairs); fflush(stdout);
    ll mu = 0, m2 = 0;
    printf("odd misses of F1 u F2 in [1, 2^%d]:", LOGX);
    for (ll N = 1; N <= X; N += 2) if (!getb(A, N) && !getb(B, N)) { mu++; if (mu <= 100) printf(" %lld", N); }
    printf("  [total %lld]\n", mu);
    printf("odd misses of F2 alone in [1, 2^%d]:", LOGX);
    for (ll N = 1; N <= X; N += 2) if (!getb(A, N)) { m2++; if (m2 <= 100) printf(" %lld", N); }
    printf("  [total %lld]\n", m2);
    printf("F1-only misses per block (2^k, 2^(k+1)): k count ratio density\n");
    ll prevm = -1;
    for (int k = 10; k < LOGX; k++) {
        ll lo = 1LL << k, hi = 1LL << (k + 1), m = 0;
        for (ll N = lo + 1; N < hi; N += 2) if (!getb(B, N)) m++;
        printf("  k=%d %lld %.4f %.6f\n", k, m, prevm > 0 ? (double)m / prevm : 0.0, (double)m / (double)(lo / 2));
        prevm = m;
    }
    printf("elapsed %.1fs\n", omp_get_wtime() - t0);
    return 0;
}
