/*
 * chanpak_f1f2_cover34_20261106.c  (opus-2026-11-06-S82)
 *
 * Odd N <= 2^LOGX attained by the one-arc join families (THM-4875 (D), audited; recertified by Session 203):
 *   F1 = { (1+u)(1+v) + g(u) g(v) - e(u) e(v) : 0 <= u <= v }              (one-apex x one-apex)
 *   F2 = { g(v) y_h - e(v) y_e + (1+v) y_R : v >= 0, y in the order-<=7 bank } (one-apex x rooted order <= 7)
 * g(0) = 1, g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2); e(v) = [v = 0 or v even].
 * g is tabulated for v < 2^GTB and evaluated beyond by the exact suffix law (Session191 one_apex/PROOF.md section 2):
 *   v = k 2^8 + r  =>  g(v) = g(k) + k C(r) + (k mod 2) eps(r) + K(r);
 * the table and the suffix law are checked against the plain recursion on random v, including the endpoints 2^LOGX - 2,
 * 2^LOGX - 1, 2^LOGX (Session 203's endpoint lesson: no out-of-range table access anywhere).
 * Output: odd misses of F1 u F2; odd misses of F2 alone; F1-only miss counts per dyadic block.
 * Usage: cover34 LOGX bankfile
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
typedef long long ll;
typedef unsigned long long ull;
static uint32_t *gtab; static ll GT;
static ll Ctab[256], Etab[256], Ktab[256];
static ll g_rec(ll v) { if (v == 0) return 1; if (v & 1) return g_rec(v >> 1); ll u = v >> 1; return g_rec(u) + u + (u & 1); }
static inline ll gval(ll v) {                 /* v >= 0; recursion depth <= (log2 v - log2 GT)/8 + 1 */
    if (v < GT) return gtab[v];
    ll k = v >> 8, r = v & 255;
    return gval(k) + k * Ctab[r] + (k & 1) * Etab[r] + Ktab[r];
}
static ull *uni, *f2b, *f1b; static ll X;
static inline void setb(ull *b, ll N) { if (N >= 1 && N <= X && (N & 1)) { ll i = N >> 1; __atomic_fetch_or(&b[i >> 6], 1ULL << (i & 63), __ATOMIC_RELAXED); } }
static inline int getb(const ull *b, ll N) { ll i = N >> 1; return (int)((b[i >> 6] >> (i & 63)) & 1); }
int main(int argc, char **argv) {
    int LOGX = atoi(argv[1]); X = 1LL << LOGX;
    FILE *f = fopen(argv[2], "r"); if (!f) { fprintf(stderr, "no bank\n"); return 1; }
    static ll yh[4096], ye[4096], yR[4096]; int nb = 0;
    while (nb < 4096 && fscanf(f, "%lld %lld %lld", &yh[nb], &ye[nb], &yR[nb]) == 3) nb++;
    fclose(f);
    printf("# chanpak_f1f2_cover34_20261106 LOGX=%d bank=%s (%d triples) threads=%d\n", LOGX, argv[2], nb, omp_get_max_threads());
    int GTB = LOGX < 27 ? LOGX + 1 : 27; GT = 1LL << GTB;
    gtab = malloc(sizeof(uint32_t) * (size_t)GT); gtab[0] = 1;
    for (ll v = 1; v < GT; v++) { if (v & 1) gtab[v] = gtab[v >> 1]; else { ll u = v >> 1; gtab[v] = (uint32_t)((ll)gtab[u] + u + (u & 1)); } }
    for (int r = 0; r < 256; r++) {           /* append the 8 bits of r, most significant first */
        ll C = 0, E = 0, K = 0, val = 0; int prev = -1;
        for (int j = 1; j <= 8; j++) {
            int bit = (r >> (8 - j)) & 1;
            if (!bit) { C += 1LL << (j - 1); if (j == 1) E = 1; else K += val + prev; }
            val = 2 * val + bit; prev = bit;
        }
        Ctab[r] = C; Etab[r] = E; Ktab[r] = K;
    }
    {   /* checks: random v <= X and the three endpoints */
        ull s = 88172645463325252ULL; ll bad = 0;
        for (int t = 0; t < 2000000; t++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            ll v = (ll)(s % (ull)(X + 1));
            if (g_rec(v) != gval(v)) bad++;
        }
        ll ends[3] = {X - 2, X - 1, X};
        for (int t = 0; t < 3; t++) if (g_rec(ends[t]) != gval(ends[t])) bad++;
        printf("g table + suffix law vs recursion: 2,000,000 random v <= X and v = X-2, X-1, X: failures %lld\n", bad);
        printf("endpoint values g(X-2) = %lld, g(X-1) = %lld, g(X) = %lld\n", gval(X - 2), gval(X - 1), gval(X));
        if (bad) return 2;
    }
    size_t words = (size_t)((X >> 1) / 64 + 2);
    uni = calloc(words, sizeof(ull)); f2b = calloc(words, sizeof(ull)); f1b = calloc(words, sizeof(ull));
    if (!uni || !f2b || !f1b) { fprintf(stderr, "alloc failed\n"); return 3; }
    double t0 = omp_get_wtime();
    ll umax = 0; while ((umax + 2) * (umax + 2) <= X + 1) umax++;
    ll pairs = 0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:pairs)
    for (ll u = 0; u <= umax; u++) {
        ll gu = gval(u), eu = (u == 0 || !(u & 1));
        for (ll v = u; (1 + u) * (1 + v) <= X; v++) {
            ll ev = (v == 0 || !(v & 1));
            ll N = (1 + u) * (1 + v) + gu * gval(v) - eu * ev; pairs++;
            if (N <= X) { setb(uni, N); setb(f1b, N); }
        }
    }
    printf("F1 done %.1fs (%lld pairs)\n", omp_get_wtime() - t0, pairs); fflush(stdout);
    ll evals = 0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:evals)
    for (int i = 0; i < nb; i++)
        for (ll v = 0; (1 + v) * yR[i] <= X; v++) {
            ll ev = (v == 0 || !(v & 1));
            ll N = gval(v) * yh[i] - ev * ye[i] + (1 + v) * yR[i]; evals++;
            if (N <= X) { setb(uni, N); setb(f2b, N); }
        }
    printf("F2 done %.1fs (%lld evaluations)\n", omp_get_wtime() - t0, evals); fflush(stdout);
    ll mu = 0, m2 = 0;
    printf("odd misses of F1 u F2 in [1, 2^%d]:", LOGX);
    for (ll N = 1; N <= X; N += 2) if (!getb(uni, N)) { mu++; if (mu <= 100) printf(" %lld", N); }
    printf("  [total %lld]\n", mu);
    printf("odd misses of F2 alone in [1, 2^%d]:", LOGX);
    for (ll N = 1; N <= X; N += 2) if (!getb(f2b, N)) { m2++; if (m2 <= 100) printf(" %lld", N); }
    printf("  [total %lld]\n", m2);
    printf("F1-only misses per block (2^k, 2^(k+1)): k count ratio density\n");
    ll prevm = -1;
    for (int k = 10; k < LOGX; k++) {
        ll lo = 1LL << k, hi = 1LL << (k + 1), m = 0;
        for (ll N = lo + 1; N < hi; N += 2) if (!getb(f1b, N)) m++;
        printf("  k=%d %lld %.4f %.6f\n", k, m, prevm > 0 ? (double)m / prevm : 0.0, (double)m / (double)(lo / 2));
        prevm = m;
    }
    printf("elapsed %.1fs\n", omp_get_wtime() - t0);
    return 0;
}
