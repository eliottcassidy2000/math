/*
 * chanpak_one_apex_join_cover_20261104.c  (opus-2026-11-04-S78)
 *
 * Which odd N <= 2^LOGX are Hamiltonian path counts of one-arc joins built from one-apex tournaments?
 *
 * One-arc join (session29 ADDENDUM sections 2 and 5; THM-4866 (B)): for J = A => B with b -> a reversed,
 *     H(J) = x_h y_h - x_e y_e + x_R y_R,
 * x = (h(A), e_a(A), R_a(A)), y = the same triple for (B^op, b); R_a = sum_{X subset V(A)-a} h(A[X]) s_a(A - X).
 * One-apex tournament T(w): chain 1 -> ... -> n plus x with x -> i iff bit i of w is 1; v = value of w. Its triple at x is
 *     (g(v), e(v), 1 + v),  g(0) = 1, g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2),  e(v) = [v = 0 or v even].
 * Families (HYP-9468):
 *   F0  plain counts: the h of every bank triple, and g(v);
 *   F1  one-apex x one-apex: (1+u)(1+v) + g(u) g(v) - e(u) e(v);
 *   F2  one-apex x bank:     g(v) y_h - e(v) y_e + (1+v) y_R  for every bank triple y = (y_h, y_e, y_R).
 * g(v) for large v uses the exact suffix law: if v = k 2^8 + r (k >= 1, r < 2^8), then
 *     g(v) = g(k) + k C(r) + (k mod 2) eps(r) + K(r),
 * with C, eps, K read off the one-apex recurrence by appending the 8 bits of r. This is checked against the plain recursion on
 * 2,000,000 random v before use.
 * Output: every odd miss of F0 + F1 + F2; the F1-only miss count per dyadic block; totals.
 * Usage: cover LOGX bankfile   (bankfile: lines "h e R", e.g. the 2,732 distinct triples of order <= 7)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
typedef long long ll;
typedef unsigned long long ull;
static uint32_t *gtab; static ll GT;               /* g(v) for v < GT */
static ll Ctab[256], Etab[256], Ktab[256];
static ll g_rec(ll v) { if (v == 0) return 1; if (v & 1) return g_rec(v >> 1); ll u = v >> 1; return g_rec(u) + u + (u & 1); }
static inline ll gval(ll v) {
    if (v < GT) return gtab[v];
    ll k = v >> 8, r = v & 255;
    return gval(k) + k * Ctab[r] + (k & 1) * Etab[r] + Ktab[r];
}
static ull *all, *f1; static ll X;
static inline void setbit(ull *b, ll N) { if (N >= 1 && N <= X && (N & 1)) { ll i = N >> 1; __atomic_fetch_or(&b[i >> 6], 1ULL << (i & 63), __ATOMIC_RELAXED); } }
static inline int getbit(const ull *b, ll N) { ll i = N >> 1; return (int)((b[i >> 6] >> (i & 63)) & 1); }
int main(int argc, char **argv) {
    int LOGX = argc > 1 ? atoi(argv[1]) : 24; X = 1LL << LOGX;
    const char *fn = argc > 2 ? argv[2] : "bank7.txt";
    FILE *f = fopen(fn, "r"); if (!f) { fprintf(stderr, "no bank\n"); return 1; }
    int cap = 4096, nb = 0; ll *yh = malloc(sizeof(ll) * cap), *ye = malloc(sizeof(ll) * cap), *yR = malloc(sizeof(ll) * cap); ll a, b, c;
    while (fscanf(f, "%lld %lld %lld", &a, &b, &c) == 3) { if (nb == cap) { cap *= 2; yh = realloc(yh, sizeof(ll) * cap); ye = realloc(ye, sizeof(ll) * cap); yR = realloc(yR, sizeof(ll) * cap); } yh[nb] = a; ye[nb] = b; yR[nb] = c; nb++; }
    fclose(f);
    printf("# chanpak_one_apex_join_cover_20261104 LOGX=%d bank=%s (%d triples) threads=%d\n", LOGX, fn, nb, omp_get_max_threads());
    /* g table */
    GT = (X < (1LL << 26)) ? X + 1 : (1LL << 26);
    gtab = malloc(sizeof(uint32_t) * (size_t)GT); gtab[0] = 1;
    for (ll v = 1; v < GT; v++) { if (v & 1) gtab[v] = gtab[v >> 1]; else { ll u = v >> 1; gtab[v] = (uint32_t)((ll)gtab[u] + u + (u & 1)); } }
    /* suffix tables: start state (h, A, p) = (0, 1, p0); append the 8 bits of r (most significant first) */
    for (int r = 0; r < 256; r++) {
        ll C = 0, E = 0, K = 0, val = 0; int prev = -1;
        for (int j = 1; j <= 8; j++) {
            int bit = (r >> (8 - j)) & 1;
            if (!bit) { C += 1LL << (j - 1); if (j == 1) E = 1; else K += val + prev; }
            val = 2 * val + bit; prev = bit;
        }
        Ctab[r] = C; Etab[r] = E; Ktab[r] = K;
    }
    /* check the suffix law on random v (including v beyond GT) */
    {
        ull s = 88172645463325252ULL; ll bad = 0;
        for (int t = 0; t < 2000000; t++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            ll v = (ll)(s % (ull)(X + 1));
            ll k = v >> 8, r = v & 255, want = g_rec(v), got = (k >= 1) ? g_rec(k) + k * Ctab[r] + (k & 1) * Etab[r] + Ktab[r] : g_rec(v);
            if (want != got || want != gval(v)) bad++;
        }
        printf("suffix law and table check on 2,000,000 random v <= X: failures %lld\n", bad);
        if (bad) return 2;
    }
    size_t words = (size_t)((X >> 1) / 64 + 2);
    all = calloc(words, sizeof(ull)); f1 = calloc(words, sizeof(ull));
    /* F0 */
    for (int i = 0; i < nb; i++) setbit(all, yh[i]);
    #pragma omp parallel for schedule(static)
    for (ll v = 0; v <= X; v++) { ll gv = gval(v); if (gv <= X) setbit(all, gv); }
    /* F1 */
    ll umax = 0; while ((umax + 2) * (umax + 2) <= X + 1) umax++;
    #pragma omp parallel for schedule(dynamic, 1)
    for (ll u = 0; u <= umax; u++) {
        ll gu = gval(u), eu = (u == 0 || !(u & 1));
        for (ll v = u; (1 + u) * (1 + v) <= X; v++) {
            ll ev = (v == 0 || !(v & 1));
            ll N = (1 + u) * (1 + v) + gu * gval(v) - eu * ev;
            if (N <= X) { setbit(all, N); setbit(f1, N); }
        }
    }
    /* F2 */
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < nb; i++) {
        for (ll v = 0; (1 + v) * yR[i] <= X; v++) {
            ll ev = (v == 0 || !(v & 1));
            ll N = gval(v) * yh[i] - ev * ye[i] + (1 + v) * yR[i];
            if (N <= X) setbit(all, N);
        }
    }
    ll miss = 0; printf("odd misses of F0+F1+F2 in [1, 2^%d]:", LOGX);
    for (ll N = 1; N <= X; N += 2) if (!getbit(all, N)) { miss++; if (miss <= 200) printf(" %lld", N); }
    printf("  [total %lld]\n", miss);
    printf("F1-only misses per block [2^k, 2^(k+1)), k, count, ratio, density:\n");
    ll prevm = -1;
    for (int k = 10; k < LOGX; k++) {
        ll lo = 1LL << k, hi = 1LL << (k + 1), m = 0;
        for (ll N = lo + 1; N < hi; N += 2) if (!getbit(f1, N)) m++;
        printf("  k=%d %lld %.4f %.6f\n", k, m, prevm > 0 ? (double)m / prevm : 0.0, (double)m / (double)(lo / 2));
        prevm = m;
    }
    return 0;
}
