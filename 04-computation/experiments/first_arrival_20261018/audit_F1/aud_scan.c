/* aud_scan.c -- audit F1: independent every-N scan of the first-arrival pair ledger (HYP-9325 octave tables).
   Own implementation (memoised stopped paths, own incremental update); definitions as in THM-4667 / HYP-9290 / THM-4690.

   For every N = 1 .. 2^(K+1)-1:
     - c_0 changes at seed(4N-3), seed(4N-1) (by w) and seed(2N-1) (by thr); all tracked quantities are updated exactly;
     - the invariant R(N) = sum_u f(u) c_V(u) is compared with H(4N) - H(N) computed from the definition of f;
     - E_V, G0 = sum h_t^2, G1 = sum h_t h_(t+1), E0, SIB (all sibling source pairs), SIB_V (pairs with a common entry),
       SIB_out (pairs with neither member in V), E+ (unsigned), E^(1) (port-diagonal energy), are tracked;
       HL = E_V - G0 - 2 G1 (statement 3 identity, checked exactly by aud_point at sample N).
   Prints octave statistics in the layout of ledger_table.out plus the cross-port term X = E_V - E^(1) and the extremes of
   E0 - 34N/7 and SIB - 5N/7 over all N.

   Usage: ./aud_scan Vfile K */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static int nV; static int64_t *Vs; static int64_t maxV; static uint8_t *inVarr;
static inline int inV(int64_t s) { return s >= 1 && s <= maxV && inVarr[s]; }
static int vindex(int64_t s) { int lo = 0, hi = nV - 1; while (lo <= hi) { int m = (lo + hi) >> 1; if (Vs[m] == s) return m; if (Vs[m] < s) lo = m + 1; else hi = m - 1; } return -1; }
static inline int64_t seed_of(int64_t m) { while ((m & 7) == 5) m = (m - 1) >> 2; return m; }
static inline int w_of(int64_t m) { return (m % 3 == 2) ? 2 : -1; }
static inline int thr_of(int64_t m) { return (m % 3 == 1) ? -2 : 1; }
static inline int is_seed(int64_t s) { return (s & 1) && ((s & 7) != 5); }
static int cmp64(const void *a, const void *b) { int64_t x = *(const int64_t *)a, y = *(const int64_t *)b; return (x > y) - (x < y); }

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s Vfile K\n", argv[0]); return 1; }
    FILE *fv = fopen(argv[1], "r"); if (!fv) { perror(argv[1]); return 1; }
    char line[512]; int cap = 1024; Vs = malloc(8 * cap);
    while (fgets(line, sizeof line, fv)) { if (line[0] == '#') continue; long long v; if (sscanf(line, "%lld", &v) == 1) { if (nV == cap) { cap *= 2; Vs = realloc(Vs, 8 * cap); } Vs[nV++] = v; } }
    fclose(fv); qsort(Vs, nV, 8, cmp64); maxV = Vs[nV - 1]; inVarr = calloc(maxV + 1, 1); for (int i = 0; i < nV; i++) inVarr[Vs[i]] = 1;
    int K = atoi(argv[2]); int64_t Nmax = ((int64_t)1 << (K + 1)) - 1, M = 4 * Nmax;
    int64_t SM = M > maxV ? M : maxV;
    /* stopping-time parities */
    uint8_t *fpar = malloc(SM + 1); { int16_t *st = malloc(2 * (SM + 1)); st[0] = st[1] = 0;
        for (int64_t n = 2; n <= SM; n++) { int64_t x = n; int k = 0; while (x >= n) { x = (x & 1) ? (3 * x + 1) / 2 : x / 2; k++; } st[n] = (int16_t)(k + st[x]); }
        for (int64_t n = 0; n <= SM; n++) fpar[n] = st[n] & 1; free(st); }
    /* per-seed stopped-path data, memoised */
    int64_t ns = (M + 1) / 2;
    int16_t *uu = malloc(2 * ns), *tt = malloc(2 * ns); int8_t *sg = malloc(ns); int64_t *pv = malloc(8 * ns);
    if (!uu || !tt || !sg || !pv) { fprintf(stderr, "oom\n"); return 5; }
    int maxclk = 0;
    for (int64_t i = 0; i < ns; i++) {
        int64_t x = 2 * i + 1; if (!is_seed(x)) { uu[i] = -1; continue; }
        int64_t y = x, last = -1; int t = 0, s = 1, done = 0;
        while (!inV(y)) {
            if (y < x) { int64_t j = (y - 1) / 2; uu[i] = uu[j]; tt[i] = (int16_t)(t + tt[j]); sg[i] = (int8_t)(s * sg[j]); pv[i] = pv[j]; done = 1; break; }
            int64_t z = 3 * y + 1; int a = 0; while (!(z & 1)) { z >>= 1; a++; }
            int r = 0; while ((z & 7) == 5) { z = (z - 1) >> 2; r++; }
            t += a + 2 * r; if (a == 1) s = -s; last = y; y = z;
        }
        if (!done) { uu[i] = (int16_t)vindex(y); tt[i] = (int16_t)t; sg[i] = (int8_t)s; pv[i] = last; }
        if (tt[i] > maxclk) maxclk = tt[i];
        if (sg[i] != ((tt[i] & 1) ? -1 : 1)) { fprintf(stderr, "sign != (-1)^clock at %lld\n", (long long)x); return 6; }
    }
    /* port ids */
    int64_t *plist = malloc(8 * ns); int64_t np = 0; for (int64_t i = 0; i < ns; i++) if (uu[i] >= 0 && pv[i] >= 0) plist[np++] = pv[i];
    qsort(plist, np, 8, cmp64); int64_t nu = 0; for (int64_t i = 0; i < np; i++) if (nu == 0 || plist[nu - 1] != plist[i]) plist[nu++] = plist[i];
    int32_t *pid = malloc(4 * ns);
    for (int64_t i = 0; i < ns; i++) { if (uu[i] < 0 || pv[i] < 0) { pid[i] = -1; continue; } int64_t lo = 0, hi = nu - 1, k = -1; while (lo <= hi) { int64_t m = (lo + hi) / 2; if (plist[m] == pv[i]) { k = m; break; } if (plist[m] < pv[i]) lo = m + 1; else hi = m - 1; } pid[i] = (int32_t)k; }
    free(pv);
    fprintf(stderr, "V=%s |V|=%d K=%d: %lld seeds, %lld ports, max clock %d\n", argv[1], nV, K, (long long)ns, (long long)nu, maxclk);
    int T = maxclk + 2;
    int64_t *cV = calloc(nV, 8), *cP = calloc(nV, 8), *Dv = calloc(nV, 8), *Ap = calloc(nu + 1, 8), *h = calloc((size_t)nV * T, 8);
    int8_t *c0 = calloc(ns, 1);
    int64_t EV = 0, G0 = 0, G1 = 0, E0 = 0, SIB = 0, SIBV = 0, SIBout = 0, Ep = 0, Eport = 0, R = 0, H4 = 0, H1 = 0;
    int ok = 1; int64_t firstbad = -1;
    enum { EVm, EVx, G0m, R1m, R1lo, R1hi, Z0m, Z0x, GAm, GAx, HLm, HLx, EPm, DIPm, E0lo, E0hi, SVlo, SVhi, Xm, Xlo, Xhi, GAOm, NS };
    static double S[64][NS]; static int64_t cnt[64];
    for (int k = 0; k < 64; k++) for (int q = 0; q < NS; q++) S[k][q] = (q == EVx || q == R1hi || q == Z0x || q == GAx || q == HLx || q == E0hi || q == SVhi || q == Xhi) ? -1e300 : ((q == R1lo || q == E0lo || q == SVlo || q == Xlo) ? 1e300 : 0);
    double e0dev_max = 0, sibdev_max = 0, e0dev_rel = 0, sibdev_rel = 0; int64_t e0dev_at = 0, sibdev_at = 0;
    double sibVdev_rel_from20 = 0, e0dev_rel_from20 = 0, sibdev_rel_from20 = 0;
    for (int64_t N = 1; N <= Nmax; N++) {
        int64_t ms[3] = {4 * N - 3, 4 * N - 1, 2 * N - 1};
        for (int e = 0; e < 3; e++) {
            int64_t m = ms[e]; int d = (e < 2) ? w_of(m) : thr_of(m);
            int64_t x = seed_of(m), i = (x - 1) / 2; int u = uu[i], t = tt[i]; int64_t sd = (int64_t)sg[i] * d;
            /* partner in the sibling pair */
            int64_t pj = -1; if ((x & 7) == 1) pj = x; else if ((x & 15) == 3) pj = ((x - 1) / 2 - 1) / 2;
            if (pj >= 0 && pj < ns) { int64_t cp = c0[pj]; SIB += d * cp; if (uu[pj] == u) SIBV += d * cp; int64_t xp = 2 * pj + 1; if (!inV(x) && !inV(xp)) SIBout += d * cp; }
            E0 += 2 * (int64_t)c0[i] * d + (int64_t)d * d; c0[i] = (int8_t)(c0[i] + d);
            EV += 2 * cV[u] * sd + sd * sd; cV[u] += sd;
            Ep += 2 * cP[u] * d + (int64_t)d * d; cP[u] += d;
            int64_t *hu = h + (size_t)u * T; G1 += sd * ((t > 0 ? hu[t - 1] : 0) + hu[t + 1]); G0 += 2 * hu[t] * sd + sd * sd; hu[t] += sd;
            if (pid[i] < 0) { Eport += 2 * Dv[u] * sd + sd * sd; Dv[u] += sd; } else { int p = pid[i]; Eport += 2 * Ap[p] * sd + sd * sd; Ap[p] += sd; }
            R += (fpar[Vs[u]] ? -1 : 1) * sd;
        }
        for (int64_t n = 4 * N - 3; n <= 4 * N; n++) H4 += (int64_t)w_of(n) * (fpar[n] ? -1 : 1);
        H1 += (int64_t)w_of(N) * (fpar[N] ? -1 : 1);
        if (R != H4 - H1) { if (ok) firstbad = N; ok = 0; }
        int k = 63 - __builtin_clzll((uint64_t)N);
        double n = (double)N, *s = S[k]; cnt[k]++;
        double ev = EV / n, r1 = G0 ? (double)G1 / G0 : 0, z0 = (G0 - E0) / n, ga = (G1 + SIBV) / n, gao = (G1 + SIBout) / n, hl = (EV - G0 - 2 * G1) / n, e0 = E0 / n, sv = SIBV / n, xx = (EV - Eport) / n;
        s[EVm] += ev; if (ev > s[EVx]) s[EVx] = ev; s[G0m] += G0 / n;
        s[R1m] += r1; if (r1 < s[R1lo]) s[R1lo] = r1; if (r1 > s[R1hi]) s[R1hi] = r1;
        s[Z0m] += z0; if (z0 > s[Z0x]) s[Z0x] = z0; s[GAm] += ga; if (ga > s[GAx]) s[GAx] = ga; s[GAOm] += gao;
        s[HLm] += hl; if (hl > s[HLx]) s[HLx] = hl; s[EPm] += Ep / n; s[DIPm] += Ep ? (double)EV / Ep : 0;
        if (e0 < s[E0lo]) s[E0lo] = e0; if (e0 > s[E0hi]) s[E0hi] = e0; if (sv < s[SVlo]) s[SVlo] = sv; if (sv > s[SVhi]) s[SVhi] = sv;
        s[Xm] += xx; if (xx < s[Xlo]) s[Xlo] = xx; if (xx > s[Xhi]) s[Xhi] = xx;
        double de = fabs(E0 - 34.0 * N / 7), ds = fabs(SIB - 5.0 * N / 7), lg = N >= 2 ? log2(n) : 1;
        if (de > e0dev_max) { e0dev_max = de; e0dev_at = N; } if (ds > sibdev_max) { sibdev_max = ds; sibdev_at = N; }
        if (N >= 4 && de / lg > e0dev_rel) e0dev_rel = de / lg; if (N >= 4 && ds / lg > sibdev_rel) sibdev_rel = ds / lg;
        if (k >= 20) { double dv = fabs(sv - 5.0 / 7); if (dv > sibVdev_rel_from20) sibVdev_rel_from20 = dv; if (de / n > e0dev_rel_from20) e0dev_rel_from20 = de / n; if (ds / n > sibdev_rel_from20) sibdev_rel_from20 = ds / n; }
    }
    printf("V=%s |V|=%d K=%d (every N < 2^%d): invariant R(N) = H(4N)-H(N) at every N: %s", argv[1], nV, K, K + 1, ok ? "OK\n" : "FAILED");
    if (!ok) printf(" (first failure N=%lld)\n", (long long)firstbad);
    printf("max_N |E0 - 34N/7| = %.3f (at N=%lld) ; max_N>=4 |E0 - 34N/7|/log2 N = %.3f\n", e0dev_max, (long long)e0dev_at, e0dev_rel);
    printf("max_N |SIB - 5N/7| = %.3f (at N=%lld) ; max_N>=4 |SIB - 5N/7|/log2 N = %.3f   (SIB = all sibling source pairs)\n", sibdev_max, (long long)sibdev_at, sibdev_rel);
    printf("max over N >= 2^20 of |SIB_V/N - 5/7| = %.3e ; of |E0/N - 34/7| = %.3e ; of |SIB/N - 5/7| = %.3e\n", sibVdev_rel_from20, e0dev_rel_from20, sibdev_rel_from20);
    printf(" oct | E_V/N mean max | G0/N mean | G1/G0 mean [lo, hi] | (G0-E0)/N mean max | Gam1/N mean max | HL/N mean max | E+/N mean | E_V/E+ mean | E0/N [lo, hi] | SIB_V/N [lo, hi] | X/N mean [lo, hi] max|X| | Gam1out/N mean\n");
    for (int k = 8; k <= K; k++) { double *s = S[k], c = (double)cnt[k];
        printf(" %2d | %.3f %.3f | %.3f | %+.4f [%+.4f, %+.4f] | %+.3f %+.3f | %+.3f %+.3f | %+.3f %+.3f | %.3f | %.4f | [%.5f, %.5f] | [%.5f, %.5f] | %+.3f [%+.3f, %+.3f] %.3f | %+.3f\n", k,
               s[EVm] / c, s[EVx], s[G0m] / c, s[R1m] / c, s[R1lo], s[R1hi], s[Z0m] / c, s[Z0x], s[GAm] / c, s[GAx], s[HLm] / c, s[HLx], s[EPm] / c, s[DIPm] / c,
               s[E0lo], s[E0hi], s[SVlo], s[SVhi], s[Xm] / c, s[Xlo], s[Xhi], fmax(fabs(s[Xlo]), fabs(s[Xhi])), s[GAOm] / c); }
    return ok ? 0 : 2;
}
