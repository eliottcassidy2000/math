/* aud_point.c -- audit F1 (independent referee), from-scratch checks of THM-4690 statements 1-4 at chosen N.
   Written independently of fap.c / fap_lean.c / fap_levels.c (only the definitions in THM-4667 / HYP-9290 / THM-4690 are used).

   Definitions used (from the theorem texts):
     seeds = odd s with s % 8 != 5; seed(x) strips B(x) = 4x+1 while x % 8 == 5;
     P(s) = seed((3s+1)/2^a), a = nu_2(3s+1) in {1,2}, step sign (-1)^a, literal step clock a + 2*(#strips);
     sources: every odd m <= 4N carries K_N(m) = w(m) + thr(m)[m < 2N], w = 2 if m%3==2 else -1, thr = 1,-2,1 (m%3=0,1,2),
     aggregated on seed(m) -> c_0(x);
     f(n) = (-1)^(first-hit number of shortcut steps T to reach 1), T(n) = n/2 or (3n+1)/2;
     H(M) = sum_{n<=M} w(n) f(n), R(N) = H(4N) - H(N).
   Everything is recomputed from scratch for each N (no incremental update), with exact integer arithmetic.

   Usage: ./aud_point Vfile N1 [N2 ...]      (prints one block of checks per N; nonzero exit code on any failure) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef __int128 i128;
static int nV; static int64_t *Vs; static int64_t maxV; static uint8_t *inVarr;
static inline int inV(int64_t s) { return s >= 1 && s <= maxV && inVarr[s]; }
static int vindex(int64_t s) { int lo = 0, hi = nV - 1; while (lo <= hi) { int m = (lo + hi) >> 1; if (Vs[m] == s) return m; if (Vs[m] < s) lo = m + 1; else hi = m - 1; } return -1; }
static inline int64_t seed_of(int64_t m) { while ((m & 7) == 5) m = (m - 1) >> 2; return m; }
static inline int w_of(int64_t m) { return (m % 3 == 2) ? 2 : -1; }
static inline int thr_of(int64_t m) { return (m % 3 == 1) ? -2 : 1; }
static inline int is_seed(int64_t s) { return (s & 1) && ((s & 7) != 5); }
static int64_t Pstep(int64_t s, int *a, int *r) {
    if (s > (INT64_MAX - 1) / 3) { fprintf(stderr, "overflow at %lld\n", (long long)s); exit(9); }
    int64_t t = 3 * s + 1; int aa = 0; while (!(t & 1)) { t >>= 1; aa++; }
    int rr = 0; while ((t & 7) == 5) { t = (t - 1) >> 2; rr++; }
    *a = aa; *r = rr; return t;
}
static int fails = 0;
#define CHECK(cond, ...) do { if (!(cond)) { fails++; printf("  !! FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

/* ---- open-addressing map int64 -> int32 node index ---- */
typedef struct { int64_t *key; int32_t *val; uint64_t cap; } Map;
static void map_init(Map *m, uint64_t want) { uint64_t c = 1; while (c < 2 * want) c <<= 1; m->cap = c; m->key = calloc(c, 8); m->val = malloc(4 * c); if (!m->key || !m->val) { fprintf(stderr, "oom map\n"); exit(8); } }
static void map_free(Map *m) { free(m->key); free(m->val); }
static inline uint64_t hsh(int64_t x) { uint64_t z = (uint64_t)x; z ^= z >> 33; z *= 0xff51afd7ed558ccdULL; z ^= z >> 33; z *= 0xc4ceb9fe1a85ec53ULL; z ^= z >> 33; return z; }
static int32_t map_get(Map *m, int64_t k) { uint64_t i = hsh(k) & (m->cap - 1); while (m->key[i]) { if (m->key[i] == k) return m->val[i]; i = (i + 1) & (m->cap - 1); } return -1; }
static int nnodes = 0, capnodes = 0;
static int64_t *nd_val; static int32_t *nd_dep, *nd_par; static int8_t *nd_sg; static int64_t *nd_J; static int32_t *nd_c0;
static int32_t map_put(Map *m, int64_t k) { /* returns index, creating node if needed */
    uint64_t i = hsh(k) & (m->cap - 1);
    while (m->key[i]) { if (m->key[i] == k) return m->val[i]; i = (i + 1) & (m->cap - 1); }
    if ((uint64_t)nnodes * 2 >= m->cap) { fprintf(stderr, "map full\n"); exit(8); }
    if (nnodes == capnodes) { capnodes = capnodes ? 2 * capnodes : 1 << 20;
        nd_val = realloc(nd_val, 8 * (size_t)capnodes); nd_dep = realloc(nd_dep, 4 * (size_t)capnodes); nd_par = realloc(nd_par, 4 * (size_t)capnodes);
        nd_sg = realloc(nd_sg, (size_t)capnodes); nd_J = realloc(nd_J, 8 * (size_t)capnodes); nd_c0 = realloc(nd_c0, 4 * (size_t)capnodes);
        if (!nd_val || !nd_dep || !nd_par || !nd_sg || !nd_J || !nd_c0) { fprintf(stderr, "oom nodes\n"); exit(8); } }
    m->key[i] = k; m->val[i] = nnodes; nd_val[nnodes] = k; nd_dep[nnodes] = -1; nd_par[nnodes] = -1; nd_sg[nnodes] = 0; nd_J[nnodes] = 0; nd_c0[nnodes] = 0;
    return nnodes++;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s Vfile N1 [N2 ...]\n", argv[0]); return 1; }
    /* ---- read V ---- */
    FILE *fv = fopen(argv[1], "r"); if (!fv) { perror(argv[1]); return 1; }
    char line[512]; int cap = 1024; Vs = malloc(8 * cap);
    while (fgets(line, sizeof line, fv)) { if (line[0] == '#') continue; long long v; if (sscanf(line, "%lld", &v) == 1) { if (nV == cap) { cap *= 2; Vs = realloc(Vs, 8 * cap); } Vs[nV++] = v; } }
    fclose(fv);
    for (int i = 1; i < nV; i++) for (int j = i; j > 0 && Vs[j - 1] > Vs[j]; j--) { int64_t t = Vs[j]; Vs[j] = Vs[j - 1]; Vs[j - 1] = t; }
    maxV = Vs[nV - 1]; inVarr = calloc(maxV + 1, 1);
    for (int i = 0; i < nV; i++) { if (!is_seed(Vs[i])) { printf("V element %lld is not a seed\n", (long long)Vs[i]); fails++; } if (inVarr[Vs[i]]) { printf("duplicate %lld\n", (long long)Vs[i]); fails++; } inVarr[Vs[i]] = 1; }
    int closed = 1, nonclosed = 0; for (int i = 0; i < nV; i++) { int a, r; int64_t y = Pstep(Vs[i], &a, &r); if (!inV(y)) { closed = 0; nonclosed++; } }
    printf("V=%s |V|=%d max=%lld forward-closed(P(V) in V): %s (%d elements leave V)\n", argv[1], nV, (long long)maxV, closed ? "yes" : "NO", nonclosed);

    for (int ai = 2; ai < argc; ai++) {
        int64_t N = atoll(argv[ai]), M = 4 * N;
        int fails0 = fails;
        printf("\n==== N = %lld (M = 4N = %lld) ====\n", (long long)N, (long long)M);
        /* ---- shortcut stopping times st(n), n <= SM ---- */
        int64_t SM = M > maxV ? M : maxV;
        int16_t *st = malloc(2 * (SM + 1)); st[0] = 0; st[1] = 0;
        for (int64_t n = 2; n <= SM; n++) { int64_t x = n; int k = 0; while (x >= n) { x = (x & 1) ? (3 * x + 1) / 2 : x / 2; k++; } st[n] = (int16_t)(k + st[x]); }
        /* ---- R(N) by the definition ---- */
        int64_t H4 = 0, H1 = 0;
        for (int64_t n = 1; n <= M; n++) { int64_t c = (int64_t)w_of(n) * ((st[n] & 1) ? -1 : 1); H4 += c; if (n <= N) H1 += c; }
        int64_t Rdef = H4 - H1;
        /* ---- source coefficients c_0 on seeds x <= M (index (x-1)/2) ---- */
        int64_t ns = (M + 1) / 2;
        int32_t *c0 = calloc(ns, 4);
        for (int64_t m = 1; m <= M; m += 2) { int K = w_of(m) + (m < 2 * N ? thr_of(m) : 0); c0[(seed_of(m) - 1) / 2] += K; }
        /* closed form of c_0: chain x, B x, ... ; L = #chain elements < 2N; ends = [B^L(x) <= 4N]; c0 = P(x mod 3, L mod 3) + ends*lambda((x+L) mod 3) */
        static const int kap[3] = {0, -3, 3}, lam[3] = {-1, -1, 2};
        int64_t cf_bad = 0, absmax = 0; int64_t valcnt[7] = {0};
        for (int64_t i = 0; i < ns; i++) { int64_t x = 2 * i + 1; if (!is_seed(x)) { if (c0[i]) { cf_bad++; } continue; }
            int L = 0; int64_t b = x; while (b < 2 * N) { L++; b = 4 * b + 1; } int ends = (b <= M);
            int r = (int)(x % 3), Pv = 0; for (int j = 0; j < L % 3; j++) Pv += kap[(r + j) % 3];
            int cf = Pv + (ends ? lam[(r + L) % 3] : 0);
            if (cf != c0[i]) cf_bad++;
            if (llabs(c0[i]) > absmax) absmax = llabs(c0[i]);
            if (c0[i] >= -3 && c0[i] <= 3) valcnt[c0[i] + 3]++; }
        printf("c_0: max|c_0| = %lld ; value counts -3..3: %lld %lld %lld %lld %lld %lld %lld ; closed-form mismatches %lld\n", (long long)absmax,
               (long long)valcnt[0], (long long)valcnt[1], (long long)valcnt[2], (long long)valcnt[3], (long long)valcnt[4], (long long)valcnt[5], (long long)valcnt[6], (long long)cf_bad);
        CHECK(absmax <= 3, "|c_0| <= 3"); CHECK(cf_bad == 0, "closed form of c_0");
        CHECK(valcnt[1] == 0 && valcnt[4] == 0, "c_0 never takes the values -2, 1");
        /* ---- stopped paths with memo (increasing x; reuse a computed y < x not in V) ---- */
        int32_t *uix = malloc(4 * ns), *tcl = malloc(4 * ns), *dep = malloc(4 * ns); int8_t *sg = malloc(ns); int64_t *port = malloc(8 * ns);
        int maxclk = 0, maxdep = 0;
        for (int64_t i = 0; i < ns; i++) {
            int64_t x = 2 * i + 1; if (!is_seed(x)) { uix[i] = -1; continue; }
            int64_t y = x, last = -1; int tt = 0, dd = 0, ss = 1, done = 0;
            while (!inV(y)) {
                if (y < x) { int64_t j = (y - 1) / 2; uix[i] = uix[j]; tcl[i] = tt + tcl[j]; dep[i] = dd + dep[j]; sg[i] = (int8_t)(ss * sg[j]); port[i] = port[j]; done = 1; break; }
                int a, r; int64_t z = Pstep(y, &a, &r); if (a < 1 || a > 2) { printf("a=%d at seed %lld\n", a, (long long)y); fails++; }
                tt += a + 2 * r; dd++; if (a == 1) ss = -ss; last = y; y = z;
            }
            if (!done) { uix[i] = vindex(y); tcl[i] = tt; dep[i] = dd; sg[i] = (int8_t)ss; port[i] = last; }
            if (tcl[i] > maxclk) maxclk = tcl[i]; if (dep[i] > maxdep) maxdep = dep[i];
        }
        /* ---- per-source checks: sign = (-1)^clock; clock vs stopping-time difference ---- */
        int64_t badsign = 0, clkdiff = 0, clkdiff_root = 0, clkdiff_other = 0, nsrc = 0; int64_t clkdiff_vals[8] = {0};
        for (int64_t i = 0; i < ns; i++) { if (uix[i] < 0 || c0[i] == 0) continue; nsrc++; int64_t x = 2 * i + 1;
            if (sg[i] != ((tcl[i] & 1) ? -1 : 1)) badsign++;
            int64_t u = Vs[uix[i]]; int64_t dst = (int64_t)st[x] - st[u];
            if (dst != tcl[i]) { clkdiff++; if (u == 1) clkdiff_root++; else clkdiff_other++; int64_t df = dst - tcl[i]; if (df >= 0 && df < 8) clkdiff_vals[df]++; } }
        printf("sources (seeds with c_0 != 0): %lld ; max depth %d ; max clock %d\n", (long long)nsrc, maxdep, maxclk);
        printf("sign(x->u) != (-1)^t(x): %lld sources\n", (long long)badsign); CHECK(badsign == 0, "sign = (-1)^clock");
        printf("literal clock t(x) != st(x) - st(u(x)): %lld sources (entry u=1: %lld, other: %lld; (st diff - t) = 2 in %lld cases)\n",
               (long long)clkdiff, (long long)clkdiff_root, (long long)clkdiff_other, (long long)clkdiff_vals[2]);
        CHECK(clkdiff_other == 0, "clock = stopping-time difference away from u = 1");
        /* ---- frozen coefficients, invariant ---- */
        int64_t *cV = calloc(nV, 8), *cP = calloc(nV, 8);
        int T = maxclk + 2; int64_t *h = calloc((size_t)nV * T, 8), *q = calloc((size_t)nV * T, 8);
        for (int64_t i = 0; i < ns; i++) { if (uix[i] < 0 || c0[i] == 0) continue; int u = uix[i];
            cV[u] += sg[i] * (int64_t)c0[i]; cP[u] += c0[i]; h[(size_t)u * T + tcl[i]] += sg[i] * (int64_t)c0[i]; q[(size_t)u * T + tcl[i]] += c0[i]; }
        int64_t Rtr = 0; i128 EV = 0, Eplus = 0;
        for (int u = 0; u < nV; u++) { Rtr += ((st[Vs[u]] & 1) ? -1 : 1) * cV[u]; EV += (i128)cV[u] * cV[u]; Eplus += (i128)cP[u] * cP[u]; }
        printf("invariant: R(N) via transport = %lld ; H(4N)-H(N) = %lld ; R^2/N = %.4f\n", (long long)Rtr, (long long)Rdef, (double)Rdef * Rdef / N);
        CHECK(Rtr == Rdef, "invariant R(N) = sum f(u) c_V(u)");
        /* ---- clock-lag Gram, all lags ---- */
        i128 *G = calloc(T, sizeof(i128));
        for (int u = 0; u < nV; u++) { int64_t *hu = h + (size_t)u * T; for (int t = 0; t < T; t++) if (hu[t]) for (int D = 0; t + D < T; D++) G[D] += (i128)hu[t] * hu[t + D]; }
        i128 sumG = G[0]; for (int D = 1; D < T; D++) sumG += 2 * G[D];
        i128 HL = 0; for (int D = 2; D < T; D++) HL += 2 * G[D];
        printf("E_V = %lld (%.4f N) ; G(0) + 2 sum_{D>=1} G(D) = %lld ; G(0)/N = %.4f ; G(1)/G(0) = %.4f ; HL/N (direct, lags >= 2) = %+.4f\n",
               (long long)EV, (double)EV / N, (long long)sumG, (double)G[0] / N, (double)G[1] / (double)G[0], (double)HL / N);
        CHECK(EV == sumG, "E_V = G(0) + 2 sum G(D)");
        printf("G(D)/N, D=0..12:"); for (int D = 0; D <= 12 && D < T; D++) printf(" %+.3f", (double)G[D] / N); printf("\n");
        /* ---- clock spectrum ---- */
        i128 Spi = 0, S0 = 0;
        for (int u = 0; u < nV; u++) { int64_t a = 0, b = 0; for (int t = 0; t < T; t++) { a += ((t & 1) ? -1 : 1) * q[(size_t)u * T + t]; b += q[(size_t)u * T + t]; } Spi += (i128)a * a; S0 += (i128)b * b; }
        printf("spectrum: S(pi) = %lld (E_V %lld) ; S(0) = %lld (E+ %lld, E+/N %.4f, E_V/E+ %.4f)\n", (long long)Spi, (long long)EV, (long long)S0, (long long)Eplus, (double)Eplus / N, (double)EV / (double)Eplus);
        CHECK(Spi == EV, "S(pi) = E_V"); CHECK(S0 == Eplus, "S(0) = E+");
        { /* mean over Lg equispaced angles (exact for a trig polynomial of degree < Lg/2) and cosine series at test angles */
            int Lg = 4 * T; double mean = 0, maxrel = 0; double *cs = malloc(sizeof(double) * T), *sn = malloc(sizeof(double) * T);
            for (int k = 0; k < Lg + 7; k++) {
                double th = (k < Lg) ? 2 * M_PI * k / Lg : (k - Lg + 0.5) * 0.4321;
                for (int t = 0; t < T; t++) { cs[t] = cos(t * th); sn[t] = sin(t * th); }
                double S = 0; for (int u = 0; u < nV; u++) { double re = 0, im = 0; int64_t *qu = q + (size_t)u * T; for (int t = 0; t < T; t++) if (qu[t]) { re += qu[t] * cs[t]; im += qu[t] * sn[t]; } S += re * re + im * im; }
                if (k < Lg) mean += S / Lg;
                double ser = (double)G[0]; for (int D = 1; D < T; D++) ser += 2.0 * ((D & 1) ? -1.0 : 1.0) * (double)G[D] * cos(D * th);
                double rel = fabs(S - ser) / (double)G[0]; if (rel > maxrel) maxrel = rel;
            }
            printf("spectrum: mean over %d angles = %.6f N vs G(0) = %.6f N ; max |S - cosine series| / G(0) = %.2e\n", Lg, mean / N, (double)G[0] / N, maxrel);
            printf("spectrum S(k pi/12)/N, k=0..12:");
            for (int k = 0; k <= 12; k++) { double th = k * M_PI / 12, S = 0; for (int u = 0; u < nV; u++) { double re = 0, im = 0; int64_t *qu = q + (size_t)u * T; for (int t = 0; t < T; t++) if (qu[t]) { re += qu[t] * cos(t * th); im += qu[t] * sin(t * th); } S += re * re + im * im; } printf(" %.3f", S / N); }
            printf("\n");
            CHECK(fabs(mean - (double)G[0]) <= 1e-6 * (double)G[0], "mean of S = G(0)"); CHECK(maxrel < 1e-9, "cosine series");
            free(cs); free(sn);
        }
        /* ---- source ledger ---- */
        i128 E0 = 0; for (int64_t i = 0; i < ns; i++) E0 += (i128)c0[i] * c0[i];
        i128 SIB = 0, SIBce = 0, SIBnv = 0, sibLag1 = 0; int64_t npairs = 0, nmeetV = 0, nmeetVnz = 0, lagbad = 0, lagbad_nv = 0, entrybad_nv = 0, ce_notnv = 0;
        for (int64_t x = 1; 2 * x + 1 <= M; x += 8) {           /* x = 1 (mod 8) */
            int64_t i = (x - 1) / 2, j = x;                       /* index of 2x+1 is x */
            int64_t pr = (int64_t)c0[i] * c0[j]; npairs++;
            SIB += pr;
            int meet = inV(x) || inV(2 * x + 1); if (meet) { nmeetV++; if (pr) nmeetVnz++; }
            int ce = (uix[i] == uix[j]);
            if (ce) SIBce += pr;
            if (!meet) { SIBnv += pr; if (!ce) entrybad_nv++; if (tcl[j] - tcl[i] != 1 || sg[i] != -sg[j]) lagbad_nv++; }
            if (ce && meet) ce_notnv++;
            if (ce && (tcl[j] - tcl[i] == 1)) sibLag1 += pr;
            if (ce && tcl[j] - tcl[i] != 1) lagbad++;
        }
        { /* THM-4690(4) case analysis for siblings x = 1 mod 8 and 2x+1: (A) x ends at level L -> sibling at level L, not ending;
             (B) x not ending at level L >= 1 -> sibling ends at level L-1. Count exceptions (expected: O(1) per level). */
            int64_t nA = 0, nB = 0, exA = 0, exB = 0, exL[64] = {0};
            for (int64_t x = 1; 2 * x + 1 <= M; x += 8) {
                int64_t xs[2] = {x, 2 * x + 1}; int Ls[2], En[2];
                for (int k = 0; k < 2; k++) { int L = 0; int64_t b = xs[k]; while (b < 2 * N) { L++; b = 4 * b + 1; } Ls[k] = L; En[k] = (b <= M); }
                if (Ls[0] == 0) continue;               /* x >= 2N: then 2x+1 > 4N, excluded */
                if (En[0]) { nA++; if (!(Ls[1] == Ls[0] && !En[1])) { exA++; exL[Ls[0] & 63]++; } }
                else       { nB++; if (!(Ls[1] == Ls[0] - 1 && En[1])) { exB++; exL[Ls[0] & 63]++; } }
            }
            int maxper = 0; for (int k = 0; k < 64; k++) if (exL[k] > maxper) maxper = (int)exL[k];
            printf("sibling case analysis: case A %lld pairs (%lld exceptions), case B %lld pairs (%lld exceptions); max exceptions per level %d\n",
                   (long long)nA, (long long)exA, (long long)nB, (long long)exB, maxper);
        }
        double e0d = (double)(E0 - 0) - 34.0 * N / 7, sbd = (double)SIB - 5.0 * N / 7;
        printf("E_0 = %lld, E_0 - 34N/7 = %+.3f ; SIB = %lld, SIB - 5N/7 = %+.3f ; SIB/N = %.6f\n", (long long)E0, e0d, (long long)SIB, sbd, (double)SIB / N);
        printf("sibling pairs (x = 1 mod 8, 2x+1 <= 4N): %lld ; meeting V: %lld (nonzero product %lld) ; pairs outside V with different entry: %lld, with lag != 1 or equal signs: %lld\n",
               (long long)npairs, (long long)nmeetV, (long long)nmeetVnz, (long long)entrybad_nv, (long long)lagbad_nv);
        printf("SIB_V (common entry) = %lld ; SIB over pairs with neither member in V = %lld ; common-entry pairs meeting V: %lld ; common-entry pairs with lag != 1: %lld\n",
               (long long)SIBce, (long long)SIBnv, (long long)ce_notnv, (long long)lagbad);
        printf("|SIB - SIB_V| = %lld <= 9 * #pairs meeting V = %lld <= 9|V| = %d (claimed 18|V| = %d)\n", (long long)llabs((long long)(SIB - SIBce)), (long long)(9 * nmeetVnz), 9 * nV, 18 * nV);
        CHECK(entrybad_nv == 0 && lagbad_nv == 0, "siblings outside V: common entry, lag 1, opposite signs");
        CHECK(llabs((long long)(SIB - SIBce)) <= 9 * nmeetVnz, "|SIB - SIB_V| <= 9 #split pairs");
        /* nonzero-product pairs meeting V (for the 'at most two pairs per seed' remark) */
        { int64_t perseed_max = 0; for (int k = 0; k < nV; k++) { int64_t v = Vs[k]; int c = 0; if ((v & 7) == 1) c++; if ((v & 15) == 3) c++; if (c > perseed_max) perseed_max = c; }
          printf("max number of sibling pairs containing one seed of V: %lld\n", (long long)perseed_max); }
        /* ledger identity, Gam1 two ways */
        i128 Gam1 = G[1] + SIBce;
        i128 Gam1_direct = G[1] + sibLag1;   /* lag-1 pairs other than sibling-source lag-1 pairs */
        i128 lhs = E0 - 2 * SIBce + (G[0] - E0) + 2 * Gam1 + HL;
        printf("ledger: E_0 - 2 SIB_V + (G0 - E0) + 2 Gam1 + HL = %lld (E_V %lld) ; terms /N: E0 %.5f  -2SIB_V %+.5f  G0-E0 %+.4f  2Gam1 %+.4f  HL %+.4f\n",
               (long long)lhs, (long long)EV, (double)E0 / N, -2.0 * (double)SIBce / N, (double)(G[0] - E0) / N, 2.0 * (double)Gam1 / N, (double)HL / N);
        printf("Gam1 = G1 + SIB_V = %lld ; lag-1 pairs other than sibling-source lag-1 pairs = %lld\n", (long long)Gam1, (long long)Gam1_direct);
        CHECK(lhs == EV, "pair ledger identity");
        /* (G0 - E0) = 2 * (same-clock pairs of distinct sources with common entry): sum over (u,t) of [(sum c0)^2 - sum c0^2]/... check via h */
        { i128 sq = 0; for (int64_t i = 0; i < ns; i++) if (uix[i] >= 0 && c0[i]) sq += (i128)c0[i] * c0[i]; CHECK(sq == E0, "all sources counted"); }
        /* ---- backward ladder (statement 1/2): node map over all source paths ---- */
        if (!getenv("NOLADDER")) {
            Map mp; map_init(&mp, (uint64_t)(3 * nsrc + 4096)); nnodes = 0;
            for (int64_t i = 0; i < ns; i++) { if (uix[i] < 0 || c0[i] == 0) continue; int64_t x = 2 * i + 1;
                int d = dep[i]; int32_t prev = -1; int64_t y = x; int s = 1; int8_t lastsg = 0;
                for (int k = d; k >= 0; k--) {
                    int32_t id = map_put(&mp, y);
                    if (nd_dep[id] < 0) nd_dep[id] = k; else if (nd_dep[id] != k) { fails++; if (fails < 20) printf("  !! dep inconsistency at %lld\n", (long long)y); }
                    if (prev >= 0) { if (nd_par[prev] < 0) { nd_par[prev] = id; nd_sg[prev] = lastsg; } else if (nd_par[prev] != id) { fails++; } }
                    nd_J[id] += (int64_t)s * c0[i];
                    if (k == d) nd_c0[id] = c0[i];
                    if (k == 0) { if (!inV(y)) { fails++; printf("  !! path end not in V\n"); } break; }
                    int a, r; int64_t z = Pstep(y, &a, &r); lastsg = (a == 2) ? 1 : -1; s *= lastsg; prev = id; y = z;
                }
            }
            /* basin form (statement 1): J(y) = f(y) R_y(N), R_y = sum over sources x through y of f(x) c_0(x) */
            { int64_t *Rn = calloc(nnodes, 8);
              for (int64_t i = 0; i < ns; i++) { if (uix[i] < 0 || c0[i] == 0) continue; int64_t x = 2 * i + 1; int fx = (st[x] & 1) ? -1 : 1;
                  int32_t id = map_get(&mp, x); while (id >= 0) { Rn[id] += (int64_t)fx * c0[i]; id = nd_par[id]; } }
              int64_t bad = 0;
              for (int id = 0; id < nnodes; id++) { int64_t y = nd_val[id], z = y; int k = 0; while (z > SM) { z = (z & 1) ? (3 * z + 1) / 2 : z / 2; k++; }
                  int fy = ((k + st[z]) & 1) ? -1 : 1; if (nd_J[id] != fy * Rn[id]) bad++; }
              printf("basin form J(y) = f(y) R_y(N): %lld violations over %d nodes\n", (long long)bad, nnodes);
              CHECK(bad == 0, "basin form"); free(Rn); }
            /* E^(k) by the formula, for every k */
            i128 *Ek = calloc(maxdep + 2, sizeof(i128)), *below = calloc(maxdep + 3, sizeof(i128));
            for (int64_t i = 0; i < ns; i++) if (uix[i] >= 0 && c0[i]) { int d = dep[i]; below[d + 1] += (i128)c0[i] * c0[i]; }
            for (int k = 1; k <= maxdep + 2; k++) below[k] += below[k - 1];
            for (int id = 0; id < nnodes; id++) Ek[nd_dep[id]] += (i128)nd_J[id] * nd_J[id];
            for (int k = 0; k <= maxdep; k++) Ek[k] += below[k];
            /* E^(0) = E_V, E^(maxdep) = E_0, E^(1) = sum_V c0^2 + sum_ports J^2 */
            i128 Vc0 = 0, portJ = 0; int64_t nports = 0;
            for (int id = 0; id < nnodes; id++) { if (nd_dep[id] == 0) Vc0 += (i128)nd_c0[id] * nd_c0[id]; if (nd_dep[id] == 1) { portJ += (i128)nd_J[id] * nd_J[id]; nports++; } }
            printf("ladder: nodes %d ; max dep %d ; E^(0) = %lld (E_V %lld) ; E^(maxdep) = %lld (E_0 %lld) ; E^(1) = %lld , sum_V c0^2 + sum_ports J^2 = %lld (%lld ports with mass)\n",
                   nnodes, maxdep, (long long)Ek[0], (long long)EV, (long long)Ek[maxdep], (long long)E0, (long long)Ek[1], (long long)(Vc0 + portJ), (long long)nports);
            CHECK(Ek[0] == EV, "E^(0) = E_V"); CHECK(Ek[maxdep] == E0, "E^(max) = E_0"); CHECK(Ek[1] == Vc0 + portJ, "E^(1) = port-diagonal energy");
            /* direct stopped simulation at selected k: stop at the first path point with dep <= k (membership in V^(k)) */
            int ksel[12] = {0, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144}; int nk = 0;
            int64_t *acc = calloc(nnodes, 8);
            for (int kk = 0; kk < 12; kk++) { int k = ksel[kk]; if (k > maxdep) continue; nk++;
                memset(acc, 0, 8 * (size_t)nnodes);
                for (int64_t i = 0; i < ns; i++) { if (uix[i] < 0 || c0[i] == 0) continue; int64_t y = 2 * i + 1; int s = 1;
                    int32_t id = map_get(&mp, y);
                    while (nd_dep[id] > k) { s *= nd_sg[id]; id = nd_par[id]; }
                    acc[id] += (int64_t)s * c0[i]; }
                i128 Ed = 0; for (int id = 0; id < nnodes; id++) Ed += (i128)acc[id] * acc[id];
                if (Ed != Ek[k]) { fails++; printf("  !! direct E^(%d) = %lld != formula %lld\n", k, (long long)Ed, (long long)Ek[k]); }
                /* frozen points: dep exactly k, or the source itself */
                int64_t badfrz = 0; for (int id = 0; id < nnodes; id++) if (acc[id] && nd_dep[id] != k && nd_c0[id] == 0) badfrz++;
                if (badfrz) { fails++; printf("  !! %lld frozen points with dep != k that are not sources\n", (long long)badfrz); }
            }
            free(acc);
            /* merge identity at every k: E^(k) - E^(k+1) = 2 sum_{dep y = k} [ sum_{z<z'} s s' J J' + c0(y) sum_z s J ] over children with dep k+1 */
            int64_t *sumA = calloc(nnodes, 8); i128 *sumA2 = calloc(nnodes, sizeof(i128));
            for (int id = 0; id < nnodes; id++) if (nd_par[id] >= 0) { int64_t a = nd_sg[id] * nd_J[id]; sumA[nd_par[id]] += a; sumA2[nd_par[id]] += (i128)a * a;
                    if (nd_dep[nd_par[id]] != nd_dep[id] - 1) { fails++; } }
            i128 *rhs = calloc(maxdep + 2, sizeof(i128));
            for (int id = 0; id < nnodes; id++) { i128 S = sumA[id]; rhs[nd_dep[id]] += (S * S - sumA2[id]) + 2 * (i128)nd_c0[id] * S; }
            int64_t merge_bad = 0; for (int k = 0; k < maxdep; k++) if (Ek[k] - Ek[k + 1] != rhs[k]) merge_bad++;
            /* J recursion J(y) = c0(y) + sum_children s J */
            int64_t jrec_bad = 0; for (int id = 0; id < nnodes; id++) if (nd_J[id] != nd_c0[id] + sumA[id]) jrec_bad++;
            printf("ladder: merge identity mismatches over k = 0..%d: %lld ; J-recursion mismatches: %lld ; direct stopped simulation agreed at %d levels\n", maxdep - 1, (long long)merge_bad, (long long)jrec_bad, nk);
            CHECK(merge_bad == 0, "merge identity"); CHECK(jrec_bad == 0, "J recursion");
            printf("ladder E^(k)/N:"); for (int k = 0; k <= maxdep; k++) if (k <= 25 || k % 10 == 0 || k == maxdep) printf(" [%d]%.4f", k, (double)Ek[k] / N); printf("\n");
            { /* where is the dissipation: level reaching half of E_0 - E_V */
                double tot = (double)(E0 - EV), half = (double)EV + tot / 2; int kh = -1; for (int k = 0; k <= maxdep; k++) if ((double)Ek[k] >= half) { kh = k; break; }
                printf("ladder: E_0 - E_V = %.4f N ; first level with E^(k) >= E_V + (E_0 - E_V)/2: %d ; E^(6)-E^(0) = %.4f N\n", tot / N, kh, maxdep >= 6 ? (double)(Ek[6] - Ek[0]) / N : 0.0); }
            /* X(N) = E_V - E^(1) */
            printf("cross-port X/N = (E_V - E^(1))/N = %+.4f\n", (double)(EV - Ek[1]) / N);
            free(Ek); free(below); free(sumA); free(sumA2); free(rhs); map_free(&mp);
        }
        printf("==== N = %lld: %s\n", (long long)N, fails == fails0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED");
        free(st); free(c0); free(uix); free(tcl); free(dep); free(sg); free(port); free(cV); free(cP); free(h); free(q); free(G);
    }
    printf("\nTOTAL FAILURES: %d\n", fails);
    return fails ? 2 : 0;
}
