/*
 * aud_labelled_burnside.c  (independent audit of THM-4869 / HYP-9461, opus-2026-11-04-S75 audit)
 *
 * Brute force over ALL labelled tournaments of order n (5..8), with class counts by Burnside's lemma:
 *   #classes with property P = (1/n!) * sum over sigma in S_n of #(labelled T fixed by sigma with P).
 * A tournament can be fixed only by odd-order permutations (all cycles odd), and then each orbit of sigma on
 * unordered pairs carries one free orientation bit, so Fix(sigma) has 2^(#pair orbits) members, enumerated directly.
 * Grouping by cycle type lambda: #classes = sum_lambda |Fix(sigma_lambda) & P| / z_lambda.
 * No canonical form, no vertex-extension generator, no use of Moon's theorem.
 *
 * Per tournament: hc by a subset DP from vertex 0, strong (two BFS), sole exit (out-degree 1), sole entry,
 * 2-strong (T - v strong for all v), and lambda >= 2 (no arc whose deletion leaves T non-strong).
 * Outputs, per n:
 *   t_k, u_k for k = 0..30 (u_k: no vertex of out-degree 1), and labelled counts;
 *   counts of 2-strong classes with hc = k, and of 2-arc-strong classes with a cut vertex with hc = k (k <= 30).
 * Usage: aud_labelled_burnside NMIN NMAX
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define KMAX 30
typedef unsigned long long u64;

static int strong_mask(int n, const uint32_t *out, uint32_t V) {
    if (!V) return 1;
    int s = __builtin_ctz(V);
    uint32_t seen = 1u << s, fr = seen;
    while (fr) { int x = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = out[x] & V & ~seen; seen |= nb; fr |= nb; }
    if (seen != V) return 0;
    /* backward reachability: y reaches x iff out[y] has x */
    seen = 1u << s; fr = seen;
    while (fr) {
        int x = __builtin_ctz(fr); fr &= fr - 1;
        for (uint32_t r = V & ~seen; r; r &= r - 1) { int y = __builtin_ctz(r); if ((out[y] >> x) & 1u) { seen |= 1u << y; fr |= 1u << y; } }
    }
    return seen == V;
}

/* number of Hamiltonian cycles, n <= 8: paths from 0 over masks of {1..n-1} */
static u64 hc_count(int n, const uint32_t *out) {
    uint32_t dp[128][8];
    int m = n - 1; uint32_t full = (1u << m) - 1;
    memset(dp, 0, sizeof(uint32_t) * 8 * (full + 1));
    for (int v = 1; v < n; v++) if ((out[0] >> v) & 1u) dp[1u << (v - 1)][v] = 1;
    for (uint32_t S = 1; S <= full; S++) {
        for (uint32_t r = S; r; r &= r - 1) {
            int b = __builtin_ctz(r), v = b + 1; uint32_t c = dp[S][v]; if (!c) continue;
            uint32_t nb = (out[v] >> 1) & full & ~S;
            while (nb) { int b2 = __builtin_ctz(nb); nb &= nb - 1; dp[S | (1u << b2)][b2 + 1] += c; }
        }
    }
    u64 h = 0;
    for (int v = 1; v < n; v++) if (out[v] & 1u) h += dp[full][v];
    return h;
}

static int lambda_ge2(int n, uint32_t *out) {
    uint32_t full = (1u << n) - 1;
    for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if ((out[a] >> b) & 1u) {
        out[a] &= ~(1u << b);
        int st = strong_mask(n, out, full);
        out[a] |= 1u << b;
        if (!st) return 0;
    }
    return 1;
}

typedef struct { u64 t[KMAX + 1], u[KMAX + 1], s2[KMAX + 1], cv[KMAX + 1], big, labstrong; } Acc;

int main(int argc, char **argv) {
    int NMIN = argc > 1 ? atoi(argv[1]) : 5, NMAX = argc > 2 ? atoi(argv[2]) : 8;
    printf("# aud_labelled_burnside NMIN=%d NMAX=%d threads=%d\n", NMIN, NMAX, omp_get_max_threads());
    for (int n = NMIN; n <= NMAX; n++) {
        u64 nfact = 1; for (int i = 2; i <= n; i++) nfact *= i;
        /* odd partitions of n */
        int parts[64][16], np[64], npart = 0;
        {
            int a[16];
            /* simple recursive generation via explicit stack: parts in non-increasing order, odd */
            int k = 0; (void)k;
            /* enumerate by brute force over multiplicities of 1,3,5,7,9,... */
            for (int m7 = 0; 7 * m7 <= n; m7++) for (int m5 = 0; 7 * m7 + 5 * m5 <= n; m5++)
            for (int m3 = 0; 7 * m7 + 5 * m5 + 3 * m3 <= n; m3++) {
                int m1 = n - 7 * m7 - 5 * m5 - 3 * m3; int c = 0;
                for (int i = 0; i < m7; i++) a[c++] = 7;
                for (int i = 0; i < m5; i++) a[c++] = 5;
                for (int i = 0; i < m3; i++) a[c++] = 3;
                for (int i = 0; i < m1; i++) a[c++] = 1;
                memcpy(parts[npart], a, sizeof(int) * c); np[npart] = c; npart++;
            }
        }
        /* totals weighted by class size n!/z_lambda */
        u64 Wt[KMAX + 1] = {0}, Wu[KMAX + 1] = {0}, Ws2[KMAX + 1] = {0}, Wcv[KMAX + 1] = {0};
        u64 Lt[KMAX + 1] = {0}, Lu[KMAX + 1] = {0};
        u64 Wall = 0, Wstrong = 0;
        for (int pi = 0; pi < npart; pi++) {
            /* build sigma */
            int sigma[16], pos = 0;
            for (int c = 0; c < np[pi]; c++) { int L = parts[pi][c]; for (int i = 0; i < L; i++) sigma[pos + i] = pos + (i + 1) % L; pos += L; }
            /* z_lambda */
            u64 z = 1; { int cnt[16] = {0}; for (int c = 0; c < np[pi]; c++) cnt[parts[pi][c]]++;
                for (int L = 1; L <= n; L++) { for (int i = 0; i < cnt[L]; i++) z *= L; for (int i = 2; i <= cnt[L]; i++) z *= i; } }
            u64 csize = nfact / z;
            /* pair orbits: representative ordered pairs */
            int orbA[64][64], orbB[64][64], orbLen[64], norb = 0; int seen[16][16]; memset(seen, 0, sizeof seen);
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (!seen[i][j]) {
                int a = i, b = j, L = 0;
                for (;;) {
                    int lo = a < b ? a : b, hi = a < b ? b : a;
                    if (seen[lo][hi]) {
                        /* the walk must close on the ordered start pair, never on its reverse */
                        if (!(a == i && b == j)) { fprintf(stderr, "orbit contains a reversed pair\n"); return 1; }
                        break;
                    }
                    seen[lo][hi] = 1; orbA[norb][L] = a; orbB[norb][L] = b; L++;
                    a = sigma[a]; b = sigma[b];
                }
                orbLen[norb] = L; norb++;
            }
            u64 NT = 1ull << norb;
            Acc tot; memset(&tot, 0, sizeof tot);
            #pragma omp parallel
            {
                Acc loc; memset(&loc, 0, sizeof loc);
                #pragma omp for schedule(dynamic, 4096)
                for (u64 code = 0; code < NT; code++) {
                    uint32_t out[16] = {0};
                    for (int o = 0; o < norb; o++) {
                        int bit = (code >> o) & 1;
                        for (int x = 0; x < orbLen[o]; x++) {
                            int a = orbA[o][x], b = orbB[o][x];
                            if (bit) out[a] |= 1u << b; else out[b] |= 1u << a;
                        }
                    }
                    uint32_t full = (1u << n) - 1;
                    if (!strong_mask(n, out, full)) { loc.t[0]++; loc.u[0] += 0; continue; }
                    loc.labstrong++;
                    u64 h = hc_count(n, out);
                    int se = 0; for (int v = 0; v < n; v++) if (__builtin_popcount(out[v]) == 1) se = 1;
                    if (h > KMAX) { loc.big++; continue; }
                    loc.t[h]++; if (!se) loc.u[h]++;
                    int s2 = 1; for (int v = 0; v < n && s2; v++) if (!strong_mask(n, out, full & ~(1u << v))) s2 = 0;
                    if (s2) loc.s2[h]++;
                    else if (lambda_ge2(n, out)) loc.cv[h]++;
                }
                #pragma omp critical
                {
                    for (int k = 0; k <= KMAX; k++) { tot.t[k] += loc.t[k]; tot.u[k] += loc.u[k]; tot.s2[k] += loc.s2[k]; tot.cv[k] += loc.cv[k]; }
                    tot.big += loc.big; tot.labstrong += loc.labstrong;
                }
            }
            for (int k = 0; k <= KMAX; k++) { Wt[k] += csize * tot.t[k]; Wu[k] += csize * tot.u[k]; Ws2[k] += csize * tot.s2[k]; Wcv[k] += csize * tot.cv[k]; }
            Wall += csize * NT; Wstrong += csize * tot.labstrong;
            if (pi == 0 || 1) {
                printf("  n=%d type", n); for (int c = 0; c < np[pi]; c++) printf(" %d", parts[pi][c]);
                printf("  |class|=%llu  pair-orbits=%d  fixed=%llu  fixed strong=%llu\n", csize, norb, NT, tot.labstrong);
            }
            if (np[pi] == n) { for (int k = 0; k <= KMAX; k++) { Lt[k] = tot.t[k]; Lu[k] = tot.u[k]; } }
        }
        int bad = 0;
        #define DIVCHK(x) ((x) % nfact ? (bad = 1, 0ull) : (x) / nfact)
        printf("n=%d  classes=%llu  strong classes=%llu  (Burnside remainders must vanish)\n", n, DIVCHK(Wall), DIVCHK(Wstrong));
        printf("n=%d  k : t_k (labelled) | u_k (labelled) | 2-strong classes | 2-arc-strong cut-vertex classes\n", n);
        for (int k = 1; k <= KMAX; k++) {
            u64 t = DIVCHK(Wt[k]), u = DIVCHK(Wu[k]), s2 = DIVCHK(Ws2[k]), cv = DIVCHK(Wcv[k]);
            if (t || s2 || cv) printf("   k=%2d  t=%llu (%llu)  u=%llu (%llu)  2-strong=%llu  cutvertex-2arcstrong=%llu\n", k, t, Lt[k], u, Lu[k], s2, cv);
        }
        printf("n=%d  Burnside divisibility %s\n", n, bad ? "FAILED" : "ok");
        fflush(stdout);
    }
    return 0;
}
