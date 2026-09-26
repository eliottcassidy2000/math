/* procgen_landing_20260926_scan.c -- exhaustive worst-case landing multiplicity.
 *
 * Lane `landing`, session collatz-procgen-20260922, 2026-09-26.
 *
 * For T_b(x) = x/2 (x even), (3x+b)/2 (x odd), scale X = 2^(k+1) - 1 (so floor(log2 X) = k)
 * and half-integer depths D = D2/2, every y in [1, X] is taken as a candidate FIRST dipper:
 * with y_0 = y, s = least index in [1, k] with y_s < y 2^(-D) (the landing index), the
 * multiplicity counted from y is the number of t in [0, s) with y_t <= X whose first index
 * below y_t 2^(-D) is s.  Every landing configuration whose dippers are <= X has its first
 * dipper in [1, X], so the maximum over y is the worst case over all segments at scale X.
 * All comparisons are exact (y_s < y 2^(-D2/2)  <=>  y_s^2 2^D2 < y^2, in 128-bit).
 * Checks per configuration: shell (2^D y_s < y_t <= 2^(D+1) y_s), odd separation between
 * consecutive dippers, halving into the landing index, and the Proposition U bound.
 * Output: one line per (b, k, D2): max multiplicity, bound, number of maximisers, one
 * maximiser, heavy counts H_M = #{y : multiplicity from y >= M}, and violation counters.
 *
 * Usage: scan b kmin kmax D2list(comma)      e.g.  scan 1 10 24 2,3,4,5,6,7,8,10,12
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

typedef __int128 i128;
#define KMAX 40
#define MMAX 64

static int below(int64_t u, int64_t v, int D2) {   /* u < v 2^(-D2/2), u, v > 0 */
    i128 lhs = (i128)u * (i128)u;
    i128 rhs = (i128)v * (i128)v;
    return (lhs << D2) < rhs;
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: scan b kmin kmax D2list\n"); return 1; }
    int b = atoi(argv[1]), kmin = atoi(argv[2]), kmax = atoi(argv[3]);
    int D2s[32], nD = 0;
    char *tok = strtok(argv[4], ",");
    while (tok && nD < 32) { D2s[nD++] = atoi(tok); tok = strtok(NULL, ","); }
    const double alpha = log2(3.0);
    int ab = b < 0 ? -b : b;
    for (int k = kmin; k <= kmax; k++) {
        int64_t X = ((int64_t)1 << (k + 1)) - 1;
        int64_t maxm[32], nmax[32], argy[32], viol[32], skipped[32];
        int64_t heavy[32][MMAX + 1];
        memset(maxm, 0, sizeof maxm); memset(nmax, 0, sizeof nmax); memset(argy, 0, sizeof argy);
        memset(viol, 0, sizeof viol); memset(skipped, 0, sizeof skipped); memset(heavy, 0, sizeof heavy);
        int64_t seg[KMAX + 2];
        for (int64_t y = 1; y <= X; y++) {
            seg[0] = y;
            int ok = 1;
            for (int t = 1; t <= k; t++) {
                int64_t v = seg[t - 1];
                seg[t] = (v & 1) ? (3 * v + b) / 2 : v / 2;
                if (seg[t] <= 0) { ok = 0; break; }
            }
            if (!ok) continue;
            for (int d = 0; d < nD; d++) {
                int D2 = D2s[d];
                int s = 0;
                for (int t = 1; t <= k; t++) if (below(seg[t], y, D2)) { s = t; break; }
                if (!s) continue;
                /* distinctness of seg[0..s] */
                int dist = 1;
                for (int p = 0; p <= s && dist; p++)
                    for (int q = p + 1; q <= s; q++) if (seg[p] == seg[q]) { dist = 0; break; }
                if (!dist) { skipped[d]++; continue; }
                int m = 0, last = -1, bad = 0;
                for (int t = 0; t < s; t++) {
                    if (seg[t] > X) continue;
                    if (!below(seg[s], seg[t], D2)) continue;
                    int first = 1;
                    for (int u = t + 1; u < s; u++) if (below(seg[u], seg[t], D2)) { first = 0; break; }
                    if (!first) continue;
                    /* shell: y_s < y_t 2^-D (true) and y_t 2^-D <= 2 y_s */
                    if (below(2 * seg[s], seg[t], D2)) bad++;
                    if (last >= 0) {
                        int odd = 0;
                        for (int u = last; u < t; u++) if (seg[u] & 1) { odd = 1; break; }
                        if (!odd) bad++;
                    }
                    last = t;
                    m++;
                }
                if (m >= 2 && seg[s - 1] != 2 * seg[s]) bad++;
                int fD = D2 / 2;
                int ub = (b > 0) ? (int)ceil((k - fD) / alpha - 1e-12)
                                 : (int)ceil((k - fD - 1) / alpha - 1e-12) + 1;
                if (m > ub && (b > 0 || seg[s] >= (int64_t)k * ab)) bad++;
                viol[d] += bad;
                if (m > MMAX) m = MMAX;
                for (int M = 1; M <= m; M++) heavy[d][M]++;
                if (b < 0 && seg[s] < (int64_t)k * ab) continue;  /* bound applies above k|b| */
                if (m > maxm[d]) { maxm[d] = m; nmax[d] = 1; argy[d] = y; }
                else if (m == maxm[d]) nmax[d]++;
            }
        }
        for (int d = 0; d < nD; d++) {
            int fD = D2s[d] / 2;
            int ub = (b > 0) ? (int)ceil((k - fD) / alpha - 1e-12)
                             : (int)ceil((k - fD - 1) / alpha - 1e-12) + 1;
            printf("b=%d k=%d D2=%d max=%lld ub=%d nmax=%lld argy=%lld viol=%lld skip=%lld H=",
                   b, k, D2s[d], (long long)maxm[d], ub, (long long)nmax[d], (long long)argy[d],
                   (long long)viol[d], (long long)skipped[d]);
            for (int M = 1; M <= maxm[d] && M <= MMAX; M++)
                printf("%s%lld", M > 1 ? "," : "", (long long)heavy[d][M]);
            printf("\n");
        }
        fflush(stdout);
    }
    return 0;
}
