/* Audit E (2026-10-14), THM-4655 1(e)-(g) on the strong transitive-plus-one family, m = 3..MMAX (default 14).
   Enumerates the valid sets A (path <-> A bijection; checked against DFS-enumerated Hamiltonian paths for every word with
   m <= 13 in tpo_dfs_bruteforce.c) and tallies exact pair counts.  Checks:
   1(e)  w_2 = 1: 1/2 <= P(z<u_1) < 2/3;  w_2 = 0: 1/3 <= P(u_2<u_1) <= 1/2 with 1/3 iff w = 1 0 1^(m-3) 0;
         at the other end (u_(m-1),u_m) or (u_m,z) is balanced; #balanced arcs on the cycle z u_1 .. u_m z >= 2;
         minimum of that number per m and the words attaining it (is 1^a 0^b among them?);
   1(f)  1/3 <= P(u_i < z) <= 2/3 for 2 <= i <= m-1;  and the sufficient condition 2X, 2Y <= W_M for 3 <= i <= (m+1)/2
         (X: u_i forced into A, Y: forced into B, W_M: u_i free), classified from each A directly;
   1(g)  every vertex has a balanced partner; the balanced graph is connected. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int bal(u64 c, u64 H) { return 3 * c >= H && 3 * c <= 2 * H; }
int main(int argc, char **argv) {
    int MMAX = argc > 1 ? atoi(argv[1]) : 14;
    long long tot_fail_e = 0, tot_fail_f = 0, tot_fail_g = 0, tot_fail_xy = 0, tot_xy_cases = 0;
    for (int m = 3; m <= MMAX; m++) {
        int n = m + 1, Z = m;
        long long fe = 0, ff = 0, fg = 0, fxy = 0, nxy = 0, words = 0;
        int mincyc = 99; long long mincyc_words = 0; int ab_attains = 0; char exw[64] = "";
        for (unsigned mid = 0; mid < (1u << (m - 2)); mid++) {
            unsigned w = 1u | (mid << 1);              /* bit i-1 = w_i ; w_1 = 1, w_m = 0 */
            words++;
            static u64 bef[24][24]; memset(bef, 0, sizeof bef);
            u64 H = 0;
            static u64 X[24], Y[24], WM[24]; memset(X, 0, sizeof X); memset(Y, 0, sizeof Y); memset(WM, 0, sizeof WM);
            unsigned full = (1u << m) - 1;
            for (unsigned A = 0; A <= full; A++) {
                unsigned B = full ^ A;
                if (A) { int mx = 31 - __builtin_clz(A); if (w >> mx & 1) continue; }
                if (B) { int mn = __builtin_ctz(B); if (!(w >> mn & 1)) continue; }
                H++;
                /* order: A increasing, z, B increasing.  bef[x][y] for all pairs */
                for (int i = 0; i < m; i++) {
                    if (A >> i & 1) { bef[i][Z]++; for (int j = 0; j < m; j++) if (j != i && ((B >> j & 1) || (j > i))) bef[i][j]++; }
                    else { bef[Z][i]++; for (int j = i + 1; j < m; j++) if (B >> j & 1) bef[i][j]++; }
                }
                /* classification for 1(f): initial segment iff A = {0..k-1} */
                int initial = ((A & (A + 1)) == 0);
                int beta = B ? __builtin_ctz(B) : -1, alpha = A ? 31 - __builtin_clz(A) : -1;   /* 0-based */
                for (int i = 1; i < m - 1; i++) {        /* u_(i+1), 0-based index i, 2 <= i+1 <= m-1 */
                    if (initial) { if (A >> i & 1) X[i]++; else Y[i]++; }
                    else {
                        if (i == alpha || i < beta) X[i]++;
                        else if (i > alpha || i == beta) Y[i]++;
                        else WM[i]++;                    /* beta < i < alpha: free (counts both choices) */
                    }
                }
            }
            /* 1(e) */
            int w2 = w >> 1 & 1;
            unsigned eqw = 1u; for (int i = 3; i <= m - 1; i++) eqw |= 1u << (i - 1);
            if (w2) { if (!(2 * bef[Z][0] >= H && 3 * bef[Z][0] < 2 * H)) fe++; }
            else {
                u64 c = bef[1][0];   /* u_2 before u_1 */
                if (!(3 * c >= H && 2 * c <= H)) fe++;
                if ((3 * c == H) != (w == eqw && m >= 4)) fe++;
                if (c + bef[1][Z] != H) fe++;            /* P(u_2<u_1) = 1 - P(u_2<z) */
            }
            if (!(bal(bef[m - 2][m - 1], H) || bal(bef[m - 1][Z], H))) fe++;
            int ncyc = bal(bef[Z][0], H) + bal(bef[m - 1][Z], H);
            for (int i = 0; i + 1 < m; i++) ncyc += bal(bef[i][i + 1], H);
            if (ncyc < 2) fe++;
            int is_ab = 0; { int k = 0; while (k < m && (w >> k & 1)) k++; int ok = 1; for (int j = k; j < m; j++) if (w >> j & 1) ok = 0; is_ab = ok && k >= 2 && m - k >= 2; }
            if (ncyc < mincyc) { mincyc = ncyc; mincyc_words = 0; ab_attains = 0; }
            if (ncyc == mincyc) { mincyc_words++; if (is_ab) ab_attains = 1; if (!exw[0]) { for (int i = 0; i < m; i++) exw[i] = '0' + (w >> i & 1); exw[m] = 0; } }
            /* 1(f) */
            for (int i = 1; i < m - 1; i++) {
                if (!bal(bef[i][Z], H)) ff++;
                if (X[i] + Y[i] + WM[i] != H || 2 * (X[i] + WM[i] / 2) != 2 * bef[i][Z]) ff++;
                int ii = i + 1;
                if (m >= 8 && ii >= 3 && 2 * ii <= m + 1) { nxy++; if (!(2 * X[i] <= WM[i] && 2 * Y[i] <= WM[i])) fxy++; }
            }
            /* 1(g) */
            int ok = 1;
            for (int v = 0; v < n && ok; v++) { int has = 0; for (int x = 0; x < n; x++) if (x != v && bal(bef[v][x], H)) has = 1; if (!has) ok = 0; }
            unsigned seen = 1, fr = 1;
            while (fr) { int v = __builtin_ctz(fr); fr &= fr - 1; for (int x = 0; x < n; x++) if (x != v && bal(bef[v][x], H) && !(seen >> x & 1)) { seen |= 1u << x; fr |= 1u << x; } }
            if (seen != (1u << n) - 1) ok = 0;
            if (!ok) fg++;
        }
        printf("m=%2d (n=%2d): %lld strong words | 1(e) failures %lld; min #balanced arcs on the Hamiltonian cycle = %d, attained by %lld words "
               "(e.g. %s), attained by some 1^a0^b (a,b>=2): %d | 1(f) failures %lld | 2X,2Y<=W_M cases %lld, failures %lld | 1(g) failures %lld\n",
               m, n, words, fe, mincyc, mincyc_words, exw, ab_attains, ff, nxy, fxy, fg);
        fflush(stdout);
        tot_fail_e += fe; tot_fail_f += ff; tot_fail_g += fg; tot_fail_xy += fxy; tot_xy_cases += nxy;
    }
    printf("TOTAL failures: 1(e) %lld, 1(f) %lld, 2X/2Y<=W_M %lld of %lld cases, 1(g) %lld\n", tot_fail_e, tot_fail_f, tot_fail_xy, tot_xy_cases, tot_fail_g);
    /* the analytic threshold of 1(f) */
    int firstbad = -1;
    for (int mm = 9; mm <= 400; mm++) { double lhs = 2.0 * mm + __builtin_pow(2.0, (mm + 1) / 2.0), rhs = __builtin_pow(2.0, mm - 3); if (lhs > rhs) { firstbad = mm; break; } }
    printf("analytic condition 2m + 2^((m+1)/2) <= 2^(m-3): first failure for 9 <= m <= 400: %d (-1 = none); at m = 8: %s\n", firstbad,
           (16 + __builtin_pow(2.0, 4.5) <= 32) ? "holds" : "fails");
    return 0;
}
