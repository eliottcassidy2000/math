/*
 * ocr_n9_n10_sampled.c — kind-pasteur-2026-03-21-S16c
 *
 * HIGH-PRECISION SAMPLING of OCR at n=9,10,11 to test the 120/131 plateau.
 * n=9: 2^36 ~ 69 billion tournaments (exhaustive impossible).
 * Strategy: sample 10M random tournaments, compute R^2 to ~6 decimal places.
 * If R^2 matches 0.916030534351... (= 120/131), the plateau extends.
 *
 * Uses Mersenne Twister via stdlib rand() with seed.
 *
 * Compile: gcc -O3 -o ocr_n9_n10_sampled ocr_n9_n10_sampled.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_N 11

typedef unsigned long long u64;

static int adj[MAX_N][MAX_N];

static int ham_paths(int n) {
    int dp[1 << MAX_N][MAX_N];
    int full = (1 << n) - 1;
    memset(dp, 0, sizeof(int) * (1 << n) * n);

    for (int v = 0; v < n; v++)
        dp[1 << v][v] = 1;

    for (int mask = 1; mask < (1 << n); mask++) {
        for (int v = 0; v < n; v++) {
            if (!(mask & (1 << v))) continue;
            int c = dp[mask][v];
            if (c == 0) continue;
            for (int w = 0; w < n; w++) {
                if (mask & (1 << w)) continue;
                if (adj[v][w])
                    dp[mask | (1 << w)][w] += c;
            }
        }
    }

    int total = 0;
    for (int v = 0; v < n; v++)
        total += dp[full][v];
    return total;
}

/* xorshift128+ for fast random numbers */
static u64 rng_s[2] = {0x12345678deadbeefULL, 0xfedcba9876543210ULL};

static u64 xorshift128plus(void) {
    u64 s1 = rng_s[0];
    u64 s0 = rng_s[1];
    rng_s[0] = s0;
    s1 ^= s1 << 23;
    rng_s[1] = s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26);
    return rng_s[1] + s0;
}

int main(int argc, char **argv) {
    int n = 9;
    long samples = 2000000;

    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) samples = atol(argv[2]);
    if (n < 3 || n > 11) { printf("n must be 3..11\n"); return 1; }

    int m = n * (n - 1) / 2;
    printf("OCR SAMPLING: n=%d, %ld samples, m=%d arcs\n\n", n, samples, m);
    fflush(stdout);

    /* Arc map */
    int arc_i[66], arc_j[66]; /* max C(11,2) = 55 */
    int pos = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) {
            arc_i[pos] = i;
            arc_j[pos] = j;
            pos++;
        }

    /* Running statistics (Welford's algorithm for numerical stability) */
    double mean_H = 0, M2_H = 0;
    double mean_S2 = 0, M2_S2 = 0;
    double mean_HS2 = 0;  /* For covariance */
    long count = 0;

    /* Also accumulate raw sums for cross-check */
    double sum_H = 0, sum_S2 = 0;
    double sum_H2 = 0, sum_S22 = 0, sum_HS2 = 0;

    time_t t0 = time(NULL);

    for (long iter = 0; iter < samples; iter++) {
        /* Generate random tournament */
        u64 bits = xorshift128plus();
        /* For n > 10, need more random bits */
        u64 bits2 = (m > 64) ? xorshift128plus() : 0;

        memset(adj, 0, sizeof(adj));
        for (int p = 0; p < m; p++) {
            int bit;
            if (p < 64) bit = (bits >> p) & 1;
            else bit = (bits2 >> (p - 64)) & 1;

            if (bit) adj[arc_i[p]][arc_j[p]] = 1;
            else adj[arc_j[p]][arc_i[p]] = 1;
        }

        /* Scores and S2 */
        double S2 = 0;
        double half = (n - 1) / 2.0;
        for (int i = 0; i < n; i++) {
            int s = 0;
            for (int j = 0; j < n; j++) s += adj[i][j];
            S2 += (s - half) * (s - half);
        }

        int H = ham_paths(n);

        /* Accumulate */
        count++;
        sum_H += H;
        sum_S2 += S2;
        sum_H2 += (double)H * H;
        sum_S22 += S2 * S2;
        sum_HS2 += H * S2;

        if (count % 500000 == 0) {
            /* Compute running R^2 */
            double mH = sum_H / count;
            double mS = sum_S2 / count;
            double vH = sum_H2 / count - mH * mH;
            double vS = sum_S22 / count - mS * mS;
            double cov = sum_HS2 / count - mH * mS;
            double r2 = (cov * cov) / (vH * vS);

            time_t now = time(NULL);
            printf("  %ld samples: R^2 = %.12f  (1-R^2 = %.12f)  [%.0fs]\n",
                   count, r2, 1.0 - r2, difftime(now, t0));
            fflush(stdout);
        }
    }

    /* Final result */
    double mH = sum_H / count;
    double mS = sum_S2 / count;
    double vH = sum_H2 / count - mH * mH;
    double vS = sum_S22 / count - mS * mS;
    double cov = sum_HS2 / count - mH * mS;
    double r2 = (cov * cov) / (vH * vS);

    printf("\nFINAL (n=%d, %ld samples):\n", n, count);
    printf("  E[H]  = %.6f (expected: %d!/2^%d = %.6f)\n",
           mH, n, n-1, tgamma(n+1) / (1 << (n-1)));
    printf("  E[S2] = %.6f (expected: n(n-1)/4 = %.6f)\n",
           mS, n * (n-1) / 4.0);
    printf("  Var(H)  = %.6f\n", vH);
    printf("  Var(S2) = %.6f\n", vS);
    printf("  Cov(H,S2) = %.6f\n", cov);
    printf("  R^2(S2,H) = %.15f\n", r2);
    printf("  1 - R^2   = %.15f\n", 1.0 - r2);
    printf("  120/131   = %.15f\n", 120.0 / 131.0);
    printf("  Difference from 120/131: %.2e\n", fabs(r2 - 120.0/131.0));

    time_t t1 = time(NULL);
    printf("\n  Total time: %.0f seconds\n", difftime(t1, t0));

    return 0;
}
