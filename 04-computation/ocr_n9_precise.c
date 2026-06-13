/*
 * ocr_n9_precise.c — kind-pasteur-2026-03-21-S16c
 *
 * Precise OCR sampling at n=9 using __int128 accumulators.
 * H values up to ~15745, S2 up to ~32.
 * H^2 up to ~248M. Over 2M samples, sum_H2 ~ 5e14 — fits in int64.
 * But sum_H^2 / N^2 subtraction needs care. Use __int128.
 *
 * Compile: gcc -O3 -o ocr_n9_precise ocr_n9_precise.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define N_MAX 11

typedef long long i64;
typedef __int128 i128;
typedef unsigned long long u64;

static int adj[N_MAX][N_MAX];

static int ham_paths(int n) {
    static int dp[1 << N_MAX][N_MAX];
    int full = (1 << n) - 1;
    memset(dp, 0, sizeof(dp));  /* Clear ENTIRE array */
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (int mask = 1; mask < (1 << n); mask++)
        for (int v = 0; v < n; v++) {
            if (!(mask & (1 << v))) continue;
            int c = dp[mask][v];
            if (c == 0) continue;
            for (int w = 0; w < n; w++) {
                if (mask & (1 << w)) continue;
                if (adj[v][w]) dp[mask | (1 << w)][w] += c;
            }
        }
    int total = 0;
    for (int v = 0; v < n; v++) total += dp[full][v];
    return total;
}

static u64 rng_s[2] = {0x12345678deadbeefULL, 0xfedcba9876543210ULL};
static u64 xorshift128plus(void) {
    u64 s1 = rng_s[0], s0 = rng_s[1];
    rng_s[0] = s0;
    s1 ^= s1 << 23;
    rng_s[1] = s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26);
    return rng_s[1] + s0;
}

static void print128(i128 val) {
    if (val < 0) { printf("-"); val = -val; }
    if (val == 0) { printf("0"); return; }
    char buf[50]; int pos = 0;
    while (val > 0) { buf[pos++] = '0' + (int)(val % 10); val /= 10; }
    for (int i = pos-1; i >= 0; i--) printf("%c", buf[i]);
}

int main(int argc, char **argv) {
    int n = 9;
    long samples = 2000000;
    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) samples = atol(argv[2]);

    int m = n * (n-1) / 2;
    printf("PRECISE OCR SAMPLING: n=%d, %ld samples, m=%d arcs\n", n, samples, m);
    fflush(stdout);

    int arc_i[66], arc_j[66];
    int pos = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) { arc_i[pos] = i; arc_j[pos] = j; pos++; }

    /* Integer S2: compute 4*S2 = sum (2*s_i - (n-1))^2 to keep everything integer */
    /* This works for both odd and even n since 2*s_i and n-1 are both integers */

    i128 sum_H = 0, sum_4S2 = 0;
    i128 sum_H2 = 0, sum_16S22 = 0, sum_4HS2 = 0;
    long count = 0;
    time_t t0 = time(NULL);

    for (long iter = 0; iter < samples; iter++) {
        u64 bits = xorshift128plus();

        memset(adj, 0, sizeof(int) * n * N_MAX);
        for (int p = 0; p < m; p++) {
            if ((bits >> p) & 1) adj[arc_i[p]][arc_j[p]] = 1;
            else adj[arc_j[p]][arc_i[p]] = 1;
        }

        /* 4*S2 = sum (2*s_i - (n-1))^2 */
        i64 S2_4 = 0;
        for (int i = 0; i < n; i++) {
            int s = 0;
            for (int j = 0; j < n; j++) s += adj[i][j];
            int ds = 2*s - (n-1);
            S2_4 += (i64)ds * ds;
        }

        i64 H = ham_paths(n);

        sum_H += H;
        sum_4S2 += S2_4;
        sum_H2 += (i128)H * H;
        sum_16S22 += (i128)S2_4 * S2_4;
        sum_4HS2 += (i128)H * S2_4;
        count++;

        if (count % 500000 == 0) {
            /* Compute R^2 using __int128 */
            i128 NN = count;
            i128 cov = NN * sum_4HS2 - sum_H * sum_4S2;
            i128 vH = NN * sum_H2 - sum_H * sum_H;
            i128 vS = NN * sum_16S22 - sum_4S2 * sum_4S2;

            /* R^2 = cov^2 / (vH * vS) — compute as double */
            double cov_d = (double)cov;
            double vH_d = (double)vH;
            double vS_d = (double)vS;
            double r2 = (cov_d * cov_d) / (vH_d * vS_d);

            time_t now = time(NULL);
            printf("  %ld: R^2=%.12f 1-R^2=%.12f [%.0fs]\n",
                   count, r2, 1.0-r2, difftime(now, t0));
            fflush(stdout);
        }
    }

    /* Final */
    i128 NN = count;
    i128 cov = NN * sum_4HS2 - sum_H * sum_4S2;
    i128 vH = NN * sum_H2 - sum_H * sum_H;
    i128 vS = NN * sum_16S22 - sum_4S2 * sum_4S2;

    double cov_d = (double)cov;
    double vH_d = (double)vH;
    double vS_d = (double)vS;
    double r2 = (cov_d * cov_d) / (vH_d * vS_d);

    printf("\nFINAL n=%d, %ld samples:\n", n, count);
    printf("  sum_H  = "); print128(sum_H); printf("\n");
    printf("  sum_4S2 = "); print128(sum_4S2); printf("\n");
    printf("  R^2 = %.15f\n", r2);
    printf("  1-R^2 = %.15f\n", 1.0-r2);
    printf("  120/131 = %.15f\n", 120.0/131.0);
    printf("  Diff = %.2e\n", fabs(r2 - 120.0/131.0));
    printf("  E[H] = %.2f (expected %.2f)\n",
           (double)sum_H/count, tgamma(n+1)/(1<<(n-1)));

    time_t t1 = time(NULL);
    printf("  Time: %.0fs\n", difftime(t1, t0));
    return 0;
}
