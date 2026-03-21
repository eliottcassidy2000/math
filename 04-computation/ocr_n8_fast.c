/*
 * ocr_n8_fast.c — kind-pasteur-2026-03-21-S16
 *
 * Compute EXACT OCR(8) via exhaustive enumeration of 2^28 = 268,435,456 tournaments.
 * Uses Held-Karp DP. Estimated runtime: ~20-60 minutes on modern CPU.
 *
 * We accumulate: sum(H), sum(S2), sum(H^2), sum(S2^2), sum(H*S2)
 * Using __int128 to prevent overflow (H values up to ~661 at n=8).
 *
 * Compile: gcc -O3 -o ocr_n8_fast ocr_n8_fast.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N 8
#define M 28  /* C(8,2) */
#define FULL ((1 << N) - 1)

typedef unsigned long long u64;
typedef __int128 i128;

static int adj[N][N];

static int ham_paths(void) {
    static int dp[1 << N][N];
    memset(dp, 0, sizeof(dp));

    for (int v = 0; v < N; v++)
        dp[1 << v][v] = 1;

    for (int mask = 1; mask < (1 << N); mask++) {
        for (int v = 0; v < N; v++) {
            if (!(mask & (1 << v))) continue;
            int c = dp[mask][v];
            if (c == 0) continue;
            for (int w = 0; w < N; w++) {
                if (mask & (1 << w)) continue;
                if (adj[v][w])
                    dp[mask | (1 << w)][w] += c;
            }
        }
    }

    int total = 0;
    for (int v = 0; v < N; v++)
        total += dp[FULL][v];
    return total;
}

/* Print __int128 */
static void print128(i128 val) {
    if (val < 0) {
        printf("-");
        val = -val;
    }
    if (val == 0) {
        printf("0");
        return;
    }
    char buf[40];
    int pos = 0;
    while (val > 0) {
        buf[pos++] = '0' + (int)(val % 10);
        val /= 10;
    }
    for (int i = pos - 1; i >= 0; i--)
        printf("%c", buf[i]);
}

int main(void) {
    printf("OCR(8) EXACT COMPUTATION\n");
    printf("Enumerating all 2^%d = %u tournaments on %d vertices...\n\n", M, 1u << M, N);
    fflush(stdout);

    /* Arc map */
    int arc_i[M], arc_j[M];
    int pos = 0;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++) {
            arc_i[pos] = i;
            arc_j[pos] = j;
            pos++;
        }

    i128 sum_H = 0, sum_S2 = 0;
    i128 sum_H2 = 0, sum_S22 = 0, sum_HS2 = 0;
    u64 total = 1ULL << M;

    time_t t0 = time(NULL);

    for (u64 bits = 0; bits < total; bits++) {
        /* Build adj */
        memset(adj, 0, sizeof(adj));
        for (int p = 0; p < M; p++) {
            if (bits & (1ULL << p))
                adj[arc_i[p]][arc_j[p]] = 1;
            else
                adj[arc_j[p]][arc_i[p]] = 1;
        }

        /* Scores and S2 */
        int S2 = 0;
        for (int i = 0; i < N; i++) {
            int s = 0;
            for (int j = 0; j < N; j++) s += adj[i][j];
            /* (n-1)/2 = 3.5, but S2 = sum(s - 3.5)^2 = sum(2s-7)^2 / 4 */
            int ds = 2*s - 7;  /* 2*(s - 3.5) */
            S2 += ds * ds;
        }
        /* S2 is now 4 * actual_S2. Use this integer form. */

        int H = ham_paths();

        sum_H += H;
        sum_S2 += S2;
        sum_H2 += (i128)H * H;
        sum_S22 += (i128)S2 * S2;
        sum_HS2 += (i128)H * S2;

        if ((bits & 0x3FFFFFF) == 0 && bits > 0) {
            time_t now = time(NULL);
            double elapsed = difftime(now, t0);
            double frac = (double)bits / total;
            double remaining = elapsed / frac * (1 - frac);
            printf("  %llu/%llu (%.1f%%) elapsed %.0fs remaining ~%.0fs\n",
                   (unsigned long long)bits, (unsigned long long)total,
                   100.0*frac, elapsed, remaining);
            fflush(stdout);
        }
    }

    time_t t1 = time(NULL);
    printf("\nDone in %.0f seconds.\n\n", difftime(t1, t0));

    printf("Results (S2 scaled by 4 to keep integer):\n");
    printf("  sum(H)    = "); print128(sum_H); printf("\n");
    printf("  sum(4S2)  = "); print128(sum_S2); printf("\n");
    printf("  sum(H^2)  = "); print128(sum_H2); printf("\n");
    printf("  sum(16S2^2) = "); print128(sum_S22); printf("\n");
    printf("  sum(4HS2) = "); print128(sum_HS2); printf("\n");

    /* R^2 = [N*sum_HS2 - sum_H*sum_S2]^2 / [(N*sum_S22 - sum_S2^2)(N*sum_H2 - sum_H^2)] */
    /* Note: S2 is multiplied by 4. To get true R^2:
     * True Cov(H, S2) = (sum_HS2/4 - sum_H*sum_S2/4/N)/N
     * = Cov(H, 4*S2) / 4
     * R^2 = Cov(H, 4S2)^2 / (Var(4S2)*Var(H))
     * = Cov(H, S2)^2 * 16 / (16*Var(S2)*Var(H))
     * = Cov(H, S2)^2 / (Var(S2)*Var(H))
     * So R^2(H, 4S2) = R^2(H, S2). Scale doesn't matter! */

    i128 NN = (i128)total;
    i128 cov = NN * sum_HS2 - sum_H * sum_S2;
    i128 vH = NN * sum_H2 - sum_H * sum_H;
    i128 vS = NN * sum_S22 - sum_S2 * sum_S2;

    printf("\n  N*Cov  = "); print128(cov); printf("\n");
    printf("  N2*VarH = "); print128(vH); printf("\n");
    printf("  N2*VarS = "); print128(vS); printf("\n");

    /* Approximate R^2 */
    double cov_d = (double)cov;
    double vH_d = (double)vH;
    double vS_d = (double)vS;
    double r2 = (cov_d * cov_d) / (vH_d * vS_d);
    printf("\n  R^2(S2, H) = %.12f\n", r2);
    printf("  1 - R^2    = %.12f\n", 1.0 - r2);

    /* For exact rational: need to compute gcd(cov^2, vH*vS) */
    /* This is too large for __int128. Report the raw values for Python post-processing. */
    printf("\n  (Use Python to compute exact rational from the sums above.)\n");

    return 0;
}
