/*
 * h21_w_n8_sample.c — Sample n=8 tournaments to find achievable w values.
 * Check which n=7 gaps persist at n=8.
 *
 * Author: opus-2026-03-07-S43
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define N 8
#define TRIALS 50000000

static unsigned int rng_state = 12345;
static inline unsigned int xorshift() {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

static int hamiltonian_paths(unsigned char out[N]) {
    int dp[1 << N][N];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++)
        dp[1 << v][v] = 1;
    int full = (1 << N) - 1;
    for (int mask = 1; mask <= full; mask++) {
        for (int v = 0; v < N; v++) {
            if (!(mask & (1 << v)) || dp[mask][v] == 0) continue;
            for (int u = 0; u < N; u++) {
                if (mask & (1 << u)) continue;
                if ((out[v] >> u) & 1)
                    dp[mask | (1 << u)][u] += dp[mask][v];
            }
        }
    }
    int total = 0;
    for (int v = 0; v < N; v++)
        total += dp[full][v];
    return total;
}

int main() {
    setvbuf(stdout, NULL, _IONBF, 0);

    /* Track which H values we've seen */
    int seen[700]; /* max H at n=8 is 661 */
    memset(seen, 0, sizeof(seen));

    unsigned char out[N];
    int max_H = 0;

    for (long t = 0; t < TRIALS; t++) {
        memset(out, 0, N);
        for (int i = 0; i < N; i++)
            for (int j = i+1; j < N; j++) {
                if (xorshift() & 1)
                    out[i] |= (1 << j);
                else
                    out[j] |= (1 << i);
            }

        int H = hamiltonian_paths(out);
        if (H < 700) seen[H] = 1;
        if (H > max_H) max_H = H;

        if ((t + 1) % 10000000 == 0) {
            /* Count missing odd values */
            int missing = 0;
            for (int h = 1; h < 700 && h <= max_H; h += 2)
                if (!seen[h]) missing++;
            fprintf(stderr, "Progress: %ldM, max_H=%d, missing_odd=%d\n",
                    (t+1)/1000000, max_H, missing);
        }
    }

    printf("n=%d (%ldM samples): max_H=%d\n", N, TRIALS/1000000, max_H);
    printf("\nMissing odd H values in [1..%d]:\n", max_H < 700 ? max_H : 699);
    for (int h = 1; h < 700 && h <= max_H; h += 2)
        if (!seen[h])
            printf("  H=%d (w=%d)\n", h, (h-1)/2);

    /* Check n=7 gaps */
    int n7_gaps[] = {7, 21, 63, 107, 119, 149, 161, 163, 165, 167, 169, 173,
                     177, 179, 181, 183, 185, 187};
    printf("\nn=7 gaps that PERSIST at n=8:\n");
    for (int i = 0; i < 18; i++) {
        int h = n7_gaps[i];
        if (h < 700 && !seen[h])
            printf("  H=%d STILL MISSING at n=8\n", h);
        else
            printf("  H=%d FILLED at n=8\n", h);
    }

    return 0;
}
