/*
 * h21_achievable_w_n7.c — Exhaustive at n=7: find ALL achievable values of
 * w = (H-1)/2 = alpha_1 + 2*alpha_2 + 4*alpha_3 + ...
 *
 * Key question: which values of w are NEVER achievable at any n?
 * We know w=3 (H=7) and w=10 (H=21) are permanent gaps.
 * What pattern do they follow?
 *
 * Author: opus-2026-03-07-S43
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define N 7
#define NBITS (N*(N-1)/2)
#define TOTAL (1U << NBITS)

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
    int max_H = 0;
    int h_count[200]; /* count how many tournaments achieve each H */
    memset(h_count, 0, sizeof(h_count));

    unsigned char out[N];
    for (unsigned int T = 0; T < TOTAL; T++) {
        memset(out, 0, N);
        int pos = 0;
        for (int i = 0; i < N; i++)
            for (int j = i+1; j < N; j++) {
                if ((T >> pos) & 1)
                    out[i] |= (1 << j);
                else
                    out[j] |= (1 << i);
                pos++;
            }
        int H = hamiltonian_paths(out);
        if (H < 200) h_count[H]++;
        if (H > max_H) max_H = H;
    }

    printf("n=%d: max_H=%d\n", N, max_H);
    printf("\nAchievable H values:\n");
    for (int h = 1; h < 200; h += 2)
        if (h_count[h] > 0)
            printf("  H=%3d: %d tournaments\n", h, h_count[h]);

    printf("\nMissing odd H values in [1..%d]:\n", max_H < 200 ? max_H : 199);
    for (int h = 1; h < 200 && h <= max_H; h += 2)
        if (h_count[h] == 0)
            printf("  H=%d (w=%d) MISSING\n", h, (h-1)/2);

    printf("\nAchievable w=(H-1)/2 values:\n");
    for (int w = 0; w < 100; w++) {
        int h = 2*w + 1;
        if (h < 200 && h_count[h] > 0)
            printf("  w=%d ", w);
    }
    printf("\n\nMissing w values in [0..%d]:\n", (max_H-1)/2 < 100 ? (max_H-1)/2 : 99);
    for (int w = 0; w <= (max_H-1)/2 && w < 100; w++) {
        int h = 2*w + 1;
        if (h < 200 && h_count[h] == 0)
            printf("  w=%d (H=%d)\n", w, h);
    }

    return 0;
}
