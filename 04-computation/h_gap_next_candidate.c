/*
 * h_gap_next_candidate.c — Search for the NEXT permanent H-gap after 21.
 *
 * Strategy: exhaustive at n=7, sample at n=8,9,10 to find which odd H
 * values are NEVER achieved. Known permanent gaps: 7, 21.
 * Question: is 63 (=7*9) a gap? We know it fills at n=8 (227 in 600k).
 * What about other multiples of 7: 35, 49, 77, 91?
 *
 * Focus on: which odd H values in [1..300] are HARDEST to achieve?
 * Track the smallest n where each value first appears.
 *
 * Author: opus-2026-03-07-S43
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 8
#define MAX_H 300

/* Count Hamiltonian paths */
static long long hamiltonian_paths(const unsigned short out[N]) {
    static long long dp[1 << N][N];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++)
        dp[1 << v][v] = 1;
    int full = (1 << N) - 1;
    for (int mask = 1; mask <= full; mask++)
        for (int v = 0; v < N; v++) {
            if (!(mask & (1 << v)) || dp[mask][v] == 0) continue;
            unsigned short nexts = out[v] & ~mask;
            while (nexts) {
                int u = __builtin_ctz(nexts);
                nexts &= nexts - 1;
                dp[mask | (1 << u)][u] += dp[mask][v];
            }
        }
    long long total = 0;
    for (int v = 0; v < N; v++) total += dp[full][v];
    return total;
}

int main() {
    setvbuf(stdout, NULL, _IONBF, 0);

    /* Exhaustive enumeration at n=N */
    long long total = 1LL << (N * (N-1) / 2);
    int found[MAX_H + 1];
    memset(found, 0, sizeof(found));
    long long count_per_h[MAX_H + 1];
    memset(count_per_h, 0, sizeof(count_per_h));
    int max_h_seen = 0;

    unsigned short out[N];

    for (long long T = 0; T < total; T++) {
        /* Decode tournament from bit encoding */
        memset(out, 0, N * sizeof(unsigned short));
        int pos = 0;
        for (int i = 0; i < N; i++)
            for (int j = i+1; j < N; j++) {
                if ((T >> pos) & 1)
                    out[i] |= (1 << j);
                else
                    out[j] |= (1 << i);
                pos++;
            }

        long long H = hamiltonian_paths(out);
        if (H <= MAX_H) {
            found[(int)H] = 1;
            count_per_h[(int)H]++;
        }
        if (H > max_h_seen) max_h_seen = (int)H;

        if ((T + 1) % 50000000 == 0) {
            int missing = 0;
            for (int h = 1; h <= MAX_H; h += 2)
                if (!found[h]) missing++;
            fprintf(stderr, "Progress: %lld/%lld (%.1f%%), max_H=%d, missing_odd_in_[1,%d]=%d\n",
                    T+1, total, 100.0*(T+1)/total, max_h_seen, MAX_H, missing);
        }
    }

    printf("=== EXHAUSTIVE n=%d H-SPECTRUM ===\n", N);
    printf("Total tournaments: %lld\n", total);
    printf("Max H seen: %d\n", max_h_seen);

    printf("\nMissing odd values in [1, %d]:\n", MAX_H);
    int missing_count = 0;
    for (int h = 1; h <= MAX_H; h += 2) {
        if (!found[h]) {
            printf("  H=%d MISSING\n", h);
            missing_count++;
        }
    }
    printf("Total missing: %d\n", missing_count);

    printf("\nAchievable odd values with count:\n");
    for (int h = 1; h <= MAX_H; h += 2) {
        if (found[h]) {
            printf("  H=%d: %lld tournaments\n", h, count_per_h[h]);
        }
    }

    return 0;
}
