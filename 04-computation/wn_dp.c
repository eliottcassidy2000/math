/*
 * wn_dp.c — Compute W(n) = Σ_{NUD perms} 2^{adj1} using bitmask DP
 * opus-2026-03-15-S89c
 *
 * Compile: gcc -O3 -o wn_dp wn_dp.c
 * Usage: ./wn_dp <n>
 *
 * Uses __int128 for large values. For n≤26 this suffices.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef unsigned __int128 u128;

void print_u128(u128 x) {
    if (x == 0) { printf("0"); return; }
    char buf[40];
    int pos = 0;
    while (x > 0) {
        buf[pos++] = '0' + (int)(x % 10);
        x /= 10;
    }
    for (int i = pos - 1; i >= 0; i--)
        putchar(buf[i]);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    if (n < 2 || n > 26) {
        fprintf(stderr, "n must be 2..26\n");
        return 1;
    }

    long long full = (1LL << n) - 1;
    long long total_masks = full + 1;

    /* Allocate dp[mask][v] as flat array.
     * For n=24: 2^24 × 24 × 16 bytes = 6.4 GB — too much!
     * Need a hash map or sparse approach.
     */

    /* For n ≤ 22: 2^22 × 22 × 8 bytes = 739 MB — tight but possible with u64.
     * For n = 24: need ~6 GB. Use a different approach.
     *
     * Alternative: process masks by popcount. At each step, only keep
     * current and next popcount levels.
     *
     * At popcount p, the number of masks is C(n,p).
     * Max is C(24,12) ≈ 2.7M. Each has n entries.
     * 2.7M × 24 × 16 bytes ≈ 1 GB. Still big.
     *
     * Better: use a hash table. Only store nonzero entries.
     */

    /* Simple approach: use arrays, but only for n ≤ 23.
     * For n=24, use hash-based sparse DP.
     *
     * Actually, let me just use arrays for the full space.
     * For n=24: 2^24 = 16M masks. 16M × 24 = 384M entries.
     * At 8 bytes each (u64): 3 GB. With u128: 6 GB. Hmm.
     *
     * For values, u64 might overflow for n=24.
     * W(23) ≈ 2.8e19, so path weights at n=24 might exceed 2^63.
     * Need u128 or careful tracking.
     *
     * Let's use u64 and hope it doesn't overflow (values sum to W(n)
     * which is ~10^20 for n=24, but individual dp[mask][v] values
     * are much smaller — they're partial sums over subsets).
     *
     * Actually, the issue is that dp[(full,v)] can be up to W(n)/n ≈ 10^18.
     * u64 max is ~1.8×10^19. Should be fine for n=24.
     *
     * Memory: 2^24 × 24 × 8 = 3.2 GB. Tight.
     * Let's try with n ≤ 24 using arrays.
     */

    if (n > 24) {
        fprintf(stderr, "Array approach only supports n <= 24\n");
        return 1;
    }

    /* For n=24: 16M * 24 * 8 bytes = 3 GB */
    long long arr_size = (1LL << n) * n;
    unsigned long long *dp = (unsigned long long *)calloc(arr_size, sizeof(unsigned long long));
    if (!dp) {
        fprintf(stderr, "Failed to allocate %lld bytes\n", arr_size * 8);
        return 1;
    }

    /* Initialize single-vertex paths */
    for (int v = 0; v < n; v++) {
        dp[((1LL << v) * n) + v] = 1;
    }

    time_t t0 = time(NULL);
    long long progress_step = 1LL << (n > 10 ? n - 10 : 0);

    for (long long mask = 1; mask <= full; mask++) {
        if (mask % progress_step == 0) {
            double pct = (double)mask / full * 100;
            fprintf(stderr, "  %.1f%% (%lds)\n", pct, time(NULL) - t0);
        }

        long long base = mask * n;
        for (int v = 0; v < n; v++) {
            unsigned long long cnt = dp[base + v];
            if (cnt == 0) continue;

            for (int u = 0; u < n; u++) {
                if (mask & (1LL << u)) continue;
                if (u == v - 1) continue; /* unit descent: skip */

                unsigned long long weight = (u == v + 1) ? 2 * cnt : cnt;
                long long new_mask = mask | (1LL << u);
                dp[new_mask * n + u] += weight;
            }
        }
    }

    /* Sum over all endpoints at full mask */
    u128 W = 0;
    long long full_base = full * n;
    for (int v = 0; v < n; v++) {
        W += dp[full_base + v];
    }

    printf("W(%d) = ", n);
    print_u128(W);
    printf("\n");

    fprintf(stderr, "Total time: %lds\n", time(NULL) - t0);

    free(dp);
    return 0;
}
