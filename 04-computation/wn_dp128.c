/*
 * wn_dp128.c — Compute W(n) using __int128 to avoid overflow
 * opus-2026-03-15-S89c
 *
 * Compile: gcc -O3 -o wn_dp128 wn_dp128.c
 * Usage: ./wn_dp128 <n>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef unsigned __int128 u128;

void print_u128(u128 x) {
    if (x == 0) { printf("0"); return; }
    char buf[50];
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
    long long arr_size = (1LL << n) * (long long)n;
    long long byte_size = arr_size * (long long)sizeof(u128);

    fprintf(stderr, "n=%d: allocating %.1f GB for dp array\n", n, byte_size / 1e9);

    u128 *dp = (u128 *)calloc(arr_size, sizeof(u128));
    if (!dp) {
        fprintf(stderr, "Failed to allocate %lld bytes\n", byte_size);
        return 1;
    }

    /* Initialize single-vertex paths */
    for (int v = 0; v < n; v++) {
        dp[((long long)(1LL << v)) * n + v] = 1;
    }

    time_t t0 = time(NULL);
    long long progress_step = 1LL << (n > 6 ? n - 6 : 0);

    for (long long mask = 1; mask <= full; mask++) {
        if (mask % progress_step == 0) {
            double pct = (double)mask / full * 100;
            fprintf(stderr, "  %.1f%% (%lds)\n", pct, time(NULL) - t0);
        }

        long long base = mask * n;
        for (int v = 0; v < n; v++) {
            u128 cnt = dp[base + v];
            if (cnt == 0) continue;

            for (int u = 0; u < n; u++) {
                if (mask & (1LL << u)) continue;
                if (u == v - 1) continue;

                u128 weight = (u == v + 1) ? 2 * cnt : cnt;
                long long new_mask = mask | (1LL << u);
                dp[(long long)new_mask * n + u] += weight;
            }
        }
    }

    u128 W = 0;
    long long full_base = (long long)full * n;
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
