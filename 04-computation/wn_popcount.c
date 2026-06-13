/*
 * wn_popcount.c — Compute W(n) using popcount-level DP (memory efficient)
 * opus-2026-03-15-S89c
 *
 * Processes masks level by level (by popcount). Only keeps current level
 * and next level in memory. Memory: 2 × C(n,n/2) × n × 16 bytes.
 *
 * For n=25: ~4.2 GB. For n=26: ~8.7 GB.
 *
 * Compile: gcc -O3 -o wn_popcount wn_popcount.c -lm
 * Usage: ./wn_popcount <n>
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

/* Count bits in a 32-bit mask */
static inline int popcount(unsigned int x) {
    return __builtin_popcount(x);
}

/* Enumerate all n-bit masks with exactly p bits set.
 * Returns array of masks and count. Caller must free. */
unsigned int *enum_masks(int n, int p, long long *count) {
    /* Count first */
    long long cnt = 0;
    unsigned int full = (1u << n) - 1;
    for (unsigned int m = 0; m <= full; m++) {
        if (popcount(m) == p) cnt++;
    }
    *count = cnt;
    unsigned int *arr = (unsigned int *)malloc(cnt * sizeof(unsigned int));
    long long idx = 0;
    for (unsigned int m = 0; m <= full; m++) {
        if (popcount(m) == p) arr[idx++] = m;
    }
    return arr;
}

/* Hash table: mask -> index in array */
/* For efficiency, use a direct-mapped array indexed by mask rank.
 * We precompute the rank of each mask within its popcount class.
 */

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    if (n < 2 || n > 28) {
        fprintf(stderr, "n must be 2..28\n");
        return 1;
    }

    unsigned int full = (1u << n) - 1;
    time_t t0 = time(NULL);

    /* For each popcount level, we need:
     * - Array of masks at this level
     * - dp values: u128 array of size count × n
     * We process level p to generate level p+1.
     */

    /* Start: level 1 (single vertices) */
    long long cnt_curr;
    unsigned int *masks_curr = enum_masks(n, 1, &cnt_curr);
    /* cnt_curr = n */
    u128 *dp_curr = (u128 *)calloc(cnt_curr * n, sizeof(u128));

    /* Initialize: dp[1<<v][v] = 1 for each v */
    for (int v = 0; v < n; v++) {
        /* Find index of mask (1<<v) in masks_curr */
        /* Since masks are enumerated in order, mask (1<<v) has index v */
        dp_curr[(long long)v * n + v] = 1;
    }

    fprintf(stderr, "n=%d, starting DP\n", n);

    for (int p = 1; p < n; p++) {
        /* Generate level p+1 from level p */
        long long cnt_next;
        unsigned int *masks_next = enum_masks(n, p + 1, &cnt_next);
        u128 *dp_next = (u128 *)calloc(cnt_next * n, sizeof(u128));

        if (!dp_next) {
            fprintf(stderr, "Failed to allocate dp_next (level %d, %lld masks, %.1f GB)\n",
                    p+1, cnt_next, (double)cnt_next * n * 16 / 1e9);
            free(dp_curr); free(masks_curr); free(masks_next);
            return 1;
        }

        /* Build mask -> index mapping for next level.
         * For efficiency, use a hash table or direct lookup.
         * Since masks_next is sorted, we can use binary search.
         * But that's O(log) per lookup. For speed, use a direct array.
         * For n≤28, masks fit in 28 bits. We can't afford a 2^28 array (256M entries).
         * Use a hash table instead.
         */

        /* Simple hash table: size = 2*cnt_next, linear probing */
        long long ht_size = cnt_next * 2;
        if (ht_size < 16) ht_size = 16;
        unsigned int *ht_keys = (unsigned int *)malloc(ht_size * sizeof(unsigned int));
        long long *ht_vals = (long long *)malloc(ht_size * sizeof(long long));
        memset(ht_keys, 0xFF, ht_size * sizeof(unsigned int)); /* 0xFFFFFFFF = empty */

        for (long long i = 0; i < cnt_next; i++) {
            unsigned int key = masks_next[i];
            long long h = (unsigned long long)key % ht_size;
            while (ht_keys[h] != 0xFFFFFFFF) {
                h = (h + 1) % ht_size;
            }
            ht_keys[h] = key;
            ht_vals[h] = i;
        }

        /* Process all masks at level p */
        for (long long mi = 0; mi < cnt_curr; mi++) {
            unsigned int mask = masks_curr[mi];
            long long base = mi * n;

            for (int v = 0; v < n; v++) {
                u128 cnt = dp_curr[base + v];
                if (cnt == 0) continue;
                if (!(mask & (1u << v))) continue;

                for (int u = 0; u < n; u++) {
                    if (mask & (1u << u)) continue;
                    if (u == v - 1) continue; /* skip unit descent */

                    u128 weight = (u == v + 1) ? 2 * cnt : cnt;
                    unsigned int new_mask = mask | (1u << u);

                    /* Look up index of new_mask in hash table */
                    long long h = (unsigned long long)new_mask % ht_size;
                    while (ht_keys[h] != new_mask) {
                        h = (h + 1) % ht_size;
                    }
                    long long ni = ht_vals[h];
                    dp_next[ni * n + u] += weight;
                }
            }
        }

        fprintf(stderr, "  level %d->%d: %lld masks (%lds)\n",
                p, p+1, cnt_next, time(NULL) - t0);

        /* Free current level, advance */
        free(dp_curr);
        free(masks_curr);
        free(ht_keys);
        free(ht_vals);

        dp_curr = dp_next;
        masks_curr = masks_next;
        cnt_curr = cnt_next;
    }

    /* Sum over endpoints at full mask (level n, single mask) */
    u128 W = 0;
    for (int v = 0; v < n; v++) {
        W += dp_curr[v];
    }

    printf("W(%d) = ", n);
    print_u128(W);
    printf("\n");

    fprintf(stderr, "Total time: %lds\n", time(NULL) - t0);

    free(dp_curr);
    free(masks_curr);
    return 0;
}
