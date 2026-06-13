/*
 * wn_rank2.c — Compute W(n), memory-optimized with packed endpoint storage
 * opus-2026-03-15-S89c
 *
 * Key optimization: only store values for endpoints that ARE set in the mask.
 * A mask with p bits has exactly p valid endpoints (not n).
 * This reduces memory by factor n/p at each level.
 *
 * At peak (n=26, level 13): C(26,13)*13*16 = 2.17 GB (vs 4.33 GB naive).
 * Two levels: 4.17 GB (vs 8.33 GB). Fits in 8 GB easily.
 *
 * Compile: gcc -O3 -o wn_rank2 wn_rank2.c
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

static inline int popcount(unsigned int x) {
    return __builtin_popcount(x);
}

static long long binom[31][31];

void init_binom() {
    for (int i = 0; i <= 30; i++) {
        binom[i][0] = 1;
        for (int j = 1; j <= i; j++)
            binom[i][j] = binom[i-1][j-1] + binom[i-1][j];
        for (int j = i+1; j <= 30; j++)
            binom[i][j] = 0;
    }
}

static inline long long mask_rank(unsigned int mask, int n, int p) {
    long long rank = 0;
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (mask & (1u << i)) {
            rank += binom[i][count + 1];
            count++;
        }
    }
    return rank;
}

static inline unsigned int mask_unrank(long long r, int n, int p) {
    unsigned int mask = 0;
    int remaining = p;
    for (int i = n - 1; i >= 0 && remaining > 0; i--) {
        long long c = binom[i][remaining];
        if (r >= c) {
            r -= c;
            mask |= (1u << i);
            remaining--;
        }
    }
    return mask;
}

/* Given a mask with p bits set, return the LOCAL index (0..p-1) of bit v within mask.
 * v must be set in mask. Returns popcount(mask & ((1<<v)-1)). */
static inline int local_index(unsigned int mask, int v) {
    return popcount(mask & ((1u << v) - 1));
}

/* Given a mask and a local index j (0..p-1), return the j-th set bit */
static inline int local_to_global(unsigned int mask, int j) {
    int count = 0;
    for (int i = 0; i < 32; i++) {
        if (mask & (1u << i)) {
            if (count == j) return i;
            count++;
        }
    }
    return -1; /* shouldn't happen */
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    if (n < 2 || n > 30) {
        fprintf(stderr, "n must be 2..30\n");
        return 1;
    }

    init_binom();
    time_t t0 = time(NULL);

    /* Level 1: p=1, each mask has 1 endpoint */
    int p_curr = 1;
    long long cnt_curr = binom[n][1]; /* = n */
    /* dp_curr: cnt_curr * p_curr entries */
    u128 *dp_curr = (u128 *)calloc(cnt_curr * p_curr, sizeof(u128));
    if (!dp_curr) { fprintf(stderr, "alloc fail\n"); return 1; }

    for (int v = 0; v < n; v++) {
        long long r = mask_rank(1u << v, n, 1);
        /* local index of v in mask (1<<v) is 0 */
        dp_curr[r * p_curr + 0] = 1;
    }

    fprintf(stderr, "n=%d, starting DP\n", n);

    for (int p = 1; p < n; p++) {
        int p_next = p + 1;
        long long cnt_next = binom[n][p_next];
        double mem_curr_gb = (double)cnt_curr * p * sizeof(u128) / 1e9;
        double mem_next_gb = (double)cnt_next * p_next * sizeof(u128) / 1e9;

        fprintf(stderr, "  level %d->%d: %lld -> %lld masks (%.1f + %.1f = %.1f GB) (%lds)...\n",
                p, p_next, cnt_curr, cnt_next,
                mem_curr_gb, mem_next_gb, mem_curr_gb + mem_next_gb,
                time(NULL) - t0);
        fflush(stderr);

        u128 *dp_next = (u128 *)calloc(cnt_next * p_next, sizeof(u128));
        if (!dp_next) {
            fprintf(stderr, "Failed to allocate dp_next (%.1f GB)\n", mem_next_gb);
            free(dp_curr);
            return 1;
        }

        for (long long mi = 0; mi < cnt_curr; mi++) {
            unsigned int mask = mask_unrank(mi, n, p);
            long long base_curr = mi * p;

            /* Iterate over valid endpoints (set bits in mask) */
            unsigned int tmp = mask;
            int li = 0;
            while (tmp) {
                int v = __builtin_ctz(tmp); /* lowest set bit */
                tmp &= tmp - 1; /* clear it */

                u128 cnt = dp_curr[base_curr + li];
                li++;
                if (cnt == 0) continue;

                for (int u = 0; u < n; u++) {
                    if (mask & (1u << u)) continue;
                    if (u == v - 1) continue;

                    u128 weight = (u == v + 1) ? 2 * cnt : cnt;
                    unsigned int new_mask = mask | (1u << u);
                    long long nr = mask_rank(new_mask, n, p_next);
                    int local_u = local_index(new_mask, u);
                    dp_next[nr * p_next + local_u] += weight;
                }
            }
        }

        fprintf(stderr, "    done (%lds)\n", time(NULL) - t0);
        free(dp_curr);
        dp_curr = dp_next;
        cnt_curr = cnt_next;
        p_curr = p_next;
    }

    /* Sum over all endpoints at full mask */
    u128 W = 0;
    for (int j = 0; j < n; j++) {
        W += dp_curr[j];
    }

    printf("W(%d) = ", n);
    print_u128(W);
    printf("\n");

    fprintf(stderr, "Total time: %lds\n", time(NULL) - t0);

    free(dp_curr);
    return 0;
}
