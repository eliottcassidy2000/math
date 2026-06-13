/*
 * wn_rank.c — Compute W(n) using popcount-level DP with rank-based indexing
 * opus-2026-03-15-S89c
 *
 * Key improvement over wn_popcount.c: No hash table needed.
 * Each mask's position in its popcount class is computed via combinatorial
 * ranking (combinadic). This saves ~50% memory at peak levels.
 *
 * Memory: C(n,p)*n*16 + C(n,p+1)*n*16 bytes (two levels only)
 * For n=26: peak at level 12->13: 9657700*26*16 + 10400600*26*16 = 8.3 GB
 *   Still tight. But we can process in a streaming fashion:
 *   - Allocate ONLY next level array
 *   - Read current level, generate next, free current
 *   - At peak: 4.3 GB for the larger level
 *
 * Actually: we need both curr and next simultaneously during transition.
 * Peak = C(n,12)*n*16 + C(n,13)*n*16 = 4.0 + 4.3 = 8.3 GB.
 * With 8 GB RAM, this might just barely work with some swap.
 *
 * Optimization: use u64 where possible, switch to u128 only when needed.
 * For n=26, values overflow u64 around level 18-19. Use u128 throughout for safety.
 *
 * Compile: gcc -O3 -o wn_rank wn_rank.c
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

/* Precomputed binomial coefficients C(n,k) for n<=30, k<=30 */
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

/* Combinatorial rank of mask among all n-bit masks with exactly p bits set.
 * Uses the combinadic: rank = sum over set bits of C(bit_position, count).
 * Masks are ranked in increasing order. */
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

/* Unrank: given rank r among C(n,p) masks, return the mask */
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

    /* Level 1: single vertices */
    long long cnt_curr = n;
    /* dp_curr[rank * n + v] = count of paths ending at v through mask-of-rank-rank */
    u128 *dp_curr = (u128 *)calloc(cnt_curr * n, sizeof(u128));
    if (!dp_curr) {
        fprintf(stderr, "Failed to allocate initial dp\n");
        return 1;
    }

    /* Initialize: mask = 1<<v has rank v (since C(v,1)=v for the combinadic) */
    for (int v = 0; v < n; v++) {
        long long r = mask_rank(1u << v, n, 1);
        dp_curr[r * n + v] = 1;
    }

    fprintf(stderr, "n=%d, starting DP\n", n);

    for (int p = 1; p < n; p++) {
        long long cnt_next = binom[n][p + 1];
        double mem_gb = (double)(cnt_curr + cnt_next) * n * sizeof(u128) / 1e9;

        fprintf(stderr, "  level %d->%d: %lld -> %lld masks (%.1f GB) (%lds)...\n",
                p, p+1, cnt_curr, cnt_next, mem_gb, time(NULL) - t0);
        fflush(stderr);

        u128 *dp_next = (u128 *)calloc(cnt_next * n, sizeof(u128));
        if (!dp_next) {
            fprintf(stderr, "Failed to allocate dp_next (%.1f GB)\n",
                    (double)cnt_next * n * sizeof(u128) / 1e9);
            free(dp_curr);
            return 1;
        }

        /* Iterate over all masks at current level */
        for (long long mi = 0; mi < cnt_curr; mi++) {
            unsigned int mask = mask_unrank(mi, n, p);
            long long base_curr = mi * n;

            for (int v = 0; v < n; v++) {
                u128 cnt = dp_curr[base_curr + v];
                if (cnt == 0) continue;
                if (!(mask & (1u << v))) continue;

                for (int u = 0; u < n; u++) {
                    if (mask & (1u << u)) continue;
                    if (u == v - 1) continue;

                    u128 weight = (u == v + 1) ? 2 * cnt : cnt;
                    unsigned int new_mask = mask | (1u << u);
                    long long nr = mask_rank(new_mask, n, p + 1);
                    dp_next[nr * n + u] += weight;
                }
            }
        }

        fprintf(stderr, "    done (%lds)\n", time(NULL) - t0);

        free(dp_curr);
        dp_curr = dp_next;
        cnt_curr = cnt_next;
    }

    /* Sum over endpoints at full mask */
    u128 W = 0;
    for (int v = 0; v < n; v++) {
        W += dp_curr[v];
    }

    printf("W(%d) = ", n);
    print_u128(W);
    printf("\n");

    fprintf(stderr, "Total time: %lds\n", time(NULL) - t0);

    free(dp_curr);
    return 0;
}
