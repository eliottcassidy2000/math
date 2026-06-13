/*
 * h21_min_h_cycle_rich_n10.c — Sample min H for cycle-rich n=10 tournaments.
 * Uses the OCF formula: H = I(Omega(T), 2) where Omega(T) is the odd-cycle graph.
 *
 * Since computing I(G,2) exactly for Omega is expensive, we use the formula:
 * H = sum over independent sets S of Omega(T) of 2^|S|
 * = sum_k alpha_k * 2^k where alpha_k = # independent sets of size k in Omega(T)
 *
 * For n=10, Omega has at most C(10,3)=120 vertices (3-cycles) plus higher odd cycles.
 * But 3-cycles dominate. We'll compute alpha_k for the 3-cycle subgraph of Omega.
 *
 * Actually, Omega(T) includes ALL odd cycles. For lower bound on H, just counting
 * 3-cycles and their independence structure suffices (I(subgraph, 2) <= I(Omega, 2)).
 * Wait, independence polynomial is for INDEPENDENT sets. More edges in Omega means
 * FEWER independent sets, so I(Omega, 2) <= I(3-cycle-subgraph, 2).
 *
 * For an UPPER bound (or exact), we'd need all odd cycles. Let's compute exact H
 * for small samples using the Hamiltonian path count approach.
 *
 * H(T) = number of Hamiltonian paths in T.
 *
 * Author: opus-2026-03-07-S43
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 10

static unsigned int rng_state = 77777;
static inline unsigned int xorshift() {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

static inline void random_tournament(unsigned short out[N]) {
    memset(out, 0, N * sizeof(unsigned short));
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++) {
            if (xorshift() & 1)
                out[i] |= (1 << j);
            else
                out[j] |= (1 << i);
        }
}

static inline int is_cycle_rich(const unsigned short out[N]) {
    unsigned short all = (1 << N) - 1;
    for (int i = 0; i < N; i++) {
        int s = __builtin_popcount(out[i]);
        if (s == 0 || s == N-1) return 0;
    }
    for (int v = 0; v < N; v++) {
        unsigned short ov = out[v];
        unsigned short iv = (~ov) & all & ~(1 << v);
        int found = 0;
        unsigned short ov_copy = ov;
        while (ov_copy && !found) {
            int a = __builtin_ctz(ov_copy);
            ov_copy &= ov_copy - 1;
            if (out[a] & iv & ~(1 << a))
                found = 1;
        }
        if (!found) return 0;
    }
    return 1;
}

/* Count Hamiltonian paths using DP on bitmask */
static long long hamiltonian_paths(const unsigned short out[N]) {
    /* dp[mask][v] = # Hamiltonian paths ending at v using vertex set mask */
    static long long dp[1 << N][N];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++)
        dp[1 << v][v] = 1;
    int full = (1 << N) - 1;
    for (int mask = 1; mask <= full; mask++) {
        for (int v = 0; v < N; v++) {
            if (!(mask & (1 << v)) || dp[mask][v] == 0) continue;
            unsigned short nexts = out[v] & ~mask;
            while (nexts) {
                int u = __builtin_ctz(nexts);
                nexts &= nexts - 1;
                dp[mask | (1 << u)][u] += dp[mask][v];
            }
        }
    }
    long long total = 0;
    for (int v = 0; v < N; v++) total += dp[full][v];
    return total;
}

/* Count 3-cycles */
static int count_3cycles(const unsigned short out[N]) {
    int count = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            for (int c = b+1; c < N; c++) {
                int ab = (out[a] >> b) & 1;
                int bc = (out[b] >> c) & 1;
                int ca = (out[c] >> a) & 1;
                if ((ab && bc && ca) || (!ab && !bc && !ca))
                    count++;
            }
    return count;
}

/* Max matching of 3-cycles (greedy) */
static int max_matching(const unsigned short out[N]) {
    int used = 0, mm = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            for (int c = b+1; c < N; c++) {
                int m = (1<<a)|(1<<b)|(1<<c);
                if (used & m) continue;
                int ab = (out[a] >> b) & 1;
                int bc = (out[b] >> c) & 1;
                int ca = (out[c] >> a) & 1;
                if ((ab && bc && ca) || (!ab && !bc && !ca)) {
                    used |= m;
                    mm++;
                }
            }
    return mm;
}

int main() {
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    long long trials = 10000000LL; /* 10M — each H computation is expensive */
    long long cycle_rich = 0;
    long long min_h = 999999999LL;
    long long min_h_t3 = 0, min_h_mm = 0;
    int min_h_scores[N];

    /* Track H distribution */
    int h_hist[200]; /* count of H values < 200 */
    memset(h_hist, 0, sizeof(h_hist));
    long long h_ge200 = 0;

    unsigned short out[N];

    for (long long t = 0; t < trials; t++) {
        random_tournament(out);
        if (!is_cycle_rich(out)) continue;
        cycle_rich++;

        long long H = hamiltonian_paths(out);
        if (H < min_h) {
            min_h = H;
            min_h_t3 = count_3cycles(out);
            min_h_mm = max_matching(out);
            for (int i = 0; i < N; i++)
                min_h_scores[i] = __builtin_popcount(out[i]);
            printf("NEW MIN H=%lld at cr=%lld (t3=%lld, mm=%lld, scores=",
                   H, cycle_rich, min_h_t3, min_h_mm);
            for (int i = 0; i < N; i++) printf("%d ", min_h_scores[i]);
            printf(")\n");
        }

        if (H < 200)
            h_hist[(int)H]++;
        else
            h_ge200++;

        if ((t + 1) % 2000000 == 0) {
            fprintf(stderr, "t=%lldM, cr=%lld, min_H=%lld\n",
                    (t+1)/1000000, cycle_rich, min_h);
        }
    }

    printf("\n=== n=%d RESULTS ===\n", N);
    printf("Trials: %lld\nCycle-rich: %lld\n", trials, cycle_rich);
    printf("Min H: %lld (t3=%lld, mm=%lld)\n", min_h, min_h_t3, min_h_mm);
    printf("Min H scores: ");
    for (int i = 0; i < N; i++) printf("%d ", min_h_scores[i]);
    printf("\n");
    printf("H distribution (small values):\n");
    for (int h = 1; h < 200; h += 2)
        if (h_hist[h] > 0)
            printf("  H=%d: %d\n", h, h_hist[h]);
    printf("  H>=200: %lld\n", h_ge200);

    return 0;
}
