/*
 * h21_5cycle_forcing.c — For cycle-rich n=8 tournaments, compute the
 * relationship between t3 (3-cycle vertex sets) and t5 (5-cycle vertex sets).
 *
 * Hypothesis: cycle-rich implies t3 + t5 > 10 (i.e., alpha_1 > 10
 * counting just 3+5 cycles), which would prove H > 21 directly.
 *
 * Author: opus-2026-03-07-S42
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 8
#define NBITS (N*(N-1)/2)
#define TOTAL (1U << NBITS)

static inline void compute_out_masks(unsigned int T, unsigned char out[N]) {
    memset(out, 0, N);
    int pos = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            if ((T >> pos) & 1)
                out[i] |= (1 << j);
            else
                out[j] |= (1 << i);
            pos++;
        }
    }
}

static inline int has_source_sink(const unsigned char out[N]) {
    for (int i = 0; i < N; i++) {
        int s = __builtin_popcount(out[i]);
        if (s == 0 || s == N-1) return 1;
    }
    return 0;
}

static inline int count_3cycles(const unsigned char out[N]) {
    int sum_c2 = 0;
    for (int i = 0; i < N; i++) {
        int s = __builtin_popcount(out[i]);
        sum_c2 += s * (s - 1) / 2;
    }
    return 56 - sum_c2;
}

static inline int all_in_3cycle(const unsigned char out[N]) {
    unsigned char all = (1 << N) - 1;
    for (int v = 0; v < N; v++) {
        int found = 0;
        unsigned char ov = out[v];
        unsigned char iv = (~out[v]) & all & ~(1 << v);
        unsigned char ov_copy = ov;
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

/* Count 3-cycle vertex sets */
static int count_3cycle_sets(const unsigned char out[N]) {
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

/* Count 5-cycle vertex sets */
static int count_5cycle_sets(const unsigned char out[N]) {
    int count = 0;
    for (int a = 0; a < N; a++)
      for (int b = a+1; b < N; b++)
        for (int c = b+1; c < N; c++)
          for (int d = c+1; d < N; d++)
            for (int e = d+1; e < N; e++) {
                int v[5] = {a,b,c,d,e};
                int dp[32][5];
                memset(dp, 0, sizeof(dp));
                dp[1][0] = 1;
                for (int mask = 1; mask < 32; mask++)
                    for (int i = 0; i < 5; i++) {
                        if (!(mask & (1<<i)) || dp[mask][i] == 0) continue;
                        for (int j = 0; j < 5; j++) {
                            if (mask & (1<<j)) continue;
                            if ((out[v[i]] >> v[j]) & 1)
                                dp[mask|(1<<j)][j] += dp[mask][i];
                        }
                    }
                int cyc = 0;
                for (int j = 1; j < 5; j++)
                    if (dp[31][j] && ((out[v[j]] >> v[0]) & 1))
                        cyc += dp[31][j];
                if (cyc > 0) count++;
            }
    return count;
}

int main() {
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    long long cycle_rich = 0;
    int min_t3_plus_t5 = 999;
    /* Track (t3, t5) pairs */
    int pair_count[20][30];  /* [t3][t5] */
    memset(pair_count, 0, sizeof(pair_count));
    unsigned char out[N];

    for (unsigned int T = 0; T < TOTAL; T++) {
        compute_out_masks(T, out);
        if (has_source_sink(out)) continue;
        int t3 = count_3cycles(out);
        if (t3 > 10) continue;
        if (!all_in_3cycle(out)) continue;

        cycle_rich++;

        int t3s = count_3cycle_sets(out);
        int t5s = count_5cycle_sets(out);
        int total = t3s + t5s;

        if (total < min_t3_plus_t5) {
            min_t3_plus_t5 = total;
            printf("NEW MIN: t3s=%d t5s=%d total=%d T=%u\n", t3s, t5s, total, T);
        }

        if (t3s < 20 && t5s < 30)
            pair_count[t3s][t5s]++;

        if ((T & 0xFFFFFF) == 0xFFFFFF) {
            fprintf(stderr, "Progress: %u/%u (%.1f%%). cr=%lld min=%d\n",
                    T + 1, TOTAL, 100.0 * (T + 1) / TOTAL, cycle_rich, min_t3_plus_t5);
        }
    }

    printf("\n=== EXHAUSTIVE n=%d RESULTS ===\n", N);
    printf("Cycle-rich: %lld\n", cycle_rich);
    printf("Min (t3_sets + t5_sets): %d\n", min_t3_plus_t5);
    printf("\n(t3_sets, t5_sets) distribution:\n");
    for (int t3 = 0; t3 < 20; t3++)
        for (int t5 = 0; t5 < 30; t5++)
            if (pair_count[t3][t5] > 0)
                printf("  t3=%d t5=%d: %d\n", t3, t5, pair_count[t3][t5]);

    return 0;
}
