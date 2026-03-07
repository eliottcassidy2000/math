/*
 * h21_min_ocf_weight.c — Compute min of alpha_1 + 2*alpha_2 for cycle-rich
 * tournaments at n=8 (exhaustive) and n=9 (sampling).
 *
 * For H=21: need alpha_1 + 2*alpha_2 + 4*alpha_3 + ... = 10.
 * So alpha_1 + 2*alpha_2 <= 10.
 * If min(alpha_1 + 2*alpha_2) > 10 for cycle-rich, H=21 is impossible.
 *
 * "Cycle-rich" = no source/sink, every vertex in 3-cycle.
 *
 * Author: opus-2026-03-07-S42
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N 8
#define NBITS (N*(N-1)/2)
#define TOTAL (1U << NBITS)

static inline void compute_out_masks(unsigned int T, unsigned char out[N]) {
    memset(out, 0, N);
    int pos = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            if ((T >> pos) & 1) {
                out[i] |= (1 << j);
            } else {
                out[j] |= (1 << i);
            }
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
    return 28 - sum_c2;  /* C(8,3) = 56, wait no C(8,3)=56 */
}

/* Actually C(8,3) = 56, not 28 */
static inline int count_3cycles_correct(const unsigned char out[N]) {
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

/* Count 3-cycle vertex SETS */
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

/* Check if two triples are disjoint */
static inline int triples_disjoint(int a1, int b1, int c1, int a2, int b2, int c2) {
    if (a1==a2 || a1==b2 || a1==c2) return 0;
    if (b1==a2 || b1==b2 || b1==c2) return 0;
    if (c1==a2 || c1==b2 || c1==c2) return 0;
    return 1;
}

/* Count alpha_2 = number of PAIRWISE-DISJOINT 3-cycle pairs
 * Actually: alpha_2 = max independent set of size 2 in Omega.
 * For decomposition purposes: i_2 = number of independent pairs.
 * Let's just count the number of disjoint pairs among 3-cycle vertex sets. */
typedef struct { unsigned char a, b, c; } Triple;

static int compute_alpha1_and_i2(const unsigned char out[N], int *alpha1_out, int *i2_out) {
    /* Find all 3-cycle vertex sets */
    Triple cyc[100];
    int nc = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            for (int c = b+1; c < N; c++) {
                int ab = (out[a] >> b) & 1;
                int bc = (out[b] >> c) & 1;
                int ca = (out[c] >> a) & 1;
                if ((ab && bc && ca) || (!ab && !bc && !ca)) {
                    if (nc < 100) {
                        cyc[nc].a = a; cyc[nc].b = b; cyc[nc].c = c;
                        nc++;
                    }
                }
            }

    /* Also count 5-cycle vertex sets */
    int t5_sets = 0;
    for (int a = 0; a < N; a++)
      for (int b = a+1; b < N; b++)
        for (int c = b+1; c < N; c++)
          for (int d = c+1; d < N; d++)
            for (int e = d+1; e < N; e++) {
                int v[5] = {a,b,c,d,e};
                /* Held-Karp for 5-vertex Hamiltonian cycle */
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
                int cyc5 = 0;
                for (int j = 1; j < 5; j++)
                    if (dp[31][j] && ((out[v[j]] >> v[0]) & 1))
                        cyc5 += dp[31][j];
                if (cyc5 > 0) t5_sets++;
            }

    /* Also count 7-cycle vertex sets */
    int t7_sets = 0;
    for (int a = 0; a < N; a++)
      for (int b = a+1; b < N; b++)
        for (int c = b+1; c < N; c++)
          for (int d = c+1; d < N; d++)
            for (int e = d+1; e < N; e++)
              for (int f = e+1; f < N; f++)
                for (int g = f+1; g < N; g++) {
                    int v[7] = {a,b,c,d,e,f,g};
                    int dp[128][7];
                    memset(dp, 0, sizeof(dp));
                    dp[1][0] = 1;
                    for (int mask = 1; mask < 128; mask++)
                        for (int i = 0; i < 7; i++) {
                            if (!(mask & (1<<i)) || dp[mask][i] == 0) continue;
                            for (int j = 0; j < 7; j++) {
                                if (mask & (1<<j)) continue;
                                if ((out[v[i]] >> v[j]) & 1)
                                    dp[mask|(1<<j)][j] += dp[mask][i];
                            }
                        }
                    int cyc7 = 0;
                    for (int j = 1; j < 7; j++)
                        if (dp[127][j] && ((out[v[j]] >> v[0]) & 1))
                            cyc7 += dp[127][j];
                    if (cyc7 > 0) t7_sets++;
                }

    int alpha1 = nc + t5_sets + t7_sets;

    /* Count disjoint pairs among ALL cycles (3+5+7) */
    /* For simplicity, just count among 3-cycles for now */
    int i2 = 0;
    for (int i = 0; i < nc; i++)
        for (int j = i+1; j < nc; j++)
            if (triples_disjoint(cyc[i].a, cyc[i].b, cyc[i].c,
                                 cyc[j].a, cyc[j].b, cyc[j].c))
                i2++;

    *alpha1_out = alpha1;
    *i2_out = i2;
    return alpha1 + 2 * i2;
}

int main() {
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    int min_weight = 999;
    long long cycle_rich = 0;
    int weight_hist[100] = {0};
    unsigned char out[N];

    for (unsigned int T = 0; T < TOTAL; T++) {
        compute_out_masks(T, out);
        if (has_source_sink(out)) continue;
        int t3 = count_3cycles_correct(out);
        if (t3 > 10) continue;
        if (!all_in_3cycle(out)) continue;

        cycle_rich++;

        int alpha1, i2;
        int w = compute_alpha1_and_i2(out, &alpha1, &i2);
        if (w < min_weight) {
            min_weight = w;
            printf("NEW MIN: w=%d (a1=%d, i2=%d) t3=%d T=%u\n", w, alpha1, i2, t3, T);
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    if (i == j) printf(".");
                    else printf("%d", (out[i] >> j) & 1);
                }
                printf("\n");
            }
            printf("\n");
        }
        if (w < 100) weight_hist[w]++;

        if ((T & 0xFFFFFF) == 0xFFFFFF) {
            fprintf(stderr, "Progress: %u/%u (%.1f%%). cycle_rich=%lld min_w=%d\n",
                    T + 1, TOTAL, 100.0 * (T + 1) / TOTAL, cycle_rich, min_weight);
        }
    }

    printf("\n=== EXHAUSTIVE n=%d RESULTS ===\n", N);
    printf("Cycle-rich tournaments: %lld\n", cycle_rich);
    printf("Min OCF weight (alpha1 + 2*i2): %d\n", min_weight);
    printf("\nWeight distribution (first 40):\n");
    for (int w = 0; w < 40; w++)
        if (weight_hist[w] > 0)
            printf("  w=%d: %d\n", w, weight_hist[w]);

    return 0;
}
