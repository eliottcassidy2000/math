/*
 * h21_min_weight_3only.c — For cycle-rich n=8 tournaments,
 * compute min of t3_sets + 2 * disjoint_3cycle_pairs.
 *
 * Since alpha_1 >= t3_sets (3-cycles are a subset of all odd cycles),
 * and disjoint 3-cycle pairs are a subset of all disjoint pairs,
 * the true OCF weight is >= this quantity.
 *
 * If min > 10, we're done.
 * If min <= 10, we need to add 5-cycles and 7-cycles.
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

typedef struct { unsigned char a, b, c; } Triple;

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
        int t3 = count_3cycles(out);
        if (t3 > 10) continue;
        if (!all_in_3cycle(out)) continue;

        cycle_rich++;

        /* Find all 3-cycle vertex sets */
        Triple cyc[60];
        int nc = 0;
        for (int a = 0; a < N; a++)
            for (int b = a+1; b < N; b++)
                for (int c = b+1; c < N; c++) {
                    int ab = (out[a] >> b) & 1;
                    int bc = (out[b] >> c) & 1;
                    int ca = (out[c] >> a) & 1;
                    if ((ab && bc && ca) || (!ab && !bc && !ca)) {
                        if (nc < 60) {
                            cyc[nc].a = a; cyc[nc].b = b; cyc[nc].c = c;
                            nc++;
                        }
                    }
                }

        /* Count disjoint pairs among 3-cycles */
        int i2 = 0;
        for (int i = 0; i < nc; i++)
            for (int j = i+1; j < nc; j++) {
                unsigned char m1 = (1<<cyc[i].a)|(1<<cyc[i].b)|(1<<cyc[i].c);
                unsigned char m2 = (1<<cyc[j].a)|(1<<cyc[j].b)|(1<<cyc[j].c);
                if ((m1 & m2) == 0) i2++;
            }

        int w = nc + 2 * i2;
        if (w < min_weight) {
            min_weight = w;
            printf("NEW MIN: w=%d (t3=%d, i2_3only=%d) T=%u\n", w, nc, i2, T);
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
    printf("Min weight (t3 + 2*i2_3only): %d\n", min_weight);
    printf("\nWeight distribution (first 30):\n");
    for (int w = 0; w < 30; w++)
        if (weight_hist[w] > 0)
            printf("  w=%d: %d\n", w, weight_hist[w]);

    return 0;
}
