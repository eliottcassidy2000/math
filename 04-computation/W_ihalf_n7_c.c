/*
 * Fast computation of H(T) and W(i/2) for all 2^21 tournaments on n=7.
 *
 * W(i/2) = P(i)/64 where P(i) = Re[sum_P prod((i+1) or (i-1))]
 * = (N_1 - N_3 + N_5) where N_f = #{perms with f forward edges}
 * (using the coefficient table from the Python analysis)
 *
 * Actually P(i) = 8*(N_1 - N_3 + N_5), so W(i/2) = (N_1 - N_3 + N_5)/8.
 *
 * We compute P(i) using bitmask DP with complex multiplication:
 *   forward edge: multiply by (1+i), i.e. (r,im) -> (r-im, r+im)
 *   backward edge: multiply by (-1+i), i.e. (r,im) -> (-r-im, r-im)
 *
 * And H using standard bitmask DP.
 *
 * Compile: gcc -O2 -o W_ihalf_n7_c W_ihalf_n7_c.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define N 7
#define M 21  /* N*(N-1)/2 */
#define FULL ((1 << N) - 1)  /* 127 */
#define TOTAL (1 << M)       /* 2097152 */

/* Precomputed pair-to-bit mapping */
int pair_bit[N][N];  /* pair_bit[i][j] for i<j = bit index */

/* DP arrays - reused per tournament */
/* For H computation */
long long dp_H[1 << N][N];

/* For P(i) computation (real and imag parts) */
long long dp_r[1 << N][N];
long long dp_im[1 << N][N];

/* Results */
int H_result[TOTAL];
long long Pi_result[TOTAL];  /* P(i) = 8*(N1-N3+N5), always real */
int c3_result[TOTAL];

void init_pairs() {
    int k = 0;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++)
            pair_bit[i][j] = k++;
}

/* Get T[v][u] from bits */
static inline int get_edge(int bits, int v, int u) {
    if (v < u) return (bits >> pair_bit[v][u]) & 1;
    else return 1 - ((bits >> pair_bit[u][v]) & 1);
}

void compute_tournament(int bits) {
    /* Clear DP arrays */
    memset(dp_H, 0, sizeof(dp_H));
    memset(dp_r, 0, sizeof(dp_r));
    memset(dp_im, 0, sizeof(dp_im));

    /* Initialize single-vertex paths */
    for (int v = 0; v < N; v++) {
        dp_H[1 << v][v] = 1;
        dp_r[1 << v][v] = 1;
        dp_im[1 << v][v] = 0;
    }

    /* DP transitions */
    for (int mask = 1; mask < (1 << N); mask++) {
        for (int v = 0; v < N; v++) {
            if (!(mask & (1 << v))) continue;

            long long h_cnt = dp_H[mask][v];
            long long r_val = dp_r[mask][v];
            long long im_val = dp_im[mask][v];

            if (h_cnt == 0 && r_val == 0 && im_val == 0) continue;

            for (int u = 0; u < N; u++) {
                if (mask & (1 << u)) continue;

                int new_mask = mask | (1 << u);
                int fwd = get_edge(bits, v, u);

                /* H DP */
                if (fwd && h_cnt > 0) {
                    dp_H[new_mask][u] += h_cnt;
                }

                /* P(i) DP */
                if (fwd) {
                    /* multiply by (1+i): (r,im) -> (r-im, r+im) */
                    dp_r[new_mask][u] += r_val - im_val;
                    dp_im[new_mask][u] += r_val + im_val;
                } else {
                    /* multiply by (-1+i): (r,im) -> (-r-im, r-im) */
                    dp_r[new_mask][u] += -r_val - im_val;
                    dp_im[new_mask][u] += r_val - im_val;
                }
            }
        }
    }

    /* Collect results */
    long long H = 0, Pi_real = 0, Pi_imag = 0;
    for (int v = 0; v < N; v++) {
        H += dp_H[FULL][v];
        Pi_real += dp_r[FULL][v];
        Pi_imag += dp_im[FULL][v];
    }

    H_result[bits] = (int)H;
    Pi_result[bits] = Pi_real;

    /* Compute c3 via out-degree formula */
    int outdeg[N] = {0};
    for (int v = 0; v < N; v++)
        for (int u = 0; u < N; u++)
            if (v != u) outdeg[v] += get_edge(bits, v, u);

    int sum_od2 = 0;
    for (int v = 0; v < N; v++) sum_od2 += outdeg[v] * outdeg[v];
    c3_result[bits] = 35 - (sum_od2 - 21) / 2;
}

int main() {
    init_pairs();

    fprintf(stderr, "Computing H and P(i) for all %d tournaments on n=%d...\n", TOTAL, N);

    clock_t t0 = clock();

    for (int bits = 0; bits < TOTAL; bits++) {
        if (bits % (TOTAL/20) == 0 && bits > 0) {
            double elapsed = (double)(clock() - t0) / CLOCKS_PER_SEC;
            double rate = bits / elapsed;
            double eta = (TOTAL - bits) / rate;
            fprintf(stderr, "  %d/%d (%.0f%%) [%.1fs, ETA %.1fs]\n",
                    bits, TOTAL, 100.0*bits/TOTAL, elapsed, eta);
        }
        compute_tournament(bits);
    }

    double elapsed = (double)(clock() - t0) / CLOCKS_PER_SEC;
    fprintf(stderr, "Completed in %.1fs\n", elapsed);

    /* === ANALYSIS === */

    printf("W(i/2) vs H(T) analysis for all %d tournaments on n=%d\n", TOTAL, N);
    printf("W(i/2) = P(i)/64, P(i) = 8*(N1 - N3 + N5)\n");
    printf("======================================================================\n");

    /* Check P(i) divisibility */
    int div8_count = 0;
    for (int i = 0; i < TOTAL; i++)
        if (Pi_result[i] % 8 == 0) div8_count++;
    printf("\nP(i) divisible by 8: %d/%d\n", div8_count, TOTAL);

    /* W*8 = P(i)/8 ... but actually P(i) = 8*(N1-N3+N5), W = (N1-N3+N5)/8 = P(i)/64 */
    /* Let W8 = P(i) / 8 (this is N1 - N3 + N5, always integer) */

    /* Build H distribution */
    int H_min = 9999, H_max = 0;
    for (int i = 0; i < TOTAL; i++) {
        if (H_result[i] < H_min) H_min = H_result[i];
        if (H_result[i] > H_max) H_max = H_result[i];
    }

    printf("\nH range: [%d, %d]\n", H_min, H_max);

    /* Count by H */
    /* H is always odd, range roughly 1-189 */
    #define MAX_H 200
    int H_count[MAX_H] = {0};
    long long W8_min[MAX_H], W8_max[MAX_H];
    long long absW8_min[MAX_H], absW8_max[MAX_H];
    double absW8_sum[MAX_H] = {0};
    int c3_min_h[MAX_H], c3_max_h[MAX_H];

    for (int h = 0; h < MAX_H; h++) {
        W8_min[h] = 999999999LL;
        W8_max[h] = -999999999LL;
        absW8_min[h] = 999999999LL;
        absW8_max[h] = 0;
        c3_min_h[h] = 9999;
        c3_max_h[h] = -1;
    }

    int W_zero_count = 0;
    int W_pos_count = 0;
    int W_neg_count = 0;

    for (int i = 0; i < TOTAL; i++) {
        int h = H_result[i];
        long long w8 = Pi_result[i] / 8;  /* = N1 - N3 + N5 */
        long long aw8 = w8 < 0 ? -w8 : w8;
        int c = c3_result[i];

        H_count[h]++;
        if (w8 < W8_min[h]) W8_min[h] = w8;
        if (w8 > W8_max[h]) W8_max[h] = w8;
        if (aw8 < absW8_min[h]) absW8_min[h] = aw8;
        if (aw8 > absW8_max[h]) absW8_max[h] = aw8;
        absW8_sum[h] += aw8;
        if (c < c3_min_h[h]) c3_min_h[h] = c;
        if (c > c3_max_h[h]) c3_max_h[h] = c;

        if (w8 == 0) W_zero_count++;
        else if (w8 > 0) W_pos_count++;
        else W_neg_count++;
    }

    printf("\nW(i/2) sign: pos=%d (%.1f%%), zero=%d (%.1f%%), neg=%d (%.1f%%)\n",
           W_pos_count, 100.0*W_pos_count/TOTAL,
           W_zero_count, 100.0*W_zero_count/TOTAL,
           W_neg_count, 100.0*W_neg_count/TOTAL);

    /* Correlations */
    double sum_H = 0, sum_W = 0, sum_aW = 0, sum_c3 = 0;
    for (int i = 0; i < TOTAL; i++) {
        sum_H += H_result[i];
        long long w8 = Pi_result[i] / 8;
        sum_W += w8;
        sum_aW += (w8 < 0 ? -w8 : w8);
        sum_c3 += c3_result[i];
    }
    double mean_H = sum_H / TOTAL;
    double mean_W = sum_W / TOTAL;
    double mean_aW = sum_aW / TOTAL;
    double mean_c3 = sum_c3 / TOTAL;

    double cov_HW = 0, cov_HaW = 0, cov_Hc3 = 0, cov_c3W = 0, cov_c3aW = 0;
    double var_H = 0, var_W = 0, var_aW = 0, var_c3 = 0;
    for (int i = 0; i < TOTAL; i++) {
        double dH = H_result[i] - mean_H;
        long long w8 = Pi_result[i] / 8;
        double dW = w8 - mean_W;
        double daW = (w8 < 0 ? -w8 : w8) - mean_aW;
        double dc3 = c3_result[i] - mean_c3;

        var_H += dH * dH;
        var_W += dW * dW;
        var_aW += daW * daW;
        var_c3 += dc3 * dc3;
        cov_HW += dH * dW;
        cov_HaW += dH * daW;
        cov_Hc3 += dH * dc3;
        cov_c3W += dc3 * dW;
        cov_c3aW += dc3 * daW;
    }

    printf("\nCORRELATIONS:\n");
    printf("  Corr(H, W8)     = %.6f\n", cov_HW / sqrt(var_H * var_W));
    printf("  Corr(H, |W8|)   = %.6f\n", cov_HaW / sqrt(var_H * var_aW));
    printf("  Corr(H, c3)     = %.6f\n", cov_Hc3 / sqrt(var_H * var_c3));
    printf("  Corr(c3, W8)    = %.6f\n", cov_c3W / sqrt(var_c3 * var_W));
    printf("  Corr(c3, |W8|)  = %.6f\n", cov_c3aW / sqrt(var_c3 * var_aW));

    /* Joint distribution table */
    printf("\nJOINT DISTRIBUTION (H vs W8 = P(i)/8 = N1-N3+N5):\n");
    printf("  W(i/2) = W8/8\n\n");
    printf("  %5s  %7s  %8s  %8s  %9s  %9s  %10s  %6s\n",
           "H", "count", "W8 min", "W8 max", "|W8| min", "|W8| max", "mean|W8|", "c3");

    for (int h = 0; h < MAX_H; h++) {
        if (H_count[h] == 0) continue;
        printf("  %5d  %7d  %8lld  %8lld  %9lld  %9lld  %10.1f  [%d,%d]\n",
               h, H_count[h], W8_min[h], W8_max[h],
               absW8_min[h], absW8_max[h],
               absW8_sum[h] / H_count[h],
               c3_min_h[h], c3_max_h[h]);
    }

    /* Monotonicity of mean |W8| */
    printf("\nMONOTONICITY OF MEAN |W8| vs H:\n");
    double prev_mean_aW = 1e18;
    int mono_violations = 0;
    for (int h = 0; h < MAX_H; h++) {
        if (H_count[h] == 0) continue;
        double m = absW8_sum[h] / H_count[h];
        if (m > prev_mean_aW) {
            mono_violations++;
            printf("  VIOLATION at H=%d: mean|W8|=%.1f > prev=%.1f\n", h, m, prev_mean_aW);
        }
        prev_mean_aW = m;
    }
    printf("  Total mean monotonicity violations: %d\n", mono_violations);

    /* W=0 analysis */
    printf("\nW(i/2) = 0 ANALYSIS:\n");
    printf("  Total: %d tournaments\n", W_zero_count);
    printf("  H distribution when W=0:\n");
    /* Recount */
    int W0_by_H[MAX_H] = {0};
    for (int i = 0; i < TOTAL; i++) {
        if (Pi_result[i] == 0) W0_by_H[H_result[i]]++;
    }
    for (int h = 0; h < MAX_H; h++) {
        if (W0_by_H[h] > 0) printf("    H=%d: %d\n", h, W0_by_H[h]);
    }

    /* Max H analysis */
    printf("\nMAX H = %d:\n", H_max);
    long long W_at_max = -999999;
    int max_count = 0;
    for (int i = 0; i < TOTAL; i++) {
        if (H_result[i] == H_max) {
            max_count++;
            W_at_max = Pi_result[i] / 8;
        }
    }
    printf("  Count: %d, W8 = %lld, W(i/2) = %.4f\n", max_count, W_at_max, W_at_max / 8.0);

    /* W^2 analysis */
    printf("\n|W(i/2)|^2 ANALYSIS:\n");
    /* Check if W8^2 is always divisible by 64 (so W^2 is integer) */
    int w2_div64 = 0;
    for (int i = 0; i < TOTAL; i++) {
        long long w8 = Pi_result[i] / 8;
        if ((w8 * w8) % 64 == 0) w2_div64++;
    }
    printf("  W8^2 always div by 64 (W^2 integer): %d/%d\n", w2_div64, TOTAL);

    /* Distinct W8 values */
    #define MAX_W8 10000
    /* Use sorting approach */
    long long *w8_vals = (long long *)malloc(TOTAL * sizeof(long long));
    for (int i = 0; i < TOTAL; i++) w8_vals[i] = Pi_result[i] / 8;

    /* Count distinct */
    /* Simple approach: just count unique W8^2 values up to some bound */
    int distinct_absw8 = 0;
    long long prev_val = -1;
    /* Sort absW8 */
    long long *aw8 = (long long *)malloc(TOTAL * sizeof(long long));
    for (int i = 0; i < TOTAL; i++) {
        long long w = w8_vals[i];
        aw8[i] = w < 0 ? -w : w;
    }

    /* Quick count of distinct |W8| values using brute force on small range */
    long long aw8_max = 0;
    for (int i = 0; i < TOTAL; i++)
        if (aw8[i] > aw8_max) aw8_max = aw8[i];

    printf("  |W8| max = %lld\n", aw8_max);

    /* Count distinct |W8| */
    char *seen = (char *)calloc(aw8_max + 1, 1);
    int n_distinct = 0;
    for (int i = 0; i < TOTAL; i++) {
        if (!seen[aw8[i]]) {
            seen[aw8[i]] = 1;
            n_distinct++;
        }
    }
    printf("  Distinct |W8| values: %d\n", n_distinct);

    free(seen);
    free(w8_vals);
    free(aw8);

    return 0;
}
