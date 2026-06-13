/*
 * h_spectrum_n8_sampled.c — kind-pasteur-2026-03-21-S17c
 *
 * H-spectrum at n=8,9,10,11 via massive sampling.
 * n=8: exhaustive is possible (2^28 = 268M, ~55min) but let's also sample
 * n=9,10,11: sample millions of random tournaments
 *
 * For gap detection at n=8 we need exhaustive. For n>=9 we track:
 * - Complete frequency histogram
 * - Which odd values are NEVER seen (candidate gaps)
 * - Score-class conditional distribution (for OCR connection)
 *
 * Compile: gcc -O3 -o h_spectrum_sampled h_spectrum_n8_sampled.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define MAX_N 12
#define MAX_H 200000

typedef unsigned long long u64;
typedef long long i64;

static int adj[MAX_N][MAX_N];

static int ham_paths(int n) {
    int dp[1 << MAX_N][MAX_N];
    int full = (1 << n) - 1;
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (int mask = 1; mask < (1 << n); mask++)
        for (int v = 0; v < n; v++) {
            if (!(mask & (1 << v))) continue;
            int c = dp[mask][v];
            if (c == 0) continue;
            for (int w = 0; w < n; w++) {
                if (mask & (1 << w)) continue;
                if (adj[v][w]) dp[mask | (1 << w)][w] += c;
            }
        }
    int total = 0;
    for (int v = 0; v < n; v++) total += dp[full][v];
    return total;
}

/* xorshift128+ RNG */
static u64 rng_s[2] = {0x12345678deadbeefULL, 0xfedcba9876543210ULL};
static u64 xorshift128plus(void) {
    u64 s1 = rng_s[0], s0 = rng_s[1];
    rng_s[0] = s0;
    s1 ^= s1 << 23;
    rng_s[1] = s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26);
    return rng_s[1] + s0;
}

/* Dynamic histogram using hash map for large H values */
#define HASH_SIZE 1000003
typedef struct entry { int key; int count; struct entry *next; } entry;
static entry *htable[HASH_SIZE];
static entry pool[5000000];
static int pool_idx = 0;

static void hist_add(int h) {
    int idx = ((unsigned)h) % HASH_SIZE;
    for (entry *e = htable[idx]; e; e = e->next) {
        if (e->key == h) { e->count++; return; }
    }
    entry *e = &pool[pool_idx++];
    e->key = h; e->count = 1; e->next = htable[idx];
    htable[idx] = e;
}

static int hist_get(int h) {
    int idx = ((unsigned)h) % HASH_SIZE;
    for (entry *e = htable[idx]; e; e = e->next)
        if (e->key == h) return e->count;
    return 0;
}

/* Collect all keys */
static int all_keys[5000000];
static int n_keys = 0;
static void hist_collect(void) {
    n_keys = 0;
    for (int i = 0; i < HASH_SIZE; i++)
        for (entry *e = htable[i]; e; e = e->next)
            all_keys[n_keys++] = e->key;
}

static int cmp_int(const void *a, const void *b) { return *(int*)a - *(int*)b; }

int main(int argc, char **argv) {
    int n = 8;
    i64 samples = 0;  /* 0 = exhaustive */
    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) samples = atoll(argv[2]);

    int m = n * (n-1) / 2;
    u64 total_possible = 1ULL << m;

    int exhaustive = (samples == 0 && m <= 28);
    i64 total_iter = exhaustive ? (i64)total_possible : samples;

    printf("H-SPECTRUM: n=%d, m=%d, %s (%lld iterations)\n",
           n, m, exhaustive ? "EXHAUSTIVE" : "SAMPLED", (long long)total_iter);
    fflush(stdout);

    int arc_i[100], arc_j[100];
    int pos = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) { arc_i[pos] = i; arc_j[pos] = j; pos++; }

    memset(htable, 0, sizeof(htable));
    pool_idx = 0;

    int max_h = 0, min_h = 999999;
    double sum_H = 0, sum_H2 = 0;
    i64 count = 0;

    /* Also track score S2 for OCR */
    double sum_S2 = 0, sum_S22 = 0, sum_HS2 = 0;

    time_t t0 = time(NULL);

    for (i64 iter = 0; iter < total_iter; iter++) {
        u64 bits;
        if (exhaustive) bits = (u64)iter;
        else bits = xorshift128plus();

        memset(adj, 0, sizeof(int)*MAX_N*MAX_N);
        for (int p = 0; p < m; p++) {
            if (bits & (1ULL << p)) adj[arc_i[p]][arc_j[p]] = 1;
            else adj[arc_j[p]][arc_i[p]] = 1;
        }

        int H = ham_paths(n);
        hist_add(H);
        count++;

        if (H > max_h) max_h = H;
        if (H < min_h) min_h = H;
        sum_H += H;
        sum_H2 += (double)H * H;

        /* S2 (integer form: 4*S2 = sum(2s-(n-1))^2) */
        int S2_4 = 0;
        for (int i = 0; i < n; i++) {
            int s = 0;
            for (int j = 0; j < n; j++) s += adj[i][j];
            int ds = 2*s - (n-1);
            S2_4 += ds * ds;
        }
        double S2 = S2_4 / 4.0;
        sum_S2 += S2;
        sum_S22 += S2 * S2;
        sum_HS2 += H * S2;

        if (count % 10000000 == 0) {
            time_t now = time(NULL);
            fprintf(stderr, "  %lld/%lld (%.1f%%) [%lds]\r",
                    (long long)count, (long long)total_iter,
                    100.0*count/total_iter, (long)(now-t0));
        }
    }

    time_t t1 = time(NULL);
    printf("Done in %lds. count=%lld\n\n", (long)(t1-t0), (long long)count);

    /* Collect and sort histogram */
    hist_collect();
    qsort(all_keys, n_keys, sizeof(int), cmp_int);

    printf("n=%d: %d distinct H values, range [%d, %d]\n", n, n_keys, min_h, max_h);

    double mean_H = sum_H / count;
    double var_H = sum_H2/count - mean_H*mean_H;
    printf("Mean H = %.4f (expected %.4f)\n", mean_H, tgamma(n+1)/(1<<(n-1)));
    printf("Std H = %.4f\n", sqrt(var_H));

    /* OCR */
    double mean_S2 = sum_S2/count;
    double var_S2 = sum_S22/count - mean_S2*mean_S2;
    double cov = sum_HS2/count - mean_H*mean_S2;
    double r2 = (cov*cov)/(var_S2*var_H);
    printf("OCR R^2 = %.8f\n\n", r2);

    /* Print frequency table */
    printf("H | count%s | frac\n", exhaustive ? "" : " (sampled)");
    int n_gaps = 0;
    int prev_h = -1;
    for (int i = 0; i < n_keys; i++) {
        int h = all_keys[i];
        int f = hist_get(h);

        /* Count gaps between prev and current */
        if (prev_h >= 0) {
            for (int g = prev_h + 2; g < h; g += 2)
                n_gaps++;
        }
        prev_h = h;

        /* Print (skip middle if too many) */
        if (n_keys <= 200 || i < 30 || i >= n_keys - 10 || f == hist_get(all_keys[0]) || h == max_h) {
            printf("%5d | %10d | %.8f\n", h, f, (double)f/count);
        } else if (i == 30) {
            printf("  ... (%d more values) ...\n", n_keys - 40);
        }
    }

    /* Gaps */
    printf("\nGaps (odd values in [1, %d] not achieved):\n", max_h);
    int gap_count = 0;
    prev_h = -1;
    for (int i = 0; i < n_keys; i++) {
        int h = all_keys[i];
        if (prev_h >= 0) {
            for (int g = prev_h + 2; g < h; g += 2) {
                if (gap_count < 50) printf("  %d", g);
                gap_count++;
            }
        }
        prev_h = h;
    }
    if (gap_count >= 50) printf(" ... (%d total gaps)", gap_count);
    printf("\n");

    /* Count permanent gap candidates (7 and 21) */
    printf("\nPermanent gap check: H=7 count=%d, H=21 count=%d\n",
           hist_get(7), hist_get(21));

    /* Frequency distribution shape */
    printf("\nSkewness and kurtosis:\n");
    double mu3 = 0, mu4 = 0;
    for (int i = 0; i < n_keys; i++) {
        int h = all_keys[i];
        int f = hist_get(h);
        double d = h - mean_H;
        mu3 += d*d*d * f;
        mu4 += d*d*d*d * f;
    }
    mu3 /= count;
    mu4 /= count;
    double sigma3 = var_H * sqrt(var_H);
    printf("  Skewness = %.6f\n", mu3/sigma3);
    printf("  Kurtosis = %.6f (Gaussian=3)\n", mu4/(var_H*var_H));

    fflush(stdout);
    return 0;
}
