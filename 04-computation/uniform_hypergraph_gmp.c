/*
 * Compute k-uniform hypergraphs on n unlabeled nodes.
 *
 * a(n) = (1/n!) * sum_{lambda} permcount(lambda) * 2^{orbits_k(lambda)}
 *
 * where orbits_k = number of orbits of a permutation with cycle type lambda
 * acting on k-element subsets of [n].
 *
 * orbits_k = (1/L) * sum_{t=0}^{L-1} fix_k(sigma^t)
 *
 * fix_k(sigma^t) = number of k-subsets that are unions of orbits of sigma^t
 *
 * For k=2: A000088 (simple graphs), 29 terms on OEIS
 * For k=3: A000665 (3-uniform hypergraphs), 29 terms on OEIS
 * For k=4: A051240 (4-uniform hypergraphs)
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o uniform_hypergraph_gmp uniform_hypergraph_gmp.c -lgmp
 *
 * Usage: ./uniform_hypergraph_gmp <k> <max_n>
 *
 * Author: opus-2026-03-09-S50
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

static int gcd_func(int a, int b) {
    while (b) { int t = b; b = a % b; a = t; }
    return a;
}

static long long lcm_func(long long a, long long b) {
    return a / gcd_func((int)a, (int)b) * b;
}

static int K;  /* uniformity */
static int N;
static mpz_t result, fact_n;

/*
 * Compute the number of k-subsets fixed by sigma^t.
 * For sigma with cycles c_1,...,c_m:
 * sigma^t has orbits: cycle i contributes gcd(c_i,t) orbits of size c_i/gcd(c_i,t).
 * A k-subset is fixed iff it's a union of orbits of sigma^t.
 * Count: number of multisubsets of orbit sizes summing to k.
 *
 * We enumerate all orbit sizes and their multiplicities,
 * then count subsets summing to exactly k using DP.
 */
static long long count_fixed_ksubsets(int *cycles, int m, long long t, int k) {
    /* Build list of orbit sizes and their counts */
    /* For each cycle of length c_i, sigma^t gives gcd(c_i,t) orbits of size c_i/gcd(c_i,t) */
    /* Collect (size, count) pairs, merging same sizes */

    int sizes[128], counts[128];
    int num_types = 0;

    for (int i = 0; i < m; i++) {
        int g = gcd_func(cycles[i], (int)(t % cycles[i]));
        if (t == 0) g = cycles[i];
        int sz = cycles[i] / g;
        int cnt = g;

        /* Merge with existing */
        int found = 0;
        for (int j = 0; j < num_types; j++) {
            if (sizes[j] == sz) {
                counts[j] += cnt;
                found = 1;
                break;
            }
        }
        if (!found) {
            sizes[num_types] = sz;
            counts[num_types] = cnt;
            num_types++;
        }
    }

    /* DP: count number of ways to select some orbits summing to exactly k */
    /* dp[s] = number of ways to pick orbits summing to s */
    /* Since k is small (typically 2-5), use simple DP */
    long long dp[64];
    memset(dp, 0, sizeof(dp));
    dp[0] = 1;

    for (int i = 0; i < num_types; i++) {
        int sz = sizes[i];
        int cnt = counts[i];
        if (sz > k) continue;

        /* Process 'cnt' identical items of size 'sz' */
        /* Use bounded knapsack: for each item type, update dp */
        /* Since k is small, we can use the standard approach */
        for (int iter = 0; iter < cnt; iter++) {
            /* Add one more orbit of size sz */
            for (int s = k; s >= sz; s--) {
                dp[s] += dp[s - sz];
            }
        }
    }

    return dp[k];
}

static void process_partition(int *pk, int *pm, int depth) {
    int num_cycles = 0;
    for (int i = 0; i < depth; i++) num_cycles += pm[i];
    int cycles[512];
    int idx = 0;
    for (int i = 0; i < depth; i++)
        for (int j = 0; j < pm[i]; j++)
            cycles[idx++] = pk[i];

    /* L = lcm of cycle lengths */
    long long L = 1;
    for (int i = 0; i < num_cycles; i++) {
        L = lcm_func(L, cycles[i]);
        if (L > 100000000LL) {
            fprintf(stderr, "Warning: L=%lld too large for partition, skipping\n", L);
            return;
        }
    }

    /* orbits_k = (1/L) * sum_{t=0}^{L-1} fix_k(sigma^t) */
    mpz_t orbits_sum;
    mpz_init(orbits_sum);
    mpz_set_ui(orbits_sum, 0);

    for (long long t = 0; t < L; t++) {
        long long fix = count_fixed_ksubsets(cycles, num_cycles, t, K);
        if (fix > 0)
            mpz_add_ui(orbits_sum, orbits_sum, (unsigned long)fix);
    }

    mpz_divexact_ui(orbits_sum, orbits_sum, (unsigned long)L);

    /* z_lambda = prod r^m * m! */
    mpz_t z_lambda;
    mpz_init(z_lambda);
    mpz_set_ui(z_lambda, 1);
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < pm[i]; j++)
            mpz_mul_ui(z_lambda, z_lambda, pk[i]);
        for (int j = 2; j <= pm[i]; j++)
            mpz_mul_ui(z_lambda, z_lambda, j);
    }

    /* contrib = (n! / z_lambda) * 2^orbits_k */
    mpz_t contrib;
    mpz_init(contrib);
    mpz_divexact(contrib, fact_n, z_lambda);

    if (mpz_fits_ulong_p(orbits_sum)) {
        unsigned long orb_val = mpz_get_ui(orbits_sum);
        mpz_t pow2;
        mpz_init(pow2);
        mpz_set_ui(pow2, 1);
        mpz_mul_2exp(pow2, pow2, orb_val);
        mpz_mul(contrib, contrib, pow2);
        mpz_clear(pow2);
    } else {
        fprintf(stderr, "Warning: orbits too large for shift\n");
    }

    mpz_add(result, result, contrib);

    mpz_clear(contrib);
    mpz_clear(z_lambda);
    mpz_clear(orbits_sum);
}

static void enumerate(int remaining, int max_part, int *pk, int *pm, int depth) {
    if (remaining == 0) {
        process_partition(pk, pm, depth);
        return;
    }
    for (int p = (max_part < remaining ? max_part : remaining); p >= 1; p--) {
        int max_m = remaining / p;
        for (int m = 1; m <= max_m; m++) {
            pk[depth] = p;
            pm[depth] = m;
            enumerate(remaining - m * p, p - 1, pk, pm, depth + 1);
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <k> <max_n>\n", argv[0]);
        fprintf(stderr, "  k=2: simple graphs (A000088)\n");
        fprintf(stderr, "  k=3: 3-uniform hypergraphs (A000665)\n");
        fprintf(stderr, "  k=4: 4-uniform hypergraphs (A051240)\n");
        return 1;
    }

    K = atoi(argv[1]);
    int max_n = atoi(argv[2]);

    fprintf(stderr, "Computing %d-uniform hypergraphs on n unlabeled nodes, n=0..%d\n\n", K, max_n);

    printf("0 1\n");
    fflush(stdout);

    for (N = 1; N <= max_n; N++) {
        struct timespec t_start, t_end;
        clock_gettime(CLOCK_MONOTONIC, &t_start);

        mpz_init(result);
        mpz_init(fact_n);
        mpz_fac_ui(fact_n, N);
        mpz_set_ui(result, 0);

        int pk[64], pm[64];
        enumerate(N, N, pk, pm, 0);
        mpz_divexact(result, result, fact_n);

        clock_gettime(CLOCK_MONOTONIC, &t_end);
        double dt = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

        gmp_printf("%d %Zd\n", N, result);
        fprintf(stderr, "n=%d: %zu digits, %.3fs\n", N, mpz_sizeinbase(result, 10), dt);
        fflush(stdout);
        fflush(stderr);

        mpz_clear(result);
        mpz_clear(fact_n);
        if (dt > 300) { fprintf(stderr, "Stopping\n"); break; }
    }
    return 0;
}
