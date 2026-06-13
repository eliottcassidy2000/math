/*
 * Compute A000568: number of tournaments on n unlabeled nodes.
 *
 * a(n) = (1/n!) * sum_{lambda partition of n into ODD parts}
 *            permcount(lambda) * 2^{t(lambda)}
 *
 * where t(lambda) = sum_{i<j} gcd(c_i,c_j) + (1/2)*sum_i(c_i - 1)
 *   c_1,...,c_m are cycle lengths (all odd)
 *
 * Only partitions into odd parts contribute (even cycles can't fix tournaments).
 *
 * OEIS: A000568 (tournaments), had 77 terms (n=0..76)
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o tournament_gmp tournament_gmp.c -lgmp
 *
 * Usage: ./tournament_gmp <max_n>
 *
 * Author: opus-2026-03-09-S50
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>

static int gcd_func(int a, int b) {
    while (b) { int t = b; b = a % b; a = t; }
    return a;
}

static int N;
static mpz_t result, fact_n;

static void process_partition(int *pk, int *pm, int depth) {
    int num_cycles = 0;
    for (int i = 0; i < depth; i++) num_cycles += pm[i];
    int cycles[512];
    int idx = 0;
    for (int i = 0; i < depth; i++)
        for (int j = 0; j < pm[i]; j++)
            cycles[idx++] = pk[i];

    /* t = sum_{i<j} gcd(c_i,c_j) + (1/2)*sum(c_i - 1) */
    int t = 0;
    for (int i = 0; i < num_cycles; i++)
        for (int j = i + 1; j < num_cycles; j++)
            t += gcd_func(cycles[i], cycles[j]);
    for (int i = 0; i < num_cycles; i++)
        t += (cycles[i] - 1) / 2;  /* c_i odd, so (c_i-1)/2 is exact */

    mpz_t z_lambda, contrib, power;
    mpz_init(z_lambda);
    mpz_init(contrib);
    mpz_init(power);
    mpz_set_ui(z_lambda, 1);
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < pm[i]; j++)
            mpz_mul_ui(z_lambda, z_lambda, pk[i]);
        for (int j = 2; j <= pm[i]; j++)
            mpz_mul_ui(z_lambda, z_lambda, j);
    }

    mpz_divexact(contrib, fact_n, z_lambda);
    mpz_set_ui(power, 2);
    mpz_pow_ui(power, power, t);
    mpz_mul(contrib, contrib, power);
    mpz_add(result, result, contrib);

    mpz_clear(contrib);
    mpz_clear(power);
    mpz_clear(z_lambda);
}

/* Enumerate partitions of 'remaining' into odd parts <= max_part */
static void enumerate(int remaining, int max_part, int *pk, int *pm, int depth) {
    if (remaining == 0) {
        process_partition(pk, pm, depth);
        return;
    }
    /* Only odd parts */
    int start = (max_part < remaining ? max_part : remaining);
    if (start % 2 == 0) start--;  /* make it odd */
    for (int p = start; p >= 1; p -= 2) {
        int max_m = remaining / p;
        for (int m = 1; m <= max_m; m++) {
            pk[depth] = p;
            pm[depth] = m;
            enumerate(remaining - m * p, p - 2, pk, pm, depth + 1);
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <max_n>\n", argv[0]);
        return 1;
    }

    int max_n = atoi(argv[1]);
    fprintf(stderr, "Computing A000568 (tournaments on n unlabeled nodes), n=0..%d\n\n", max_n);

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
        if (dt > 3600) { fprintf(stderr, "Stopping\n"); break; }
    }
    return 0;
}
