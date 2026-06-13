/*
 * Compute A000612: number of hypergraphs on n unlabeled nodes.
 *
 * a(n) = (1/n!) * sum_{lambda partition of n} permcount(lambda) * 2^{f(lambda)}
 * where f(lambda) = (1/L) * sum_{t=0}^{L-1} 2^{sum_i gcd(c_i, t)} - 1
 *   L = lcm of cycle lengths c_1,...,c_m
 *   The inner sum counts orbits on ALL subsets, minus 1 for empty set
 *
 * Also computes A003180 = 2 * A000612.
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o a000612_gmp a000612_gmp.c -lgmp
 *
 * Usage: ./a000612_gmp <max_n>
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
    return a / gcd_func(a, b) * b;
}

static int N;
static mpz_t result;
static mpz_t fact_n;

static void process_partition(int *pk, int *pm, int depth) {
    /* Expand cycle type */
    int num_cycles = 0;
    for (int i = 0; i < depth; i++) num_cycles += pm[i];
    int *cycles = (int *)malloc(num_cycles * sizeof(int));
    int idx = 0;
    for (int i = 0; i < depth; i++)
        for (int j = 0; j < pm[i]; j++)
            cycles[idx++] = pk[i];

    /* L = lcm of cycle lengths */
    long long L = 1;
    for (int i = 0; i < num_cycles; i++) {
        L = lcm_func(L, cycles[i]);
        if (L > 100000000LL) {
            fprintf(stderr, "Warning: L=%lld too large, skipping\n", L);
            free(cycles);
            return;
        }
    }

    /* Compute orbits on all subsets:
     * orbits = (1/L) * sum_{t=0}^{L-1} 2^{sum_i gcd(c_i, t)}
     * This is a huge number, but we only need the exponent for 2^(orbits-1)
     * orbits fits in a normal integer for reasonable n
     */
    mpz_t orbits_sum, power;
    mpz_init(orbits_sum);
    mpz_init(power);
    mpz_set_ui(orbits_sum, 0);

    for (long long t = 0; t < L; t++) {
        int exp = 0;
        for (int i = 0; i < num_cycles; i++) {
            /* gcd(c_i, t): since t < L and L = lcm(all c_i), t fits in int
             * gcd(c_i, 0) = c_i by convention (identity has all n fixed points) */
            exp += gcd_func(cycles[i], (int)(t % cycles[i]));
        }
        /* orbits_sum += 2^exp */
        mpz_set_ui(power, 1);
        mpz_mul_2exp(power, power, exp);
        mpz_add(orbits_sum, orbits_sum, power);
    }

    /* orbits = orbits_sum / L */
    mpz_divexact_ui(orbits_sum, orbits_sum, L);

    /* f = orbits - 1 (remove empty set orbit) */
    /* Need 2^f as the number of distinct hypergraphs fixed by this perm type */

    /* But orbits can be very large (2^n at identity), so 2^(orbits-1) is enormous */
    /* We need: contribution = (n!/z_lambda) * 2^(orbits-1) */

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

    /* contrib = (n! / z_lambda) * 2^(orbits-1) */
    mpz_t contrib;
    mpz_init(contrib);
    mpz_divexact(contrib, fact_n, z_lambda);

    /* orbits - 1 */
    mpz_sub_ui(orbits_sum, orbits_sum, 1);

    /* 2^(orbits-1) */
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
    mpz_clear(power);
    free(cycles);
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
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <max_n>\n", argv[0]);
        return 1;
    }

    int max_n = atoi(argv[1]);
    fprintf(stderr, "Computing A000612 (hypergraphs on n unlabeled nodes), n=0..%d\n\n", max_n);

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
        fprintf(stderr, "n=%d: %zu digits, %.3fs\n",
                N, mpz_sizeinbase(result, 10), dt);
        fflush(stdout);
        fflush(stderr);

        mpz_clear(result);
        mpz_clear(fact_n);

        if (dt > 300) {
            fprintf(stderr, "Stopping: too slow\n");
            break;
        }
    }

    return 0;
}
