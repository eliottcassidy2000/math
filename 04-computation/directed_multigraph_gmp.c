/*
 * Compute number of directed k-multigraphs on n unlabeled nodes.
 *
 * a(n) = (1/n!) * sum_{lambda} permcount(lambda) * base^{directed_pair_orbits(lambda)}
 * where directed_pair_orbits = 2*sum_{i<j} gcd(c_i,c_j) + sum(c_i - 1)
 *
 * Known OEIS sequences:
 *   base=2: A000273 (digraphs)
 *   base=3: A053467 (directed 2-multigraphs)
 *   base=4: A053468 (directed 3-multigraphs)
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o directed_multigraph_gmp directed_multigraph_gmp.c -lgmp
 *
 * Usage: ./directed_multigraph_gmp <base> <max_n>
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

static int BASE, N;
static mpz_t result, fact_n;

static void process_partition(int *pk, int *pm, int depth) {
    int num_cycles = 0;
    for (int i = 0; i < depth; i++) num_cycles += pm[i];
    int cycles[128];
    int idx = 0;
    for (int i = 0; i < depth; i++)
        for (int j = 0; j < pm[i]; j++)
            cycles[idx++] = pk[i];

    /* Directed pair orbits = 2*sum_{i<j} gcd + sum(c_i - 1) */
    int dir_orbits = 0;
    for (int i = 0; i < num_cycles; i++)
        for (int j = i + 1; j < num_cycles; j++)
            dir_orbits += 2 * gcd_func(cycles[i], cycles[j]);
    for (int i = 0; i < num_cycles; i++)
        dir_orbits += cycles[i] - 1;

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
    mpz_set_ui(power, BASE);
    mpz_pow_ui(power, power, dir_orbits);
    mpz_mul(contrib, contrib, power);
    mpz_add(result, result, contrib);

    mpz_clear(contrib);
    mpz_clear(power);
    mpz_clear(z_lambda);
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
        fprintf(stderr, "Usage: %s <base> <max_n>\n", argv[0]);
        fprintf(stderr, "  base=2: digraphs (A000273)\n");
        fprintf(stderr, "  base=3: directed 2-multigraphs (A053467)\n");
        fprintf(stderr, "  base=4: directed 3-multigraphs (A053468)\n");
        return 1;
    }

    BASE = atoi(argv[1]);
    int max_n = atoi(argv[2]);

    fprintf(stderr, "Computing directed %d-multigraphs (base=%d), n=0..%d\n\n",
            BASE - 1, BASE, max_n);

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
