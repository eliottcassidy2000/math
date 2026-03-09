/*
 * Compute graphs with self-loops on n unlabeled nodes.
 *
 * A000666: a(n) = (1/n!) * sum_{lambda} permcount(lambda) * 2^{pair_orbits + num_cycles}
 *   pair_orbits = sum_{i<j} gcd(c_i,c_j) + sum floor(c_i/2)
 *   num_cycles = number of cycles in the partition
 *   (self-loops add one orbit per cycle)
 *
 * Also computes:
 *   A000666 (symmetric relations / graphs with self-loops, base=2)
 *   A000595 (binary relations / digraphs with self-loops, base=2, directed pair orbits + num_cycles)
 *
 * Mode 1 (undirected, symmetric relations):
 *   orbits = sum_{i<j} gcd(c_i,c_j) + sum floor(c_i/2) + num_cycles
 *   weight = 2^orbits
 *
 * Mode 2 (directed, reflexive relations):
 *   orbits = 2*sum_{i<j} gcd(c_i,c_j) + sum(c_i-1) + num_cycles
 *   weight = 2^orbits
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o selfloop_graph_gmp selfloop_graph_gmp.c -lgmp
 *
 * Usage: ./selfloop_graph_gmp <mode> <max_n>
 *   mode=1: A000666 (symmetric relations)
 *   mode=2: A000595 (binary relations / digraphs with self-loops)
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

static int MODE, N;
static mpz_t result, fact_n;

static void process_partition(int *pk, int *pm, int depth) {
    int num_cycles = 0;
    for (int i = 0; i < depth; i++) num_cycles += pm[i];
    int cycles[128];
    int idx = 0;
    for (int i = 0; i < depth; i++)
        for (int j = 0; j < pm[i]; j++)
            cycles[idx++] = pk[i];

    int orbits = 0;
    if (MODE == 1) {
        /* Undirected: pair_orbits + self-loop orbits */
        for (int i = 0; i < num_cycles; i++)
            for (int j = i + 1; j < num_cycles; j++)
                orbits += gcd_func(cycles[i], cycles[j]);
        for (int i = 0; i < num_cycles; i++)
            orbits += cycles[i] / 2;
        orbits += num_cycles;  /* self-loop orbits */
    } else {
        /* Directed: directed pair_orbits + self-loop orbits */
        for (int i = 0; i < num_cycles; i++)
            for (int j = i + 1; j < num_cycles; j++)
                orbits += 2 * gcd_func(cycles[i], cycles[j]);
        for (int i = 0; i < num_cycles; i++)
            orbits += cycles[i] - 1;
        orbits += num_cycles;  /* self-loop orbits */
    }

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
    mpz_set_ui(power, 1);
    mpz_mul_2exp(power, power, orbits);  /* 2^orbits via bit shift */
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
        fprintf(stderr, "Usage: %s <mode> <max_n>\n", argv[0]);
        fprintf(stderr, "  mode=1: A000666 (symmetric relations)\n");
        fprintf(stderr, "  mode=2: A002854 (reflexive digraphs)\n");
        return 1;
    }

    MODE = atoi(argv[1]);
    int max_n = atoi(argv[2]);

    const char *name = (MODE == 1) ? "A000666 (symmetric relations)" :
                                      "A002854 (reflexive digraphs)";
    fprintf(stderr, "Computing %s, n=0..%d\n\n", name, max_n);

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
