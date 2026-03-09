/*
 * Compute A002724: number of inequivalent n x n binary matrices
 * under row and column permutations.
 *
 * a(n) = (1/(n!)^2) * sum_{lambda,mu partitions of n}
 *            |C(lambda)| * |C(mu)| * 2^{sum_{i,j} gcd(c_i, d_j)}
 *
 * where C(lambda) is the conjugacy class of cycle type lambda,
 * |C(lambda)| = n! / z_lambda, z_lambda = prod r^m * m!.
 *
 * Simplification: factor the sum as
 * a(n) = sum_{lambda} sum_{mu} (1/z_lambda) * (1/z_mu) * 2^{G(lambda,mu)}
 * where G(lambda,mu) = sum_{i,j} gcd(c_i, d_j)
 *     = sum_{r,s} m_r * m_s * gcd(r,s)   [with multiplicities]
 *
 * We iterate over all pairs of partitions of n.
 *
 * OEIS: A002724, had 51 terms (n=0..50)
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o binary_matrix_gmp binary_matrix_gmp.c -lgmp
 *
 * Usage: ./binary_matrix_gmp <max_n>
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
static mpz_t result, fact_n, fact_n_sq;

/* Store all partitions of N for enumeration */
#define MAX_PARTS 128
#define MAX_PARTITIONS 100000

static int part_pk[MAX_PARTITIONS][MAX_PARTS];
static int part_pm[MAX_PARTITIONS][MAX_PARTS];
static int part_depth[MAX_PARTITIONS];
static int num_partitions;

/* z_lambda values */
static mpz_t z_values[MAX_PARTITIONS];

static void store_partition(int *pk, int *pm, int depth) {
    if (num_partitions >= MAX_PARTITIONS) {
        fprintf(stderr, "Too many partitions!\n");
        return;
    }
    int idx = num_partitions++;
    for (int i = 0; i < depth; i++) {
        part_pk[idx][i] = pk[i];
        part_pm[idx][i] = pm[i];
    }
    part_depth[idx] = depth;

    /* Compute z_lambda */
    mpz_init(z_values[idx]);
    mpz_set_ui(z_values[idx], 1);
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < pm[i]; j++)
            mpz_mul_ui(z_values[idx], z_values[idx], pk[i]);
        for (int j = 2; j <= pm[i]; j++)
            mpz_mul_ui(z_values[idx], z_values[idx], j);
    }
}

static void enumerate(int remaining, int max_part, int *pk, int *pm, int depth) {
    if (remaining == 0) {
        store_partition(pk, pm, depth);
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
    fprintf(stderr, "Computing A002724 (n x n binary matrices under row/col perm), n=0..%d\n\n", max_n);

    printf("0 1\n");
    fflush(stdout);

    for (N = 1; N <= max_n; N++) {
        struct timespec t_start, t_end;
        clock_gettime(CLOCK_MONOTONIC, &t_start);

        mpz_init(result);
        mpz_init(fact_n);
        mpz_init(fact_n_sq);
        mpz_fac_ui(fact_n, N);
        mpz_mul(fact_n_sq, fact_n, fact_n);
        mpz_set_ui(result, 0);

        /* Enumerate all partitions of N */
        num_partitions = 0;
        int pk[64], pm[64];
        enumerate(N, N, pk, pm, 0);

        /* For each pair of partitions, compute contribution */
        mpz_t contrib, power, z_prod;
        mpz_init(contrib);
        mpz_init(power);
        mpz_init(z_prod);

        for (int a = 0; a < num_partitions; a++) {
            for (int b = 0; b < num_partitions; b++) {
                /* G = sum over all cycle pairs of gcd */
                /* Expand: cycles of lambda are pk[a][i] repeated pm[a][i] times,
                   cycles of mu are pk[b][j] repeated pm[b][j] times.
                   G = sum_{i,j} pm[a][i] * pm[b][j] * gcd(pk[a][i], pk[b][j]) */
                int G = 0;
                for (int i = 0; i < part_depth[a]; i++)
                    for (int j = 0; j < part_depth[b]; j++)
                        G += part_pm[a][i] * part_pm[b][j] * gcd_func(part_pk[a][i], part_pk[b][j]);

                /* contrib = (n!^2 / (z_a * z_b)) * 2^G */
                mpz_mul(z_prod, z_values[a], z_values[b]);
                mpz_divexact(contrib, fact_n_sq, z_prod);
                mpz_set_ui(power, 1);
                mpz_mul_2exp(power, power, G);
                mpz_mul(contrib, contrib, power);
                mpz_add(result, result, contrib);
            }
        }

        mpz_divexact(result, result, fact_n_sq);

        clock_gettime(CLOCK_MONOTONIC, &t_end);
        double dt = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

        gmp_printf("%d %Zd\n", N, result);
        fprintf(stderr, "n=%d: %zu digits, %.3fs (%d partitions)\n", N, mpz_sizeinbase(result, 10), dt, num_partitions);
        fflush(stdout);
        fflush(stderr);

        /* Cleanup */
        for (int i = 0; i < num_partitions; i++)
            mpz_clear(z_values[i]);
        mpz_clear(contrib);
        mpz_clear(power);
        mpz_clear(z_prod);
        mpz_clear(result);
        mpz_clear(fact_n);
        mpz_clear(fact_n_sq);
        if (dt > 300) { fprintf(stderr, "Stopping\n"); break; }
    }
    return 0;
}
