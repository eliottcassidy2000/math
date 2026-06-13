/*
 * Compute number of inequivalent m x n binary matrices
 * under row and column permutations.
 *
 * a(m,n) = (1/(m!*n!)) * sum_{lambda part of m} sum_{mu part of n}
 *            |C(lambda)| * |C(mu)| * 2^{sum_{i,j} gcd(c_i, d_j)}
 *
 * For m=n: A002724 (square binary matrices)
 * For the diagonal m=n of A028657
 *
 * Also computes the full triangle A028657(m,n) for fixed m.
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o binary_matrix_mn_gmp binary_matrix_mn_gmp.c -lgmp
 *
 * Usage: ./binary_matrix_mn_gmp <max_m> <max_n>
 *   Computes a(m,n) for all 1 <= m <= max_m, 1 <= n <= max_n
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

#define MAX_PARTS 128
#define MAX_PARTITIONS 100000

typedef struct {
    int pk[MAX_PARTS];
    int pm[MAX_PARTS];
    int depth;
    mpz_t z_val;
} partition_t;

static partition_t *partitions_m = NULL;
static partition_t *partitions_n = NULL;
static int num_parts_m, num_parts_n;

static void enum_parts(int remaining, int max_part, int *pk, int *pm, int depth,
                       partition_t *store, int *count) {
    if (remaining == 0) {
        int idx = (*count)++;
        for (int i = 0; i < depth; i++) {
            store[idx].pk[i] = pk[i];
            store[idx].pm[i] = pm[i];
        }
        store[idx].depth = depth;
        mpz_init(store[idx].z_val);
        mpz_set_ui(store[idx].z_val, 1);
        for (int i = 0; i < depth; i++) {
            for (int j = 0; j < pm[i]; j++)
                mpz_mul_ui(store[idx].z_val, store[idx].z_val, pk[i]);
            for (int j = 2; j <= pm[i]; j++)
                mpz_mul_ui(store[idx].z_val, store[idx].z_val, j);
        }
        return;
    }
    for (int p = (max_part < remaining ? max_part : remaining); p >= 1; p--) {
        int max_mul = remaining / p;
        for (int m = 1; m <= max_mul; m++) {
            pk[depth] = p;
            pm[depth] = m;
            enum_parts(remaining - m * p, p - 1, pk, pm, depth + 1, store, count);
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <max_m> <max_n>\n", argv[0]);
        fprintf(stderr, "  Computes A028657(m,n) for m=0..max_m, n=0..max_n\n");
        fprintf(stderr, "  Square matrices (m=n): A002724\n");
        return 1;
    }

    int max_m = atoi(argv[1]);
    int max_n = atoi(argv[2]);

    fprintf(stderr, "Computing binary matrices under row/col perm, m=0..%d, n=0..%d\n\n", max_m, max_n);

    /* For each (m,n) pair, we need partitions of m and partitions of n */
    /* Pre-allocate storage for max size */
    int alloc_sz = MAX_PARTITIONS;
    partitions_m = calloc(alloc_sz, sizeof(partition_t));
    partitions_n = calloc(alloc_sz, sizeof(partition_t));

    mpz_t result, fact_m, fact_n, fact_mn, contrib, power, z_prod;
    mpz_init(result); mpz_init(fact_m); mpz_init(fact_n);
    mpz_init(fact_mn); mpz_init(contrib); mpz_init(power); mpz_init(z_prod);

    for (int m = 0; m <= max_m; m++) {
        for (int n = 0; n <= max_n; n++) {
            struct timespec t_start, t_end;
            clock_gettime(CLOCK_MONOTONIC, &t_start);

            if (m == 0 || n == 0) {
                gmp_printf("%d %d 1\n", m, n);
                fflush(stdout);
                continue;
            }

            /* Enumerate partitions of m */
            num_parts_m = 0;
            int pk[64], pm[64];
            enum_parts(m, m, pk, pm, 0, partitions_m, &num_parts_m);

            /* Enumerate partitions of n */
            num_parts_n = 0;
            enum_parts(n, n, pk, pm, 0, partitions_n, &num_parts_n);

            mpz_fac_ui(fact_m, m);
            mpz_fac_ui(fact_n, n);
            mpz_mul(fact_mn, fact_m, fact_n);

            mpz_set_ui(result, 0);

            for (int a = 0; a < num_parts_m; a++) {
                for (int b = 0; b < num_parts_n; b++) {
                    int G = 0;
                    for (int i = 0; i < partitions_m[a].depth; i++)
                        for (int j = 0; j < partitions_n[b].depth; j++)
                            G += partitions_m[a].pm[i] * partitions_n[b].pm[j] *
                                 gcd_func(partitions_m[a].pk[i], partitions_n[b].pk[j]);

                    mpz_mul(z_prod, partitions_m[a].z_val, partitions_n[b].z_val);
                    mpz_mul(contrib, fact_m, fact_n);
                    mpz_divexact(contrib, contrib, z_prod);
                    mpz_set_ui(power, 1);
                    mpz_mul_2exp(power, power, G);
                    mpz_mul(contrib, contrib, power);
                    mpz_add(result, result, contrib);
                }
            }

            mpz_divexact(result, result, fact_mn);

            clock_gettime(CLOCK_MONOTONIC, &t_end);
            double dt = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

            gmp_printf("%d %d %Zd\n", m, n, result);
            if (m == n)
                fprintf(stderr, "a(%d,%d): %zu digits, %.3fs\n", m, n, mpz_sizeinbase(result, 10), dt);
            fflush(stdout);
            fflush(stderr);

            /* Cleanup partition z_values */
            for (int i = 0; i < num_parts_m; i++) mpz_clear(partitions_m[i].z_val);
            for (int i = 0; i < num_parts_n; i++) mpz_clear(partitions_n[i].z_val);

            if (dt > 300) { fprintf(stderr, "Stopping\n"); goto done; }
        }
    }
done:
    mpz_clear(result); mpz_clear(fact_m); mpz_clear(fact_n);
    mpz_clear(fact_mn); mpz_clear(contrib); mpz_clear(power); mpz_clear(z_prod);
    free(partitions_m); free(partitions_n);
    return 0;
}
