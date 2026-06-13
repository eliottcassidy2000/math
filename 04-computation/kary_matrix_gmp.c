/*
 * Compute number of inequivalent n x n k-ary matrices
 * under row and column permutations.
 *
 * a(n) = (1/(n!)^2) * sum_{lambda,mu partitions of n}
 *            |C(lambda)| * |C(mu)| * k^{sum_{i,j} gcd(c_i, d_j)}
 *
 * k=2: A002724 (binary matrices), OEIS had 51 terms
 * k=3: A006380 (ternary matrices)?
 * k=4: A006381?
 * k=5: ?
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o kary_matrix_gmp kary_matrix_gmp.c -lgmp
 *
 * Usage: ./kary_matrix_gmp <k> <max_n>
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

static int K, N;
static mpz_t result, fact_n, fact_n_sq;

#define MAX_PARTS 128
#define MAX_PARTITIONS 100000

static int ppk[MAX_PARTITIONS][MAX_PARTS];
static int ppm[MAX_PARTITIONS][MAX_PARTS];
static int pdepth[MAX_PARTITIONS];
static int nparts;
static mpz_t zvals[MAX_PARTITIONS];

static void store_part(int *pk, int *pm, int depth) {
    if (nparts >= MAX_PARTITIONS) return;
    int idx = nparts++;
    for (int i = 0; i < depth; i++) {
        ppk[idx][i] = pk[i];
        ppm[idx][i] = pm[i];
    }
    pdepth[idx] = depth;
    mpz_init(zvals[idx]);
    mpz_set_ui(zvals[idx], 1);
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < pm[i]; j++)
            mpz_mul_ui(zvals[idx], zvals[idx], pk[i]);
        for (int j = 2; j <= pm[i]; j++)
            mpz_mul_ui(zvals[idx], zvals[idx], j);
    }
}

static void enumerate(int remaining, int max_part, int *pk, int *pm, int depth) {
    if (remaining == 0) { store_part(pk, pm, depth); return; }
    for (int p = (max_part < remaining ? max_part : remaining); p >= 1; p--) {
        int max_m = remaining / p;
        for (int m = 1; m <= max_m; m++) {
            pk[depth] = p; pm[depth] = m;
            enumerate(remaining - m * p, p - 1, pk, pm, depth + 1);
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <k> <max_n>\n", argv[0]);
        fprintf(stderr, "  k=2: A002724 (binary matrices)\n");
        fprintf(stderr, "  k=3,4,...: k-ary matrices under row/col perm\n");
        return 1;
    }

    K = atoi(argv[1]);
    int max_n = atoi(argv[2]);
    fprintf(stderr, "Computing %d-ary n x n matrices under row/col perm, n=0..%d\n\n", K, max_n);

    printf("0 1\n");
    fflush(stdout);

    for (N = 1; N <= max_n; N++) {
        struct timespec t_start, t_end;
        clock_gettime(CLOCK_MONOTONIC, &t_start);

        mpz_init(result); mpz_init(fact_n); mpz_init(fact_n_sq);
        mpz_fac_ui(fact_n, N);
        mpz_mul(fact_n_sq, fact_n, fact_n);
        mpz_set_ui(result, 0);

        nparts = 0;
        int pk[64], pm[64];
        enumerate(N, N, pk, pm, 0);

        mpz_t contrib, power, z_prod;
        mpz_init(contrib); mpz_init(power); mpz_init(z_prod);

        for (int a = 0; a < nparts; a++) {
            for (int b = 0; b < nparts; b++) {
                int G = 0;
                for (int i = 0; i < pdepth[a]; i++)
                    for (int j = 0; j < pdepth[b]; j++)
                        G += ppm[a][i] * ppm[b][j] * gcd_func(ppk[a][i], ppk[b][j]);

                mpz_mul(z_prod, zvals[a], zvals[b]);
                mpz_divexact(contrib, fact_n_sq, z_prod);
                mpz_set_ui(power, K);
                mpz_pow_ui(power, power, G);
                mpz_mul(contrib, contrib, power);
                mpz_add(result, result, contrib);
            }
        }

        mpz_divexact(result, result, fact_n_sq);

        clock_gettime(CLOCK_MONOTONIC, &t_end);
        double dt = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

        gmp_printf("%d %Zd\n", N, result);
        fprintf(stderr, "n=%d: %zu digits, %.3fs\n", N, mpz_sizeinbase(result, 10), dt);
        fflush(stdout); fflush(stderr);

        for (int i = 0; i < nparts; i++) mpz_clear(zvals[i]);
        mpz_clear(contrib); mpz_clear(power); mpz_clear(z_prod);
        mpz_clear(result); mpz_clear(fact_n); mpz_clear(fact_n_sq);
        if (dt > 300) { fprintf(stderr, "Stopping\n"); break; }
    }
    return 0;
}
