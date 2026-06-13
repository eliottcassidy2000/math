/*
 * Compute FOUR sequences simultaneously:
 *   A002724: binary n×n matrices under row/col perm
 *   A006383: binary n×n matrices under row/col perm + col complement (B_n)
 *   A122082: bicolored graphs invariant under color swap
 *   A007139: bicolored bipartite graphs = (A002724 + A122082)/2
 *
 * All use single-partition GF approach for speed.
 *
 * A002724: c_t = 2^{K(q,t)}
 * A006383: c_t = 2^{K(q,t) - [t % e == 0]} where e = 2^{min v_2(q_i)}
 * A122082: single partition formula: (1/n!) * sum permcount(q) * 2^{edges(q)}
 *   where edges(q) = sum_{i<j} gcd(c_i,c_j) + sum_i floor((c_i+1)/2)
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o matrix_quad_gmp matrix_quad_gmp.c -lgmp
 *
 * Usage: ./matrix_quad_gmp <max_n>
 * Output: n A002724 A006383 A122082 A007139
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

static int v2(int x) {
    if (x == 0) return 0;
    int v = 0;
    while (x % 2 == 0) { x /= 2; v++; }
    return v;
}

static int N;
static mpz_t total_724, total_383, total_122;
static mpz_t fact_n, fact_n_sq;

static int pk[64], pm[64];
static mpz_t S724[512], S383[512], c724[512], c383[512];
static mpz_t temp, z_lambda, pcount;

static void process_partition(int depth) {
    /* For A122082: single partition formula */
    int num_cycles = 0;
    for (int i = 0; i < depth; i++) num_cycles += pm[i];
    int cycles[512];
    int idx = 0;
    for (int i = 0; i < depth; i++)
        for (int j = 0; j < pm[i]; j++)
            cycles[idx++] = pk[i];

    int edges = 0;
    for (int i = 0; i < num_cycles; i++)
        for (int j = i + 1; j < num_cycles; j++)
            edges += gcd_func(cycles[i], cycles[j]);
    for (int i = 0; i < num_cycles; i++)
        edges += (cycles[i] + 1) / 2;

    /* z_lambda and permcount */
    mpz_set_ui(z_lambda, 1);
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < pm[i]; j++)
            mpz_mul_ui(z_lambda, z_lambda, pk[i]);
        for (int j = 2; j <= pm[i]; j++)
            mpz_mul_ui(z_lambda, z_lambda, j);
    }
    mpz_divexact(pcount, fact_n, z_lambda);

    /* A122082 contribution */
    mpz_set_ui(temp, 1);
    mpz_mul_2exp(temp, temp, edges);
    mpz_mul(temp, temp, pcount);
    mpz_add(total_122, total_122, temp);

    /* For A002724 and A006383: GF approach */
    int min_v2_val = 100;
    for (int i = 0; i < depth; i++) {
        int vv = v2(pk[i]);
        if (vv < min_v2_val) min_v2_val = vv;
    }
    int e = 1 << min_v2_val;

    for (int t = 1; t <= N; t++) {
        int K = 0;
        for (int i = 0; i < depth; i++)
            K += pm[i] * gcd_func(t, pk[i]);
        mpz_set_ui(c724[t], 1);
        mpz_mul_2exp(c724[t], c724[t], K);
        int flag = (t % e == 0) ? 1 : 0;
        mpz_set_ui(c383[t], 1);
        mpz_mul_2exp(c383[t], c383[t], K - flag);
    }

    /* S recurrences for both A002724 and A006383 */
    mpz_set_ui(S724[0], 1);
    mpz_set_ui(S383[0], 1);
    for (int k = 1; k <= N; k++) {
        mpz_set_ui(S724[k], 0);
        mpz_set_ui(S383[k], 0);
        mpz_t falling;
        mpz_init(falling);
        mpz_set_ui(falling, 1);
        for (int j = 1; j <= k; j++) {
            if (j >= 2) mpz_mul_ui(falling, falling, k - j + 1);

            mpz_mul(temp, c724[j], S724[k - j]);
            mpz_mul(temp, temp, falling);
            mpz_add(S724[k], S724[k], temp);

            mpz_mul(temp, c383[j], S383[k - j]);
            mpz_mul(temp, temp, falling);
            mpz_add(S383[k], S383[k], temp);
        }
        mpz_clear(falling);
    }

    /* Add to totals */
    mpz_mul(temp, pcount, S724[N]);
    mpz_add(total_724, total_724, temp);

    mpz_mul(temp, pcount, S383[N]);
    mpz_add(total_383, total_383, temp);
}

static void enumerate(int remaining, int max_part, int depth) {
    if (remaining == 0) { process_partition(depth); return; }
    for (int p = (max_part < remaining ? max_part : remaining); p >= 1; p--) {
        int max_m = remaining / p;
        for (int m = 1; m <= max_m; m++) {
            pk[depth] = p; pm[depth] = m;
            enumerate(remaining - m * p, p - 1, depth + 1);
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <max_n>\n", argv[0]);
        return 1;
    }

    int max_n = atoi(argv[1]);
    fprintf(stderr, "Computing A002724, A006383, A122082, A007139 simultaneously, n=0..%d\n\n", max_n);

    mpz_init(total_724); mpz_init(total_383); mpz_init(total_122);
    mpz_init(fact_n); mpz_init(fact_n_sq);
    mpz_init(temp); mpz_init(z_lambda); mpz_init(pcount);
    for (int i = 0; i < 512; i++) {
        mpz_init(S724[i]); mpz_init(S383[i]);
        mpz_init(c724[i]); mpz_init(c383[i]);
    }

    printf("0 1 1 1 1\n");
    fflush(stdout);

    for (N = 1; N <= max_n; N++) {
        struct timespec t_start, t_end;
        clock_gettime(CLOCK_MONOTONIC, &t_start);

        mpz_fac_ui(fact_n, N);
        mpz_mul(fact_n_sq, fact_n, fact_n);
        mpz_set_ui(total_724, 0);
        mpz_set_ui(total_383, 0);
        mpz_set_ui(total_122, 0);

        enumerate(N, N, 0);

        mpz_t r724, r383, r122, r139;
        mpz_init(r724); mpz_init(r383); mpz_init(r122); mpz_init(r139);

        mpz_divexact(r724, total_724, fact_n_sq);
        mpz_divexact(r383, total_383, fact_n_sq);
        mpz_divexact(r122, total_122, fact_n);

        /* A007139 = (A002724 + A122082) / 2 */
        mpz_add(r139, r724, r122);
        mpz_divexact_ui(r139, r139, 2);

        clock_gettime(CLOCK_MONOTONIC, &t_end);
        double dt = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

        gmp_printf("%d %Zd %Zd %Zd %Zd\n", N, r724, r383, r122, r139);
        fprintf(stderr, "n=%d: %.3fs\n", N, dt);
        fflush(stdout); fflush(stderr);

        mpz_clear(r724); mpz_clear(r383); mpz_clear(r122); mpz_clear(r139);
        if (dt > 3600) { fprintf(stderr, "Stopping\n"); break; }
    }

    mpz_clear(total_724); mpz_clear(total_383); mpz_clear(total_122);
    mpz_clear(fact_n); mpz_clear(fact_n_sq);
    mpz_clear(temp); mpz_clear(z_lambda); mpz_clear(pcount);
    for (int i = 0; i < 512; i++) {
        mpz_clear(S724[i]); mpz_clear(S383[i]);
        mpz_clear(c724[i]); mpz_clear(c383[i]);
    }
    return 0;
}
