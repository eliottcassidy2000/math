/*
 * Compute A002724: inequivalent n x n binary matrices under row/col perm.
 * FAST version using single-partition + generating function approach.
 *
 * a(n) = (1/n!) * sum over partitions q of n:
 *           permcount(q) * [x^n] exp(sum_{t=1}^{n} 2^{K(q,t)} * x^t / t)
 *
 * where K(q,t) = sum_i gcd(t, q_i).
 *
 * Uses recurrence: S(0) = 1,
 *   S(k) = sum_{j=1}^{k} c_j * S(k-j) * (k-1)!/(k-j)!
 * with c_t = 2^{K(q,t)}.
 *
 * This is O(p(n)*n^2) instead of O(p(n)^2*d^2) for the pair approach.
 *
 * OEIS: A002724, had 51 terms
 *
 * Also outputs A006383 (with column complementation) and
 * A007139 = (A002724 + A122082) / 2 (bicolored bipartite graphs).
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o a002724_fast_gmp a002724_fast_gmp.c -lgmp
 *
 * Usage: ./a002724_fast_gmp <max_n>
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
static mpz_t total, fact_n, fact_n_sq;

static int pk[64], pm[64];

static mpz_t S_arr[512];
static mpz_t c_arr[512];
static mpz_t temp, z_lambda, pcount;

static void process_partition(int depth) {
    /* c_t = 2^{K(q,t)} for standard binary matrices */
    for (int t = 1; t <= N; t++) {
        int K = 0;
        for (int i = 0; i < depth; i++)
            K += pm[i] * gcd_func(t, pk[i]);
        mpz_set_ui(c_arr[t], 1);
        mpz_mul_2exp(c_arr[t], c_arr[t], K);
    }

    /* S(k) recurrence */
    mpz_set_ui(S_arr[0], 1);
    for (int k = 1; k <= N; k++) {
        mpz_set_ui(S_arr[k], 0);
        mpz_t falling;
        mpz_init(falling);
        mpz_set_ui(falling, 1);
        for (int j = 1; j <= k; j++) {
            if (j >= 2) mpz_mul_ui(falling, falling, k - j + 1);
            mpz_mul(temp, c_arr[j], S_arr[k - j]);
            mpz_mul(temp, temp, falling);
            mpz_add(S_arr[k], S_arr[k], temp);
        }
        mpz_clear(falling);
    }

    /* permcount = n!/z */
    mpz_set_ui(z_lambda, 1);
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < pm[i]; j++)
            mpz_mul_ui(z_lambda, z_lambda, pk[i]);
        for (int j = 2; j <= pm[i]; j++)
            mpz_mul_ui(z_lambda, z_lambda, j);
    }
    mpz_divexact(pcount, fact_n, z_lambda);

    mpz_mul(temp, pcount, S_arr[N]);
    mpz_add(total, total, temp);
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
    fprintf(stderr, "Computing A002724 (n x n binary matrices, fast), n=0..%d\n\n", max_n);

    mpz_init(total); mpz_init(fact_n); mpz_init(fact_n_sq);
    mpz_init(temp); mpz_init(z_lambda); mpz_init(pcount);
    for (int i = 0; i < 512; i++) { mpz_init(S_arr[i]); mpz_init(c_arr[i]); }

    printf("0 1\n");
    fflush(stdout);

    for (N = 1; N <= max_n; N++) {
        struct timespec t_start, t_end;
        clock_gettime(CLOCK_MONOTONIC, &t_start);

        mpz_fac_ui(fact_n, N);
        mpz_mul(fact_n_sq, fact_n, fact_n);
        mpz_set_ui(total, 0);

        enumerate(N, N, 0);

        mpz_t result;
        mpz_init(result);
        mpz_divexact(result, total, fact_n_sq);

        clock_gettime(CLOCK_MONOTONIC, &t_end);
        double dt = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

        gmp_printf("%d %Zd\n", N, result);
        fprintf(stderr, "n=%d: %zu digits, %.3fs\n", N, mpz_sizeinbase(result, 10), dt);
        fflush(stdout); fflush(stderr);

        mpz_clear(result);
        if (dt > 3600) { fprintf(stderr, "Stopping\n"); break; }
    }

    mpz_clear(total); mpz_clear(fact_n); mpz_clear(fact_n_sq);
    mpz_clear(temp); mpz_clear(z_lambda); mpz_clear(pcount);
    for (int i = 0; i < 512; i++) { mpz_clear(S_arr[i]); mpz_clear(c_arr[i]); }
    return 0;
}
