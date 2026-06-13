/*
 * Compute A006383: equivalence classes of n x n binary matrices
 * under row permutations, column permutations, and independent
 * column complementation (hyperoctahedral group B_n on columns).
 *
 * Formula (from A363349 PARI code):
 * a(n) = (1/n!) * sum over row partitions q of n:
 *           permcount(q) * [x^n] exp(sum_{t=1}^{n} c_t * x^t / t)
 *
 * where c_t = 2^{K(q,t) - [t % e == 0]}
 *       K(q,t) = sum_i gcd(t, q_i)  (summed over all cycle lengths)
 *       e = 2^{min_i v_2(q_i)}  (2 to the min 2-adic valuation of parts)
 *
 * Uses integer recurrence: S(0) = 1,
 *   S(k) = sum_{j=1}^{k} c_j * S(k-j) * (k-1)! / (k-j)!
 * Then contribution = permcount(q) * S(n), and a(n) = total / (n!)^2.
 *
 * OEIS: A006383, had 50 terms
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o a006383_gmp a006383_gmp.c -lgmp
 *
 * Usage: ./a006383_gmp <max_n>
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
static mpz_t total, fact_n, fact_n_sq;

/* Partition storage */
#define MAX_PARTS 128
static int pk[64], pm[64];

/* S array and c_t array for the recurrence */
static mpz_t S_arr[512];   /* S(0), ..., S(n) */
static mpz_t c_arr[512];   /* c_1, ..., c_n */
static mpz_t temp, z_lambda, pcount;

static void process_partition(int depth) {
    /* Compute e = 2^{min v_2(pk[i])} */
    int min_v2 = 100;
    for (int i = 0; i < depth; i++) {
        int v = v2(pk[i]);
        if (v < min_v2) min_v2 = v;
    }
    int e = 1 << min_v2;  /* 2^min_v2 */

    /* Compute c_t for t = 1, ..., N */
    for (int t = 1; t <= N; t++) {
        int K = 0;
        for (int i = 0; i < depth; i++)
            K += pm[i] * gcd_func(t, pk[i]);
        int flag = (t % e == 0) ? 1 : 0;
        int exp_val = K - flag;
        mpz_set_ui(c_arr[t], 1);
        mpz_mul_2exp(c_arr[t], c_arr[t], exp_val);
    }

    /* Compute S(k) = sum_{j=1}^{k} c_j * S(k-j) * (k-1)!/(k-j)! */
    mpz_set_ui(S_arr[0], 1);
    for (int k = 1; k <= N; k++) {
        mpz_set_ui(S_arr[k], 0);
        /* binom_coeff = (k-1)!/(k-j)! = product of (k-j+1) to (k-1) */
        /* For j=1: (k-1)!/(k-1)! = 1 */
        /* For j=2: k-1 */
        /* For j=l: (k-1)(k-2)...(k-l+1) */
        mpz_t falling;
        mpz_init(falling);
        mpz_set_ui(falling, 1);  /* (k-1)!/(k-1)! for j=1 */
        for (int j = 1; j <= k; j++) {
            /* falling = (k-1)!/(k-j)! */
            if (j >= 2) {
                mpz_mul_ui(falling, falling, k - j + 1);
            }
            /* S(k) += c_j * S(k-j) * falling */
            mpz_mul(temp, c_arr[j], S_arr[k - j]);
            mpz_mul(temp, temp, falling);
            mpz_add(S_arr[k], S_arr[k], temp);
        }
        mpz_clear(falling);
    }

    /* permcount = n! / z_lambda */
    mpz_set_ui(z_lambda, 1);
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < pm[i]; j++)
            mpz_mul_ui(z_lambda, z_lambda, pk[i]);
        for (int j = 2; j <= pm[i]; j++)
            mpz_mul_ui(z_lambda, z_lambda, j);
    }
    mpz_divexact(pcount, fact_n, z_lambda);

    /* total += permcount * S(n) */
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
    fprintf(stderr, "Computing A006383 (binary matrices, row/col perm + col complement), n=0..%d\n\n", max_n);

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
