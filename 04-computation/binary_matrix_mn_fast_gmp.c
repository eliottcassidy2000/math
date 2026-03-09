/*
 * Compute A028657: inequivalent m x n binary matrices under row/col perm.
 * FAST version using single-partition GF approach on columns.
 *
 * a(m,n) = (1/m!) * sum over row partitions q of m:
 *             permcount(q) * [x^n] exp(sum_{t=1}^{n} 2^{K(q,t)} * x^t / t)
 *
 * where K(q,t) = sum_i gcd(t, q_i) over all cycle lengths q_i of partition q.
 *
 * For the triangle A028657, this outputs (m,n,value) for all pairs.
 * Main diagonal (m=n) gives A002724.
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o binary_matrix_mn_fast_gmp binary_matrix_mn_fast_gmp.c -lgmp
 *
 * Usage: ./binary_matrix_mn_fast_gmp <max_m> <max_n>
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

#define MAX_N 512
#define MAX_PARTITIONS 100000

typedef struct {
    int pk[128];
    int pm[128];
    int depth;
    mpz_t z_val;
    mpz_t pcount;
} partition_t;

static partition_t parts[MAX_PARTITIONS];
static int num_parts;
static int pk_temp[64], pm_temp[64];

static void store_part(int depth, mpz_t fact_m) {
    int idx = num_parts++;
    for (int i = 0; i < depth; i++) {
        parts[idx].pk[i] = pk_temp[i];
        parts[idx].pm[i] = pm_temp[i];
    }
    parts[idx].depth = depth;
    mpz_init(parts[idx].z_val);
    mpz_set_ui(parts[idx].z_val, 1);
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < pm_temp[i]; j++)
            mpz_mul_ui(parts[idx].z_val, parts[idx].z_val, pk_temp[i]);
        for (int j = 2; j <= pm_temp[i]; j++)
            mpz_mul_ui(parts[idx].z_val, parts[idx].z_val, j);
    }
    mpz_init(parts[idx].pcount);
    mpz_divexact(parts[idx].pcount, fact_m, parts[idx].z_val);
}

static mpz_t fact_temp;

static void enumerate(int remaining, int max_part, int depth) {
    if (remaining == 0) { store_part(depth, fact_temp); return; }
    for (int p = (max_part < remaining ? max_part : remaining); p >= 1; p--) {
        int max_m = remaining / p;
        for (int m = 1; m <= max_m; m++) {
            pk_temp[depth] = p; pm_temp[depth] = m;
            enumerate(remaining - m * p, p - 1, depth + 1);
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <max_m> <max_n>\n", argv[0]);
        return 1;
    }

    int max_m = atoi(argv[1]);
    int max_n = atoi(argv[2]);
    fprintf(stderr, "Computing A028657 (m x n binary matrices, fast), m=0..%d, n=0..%d\n\n", max_m, max_n);

    mpz_t total, fact_n_sq, result, temp;
    mpz_t S_arr[MAX_N], c_arr[MAX_N];
    mpz_init(total); mpz_init(fact_n_sq); mpz_init(result); mpz_init(temp);
    mpz_init(fact_temp);
    for (int i = 0; i < MAX_N; i++) { mpz_init(S_arr[i]); mpz_init(c_arr[i]); }

    for (int m = 0; m <= max_m; m++) {
        /* Handle m=0 */
        if (m == 0) {
            for (int n = 0; n <= max_n; n++) {
                gmp_printf("%d %d 1\n", m, n);
            }
            fflush(stdout);
            continue;
        }

        /* Enumerate partitions of m */
        num_parts = 0;
        mpz_fac_ui(fact_temp, m);
        enumerate(m, m, 0);

        for (int n = 0; n <= max_n; n++) {
            if (n == 0) {
                gmp_printf("%d %d 1\n", m, n);
                fflush(stdout);
                continue;
            }

            struct timespec t_start, t_end;
            clock_gettime(CLOCK_MONOTONIC, &t_start);

            mpz_t fact_n;
            mpz_init(fact_n);
            mpz_fac_ui(fact_n, n);

            mpz_set_ui(total, 0);

            for (int a = 0; a < num_parts; a++) {
                /* For each row partition, compute GF for column contribution */
                for (int t = 1; t <= n; t++) {
                    int K = 0;
                    for (int i = 0; i < parts[a].depth; i++)
                        K += parts[a].pm[i] * gcd_func(t, parts[a].pk[i]);
                    mpz_set_ui(c_arr[t], 1);
                    mpz_mul_2exp(c_arr[t], c_arr[t], K);
                }

                /* S recurrence */
                mpz_set_ui(S_arr[0], 1);
                for (int k = 1; k <= n; k++) {
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

                /* total += permcount(q) * S(n) */
                mpz_mul(temp, parts[a].pcount, S_arr[n]);
                mpz_add(total, total, temp);
            }

            /* Divide by (m! * n!) but total already has m!/z * S(n).
               We need total / (n!^2) because S(n) = n! * P(n) and
               we want (1/(m!*n!)) * sum permcount * n! * P(n) / n! = total / n!^2...

               Actually: a(m,n) = (1/m!) * sum_q permcount(q) * P_q(n)
               where P_q(n) = S_q(n) / n!
               total = sum_q permcount(q) * S_q(n)
               a(m,n) = total / (m! * n!)
               But permcount(q) = m!/z(q), so total = sum m!/z(q) * S(n)
               a(m,n) = total / (m! * n!)
               Hmm wait, (1/m!) * sum (m!/z) * S/n! = sum S/(z*n!)
               total/(m!*n!) = sum S*m!/(z*m!*n!) = sum S/(z*n!)
               These are the same. So a(m,n) = total / (m! * n!).
            */
            mpz_mul(fact_n_sq, fact_temp, fact_n);  /* m! * n! */
            mpz_divexact(result, total, fact_n_sq);

            clock_gettime(CLOCK_MONOTONIC, &t_end);
            double dt = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

            gmp_printf("%d %d %Zd\n", m, n, result);
            if (m == n)
                fprintf(stderr, "a(%d,%d): %zu digits, %.3fs\n", m, n, mpz_sizeinbase(result, 10), dt);
            fflush(stdout); fflush(stderr);

            mpz_clear(fact_n);
            if (dt > 300) goto done;
        }

        /* Clean up partitions */
        for (int i = 0; i < num_parts; i++) {
            mpz_clear(parts[i].z_val);
            mpz_clear(parts[i].pcount);
        }
    }

done:
    mpz_clear(total); mpz_clear(fact_n_sq); mpz_clear(result); mpz_clear(temp);
    mpz_clear(fact_temp);
    for (int i = 0; i < MAX_N; i++) { mpz_clear(S_arr[i]); mpz_clear(c_arr[i]); }
    return 0;
}
