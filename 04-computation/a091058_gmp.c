/*
 * Compute A091058: n x n matrices over {1,...,n} under row/col/symbol perm.
 *
 * a(n) = (1/(n!)^3) * sum_{s,u partitions of n}
 *            |C(s)|*|C(u)| * sum_t |C(t)| * prod_{i,j} sigma_u(lcm(i,j))^{gcd(i,j)*s_i*t_j}
 *
 * Using GF trick on t (columns):
 * For each (s, u) pair, compute w_j = prod_i sigma_u(lcm(i,j))^{gcd(i,j)*s_i}
 * then [x^n] exp(sum_j w_j*x^j/j) gives the t-sum efficiently.
 *
 * sigma_u(L) = sum_{d|L} d * u_d  (sum of d*multiplicity over parts d dividing L)
 *
 * OEIS: A091058, had 15 terms (n=0..14)
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o a091058_gmp a091058_gmp.c -lgmp
 *
 * Usage: ./a091058_gmp <max_n>
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

static int lcm_func(int a, int b) {
    return a / gcd_func(a, b) * b;
}

#define MAX_N 256
#define MAX_PARTITIONS 200000

typedef struct {
    int pk[64];  /* distinct part values */
    int pm[64];  /* multiplicities */
    int depth;
    mpz_t z_val;
    mpz_t pcount;  /* n! / z */
} part_t;

static part_t *row_parts;  /* partitions for s (rows) */
static part_t *sym_parts;  /* partitions for u (symbols) */
static int num_row, num_sym;

static int pk_temp[64], pm_temp[64];
static mpz_t fact_n;

static void store_part(int depth, part_t *store, int *count) {
    int idx = (*count)++;
    for (int i = 0; i < depth; i++) {
        store[idx].pk[i] = pk_temp[i];
        store[idx].pm[i] = pm_temp[i];
    }
    store[idx].depth = depth;
    mpz_init(store[idx].z_val);
    mpz_set_ui(store[idx].z_val, 1);
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < pm_temp[i]; j++)
            mpz_mul_ui(store[idx].z_val, store[idx].z_val, pk_temp[i]);
        for (int j = 2; j <= pm_temp[i]; j++)
            mpz_mul_ui(store[idx].z_val, store[idx].z_val, j);
    }
    mpz_init(store[idx].pcount);
    mpz_divexact(store[idx].pcount, fact_n, store[idx].z_val);
}

static void enumerate(int remaining, int max_part, int depth, part_t *store, int *count) {
    if (remaining == 0) { store_part(depth, store, count); return; }
    for (int p = (max_part < remaining ? max_part : remaining); p >= 1; p--) {
        int max_m = remaining / p;
        for (int m = 1; m <= max_m; m++) {
            pk_temp[depth] = p; pm_temp[depth] = m;
            enumerate(remaining - m * p, p - 1, depth + 1, store, count);
        }
    }
}

/* Compute sigma_u(L) = sum_{d | L, d is part of u} d * u_d */
static long sigma_u(part_t *u, int L) {
    long s = 0;
    for (int i = 0; i < u->depth; i++) {
        if (L % u->pk[i] == 0) {
            s += (long)u->pk[i] * u->pm[i];
        }
    }
    return s;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <max_n>\n", argv[0]);
        return 1;
    }

    int max_n = atoi(argv[1]);
    fprintf(stderr, "Computing A091058 (n x n matrices over n symbols), n=0..%d\n\n", max_n);

    printf("0 1\n");
    fflush(stdout);

    row_parts = calloc(MAX_PARTITIONS, sizeof(part_t));
    sym_parts = calloc(MAX_PARTITIONS, sizeof(part_t));

    mpz_t total, fact_n_cube, result, temp;
    mpz_t S_arr[MAX_N], c_arr[MAX_N], w_arr[MAX_N];
    mpz_init(total); mpz_init(fact_n_cube); mpz_init(result); mpz_init(temp);
    mpz_init(fact_n);
    for (int i = 0; i < MAX_N; i++) {
        mpz_init(S_arr[i]); mpz_init(c_arr[i]); mpz_init(w_arr[i]);
    }

    for (int N = 1; N <= max_n; N++) {
        struct timespec t_start, t_end;
        clock_gettime(CLOCK_MONOTONIC, &t_start);

        mpz_fac_ui(fact_n, N);
        mpz_mul(fact_n_cube, fact_n, fact_n);
        mpz_mul(fact_n_cube, fact_n_cube, fact_n);

        /* Enumerate partitions */
        num_row = 0; num_sym = 0;
        enumerate(N, N, 0, row_parts, &num_row);
        /* sym_parts are the same set */
        for (int i = 0; i < num_row; i++) {
            sym_parts[i] = row_parts[i];
            /* Don't re-init GMP fields; share them via pointer copy would be dangerous.
               Instead, just copy the data and re-point. */
        }
        num_sym = num_row;
        /* Actually, sym_parts shares the same partition data as row_parts.
           We can iterate sym over the same array. */

        mpz_set_ui(total, 0);

        for (int si = 0; si < num_row; si++) {
            part_t *s = &row_parts[si];
            for (int ui = 0; ui < num_row; ui++) {
                part_t *u = &row_parts[ui];

                /* For each j = 1..N, compute w_j = prod_i sigma_u(lcm(i,j))^{gcd(i,j)*s_i} */
                for (int j = 1; j <= N; j++) {
                    mpz_set_ui(w_arr[j], 1);
                    int w_is_zero = 0;
                    for (int ii = 0; ii < s->depth && !w_is_zero; ii++) {
                        int ci = s->pk[ii];  /* row cycle length */
                        int mi = s->pm[ii];  /* multiplicity */
                        int g = gcd_func(ci, j);
                        int L = lcm_func(ci, j);
                        long sig = sigma_u(u, L);
                        int exp_val = g * mi;
                        if (sig == 0 && exp_val > 0) {
                            /* 0^positive = 0, w_j = 0 for this j */
                            mpz_set_ui(w_arr[j], 0);
                            w_is_zero = 1;
                        } else if (sig > 0 && exp_val > 0) {
                            mpz_set_ui(temp, sig);
                            mpz_pow_ui(temp, temp, exp_val);
                            mpz_mul(w_arr[j], w_arr[j], temp);
                        }
                    }
                }

                /* S recurrence: S(0) = 1, S(k) = sum_{j=1}^{k} w_j * S(k-j) * (k-1)!/(k-j)! */
                mpz_set_ui(S_arr[0], 1);
                for (int k = 1; k <= N; k++) {
                    mpz_set_ui(S_arr[k], 0);
                    mpz_t falling;
                    mpz_init(falling);
                    mpz_set_ui(falling, 1);
                    for (int j = 1; j <= k; j++) {
                        if (j >= 2) mpz_mul_ui(falling, falling, k - j + 1);
                        mpz_mul(temp, w_arr[j], S_arr[k - j]);
                        mpz_mul(temp, temp, falling);
                        mpz_add(S_arr[k], S_arr[k], temp);
                    }
                    mpz_clear(falling);
                }

                /* total += permcount(s) * permcount(u) * S(N) */
                mpz_mul(temp, s->pcount, u->pcount);
                mpz_mul(temp, temp, S_arr[N]);
                mpz_add(total, total, temp);
            }
        }

        mpz_divexact(result, total, fact_n_cube);

        clock_gettime(CLOCK_MONOTONIC, &t_end);
        double dt = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

        gmp_printf("%d %Zd\n", N, result);
        fprintf(stderr, "n=%d: %zu digits, %.3fs (%d partitions)\n",
                N, mpz_sizeinbase(result, 10), dt, num_row);
        fflush(stdout); fflush(stderr);

        /* Cleanup partition GMP fields */
        for (int i = 0; i < num_row; i++) {
            mpz_clear(row_parts[i].z_val);
            mpz_clear(row_parts[i].pcount);
        }

        if (dt > 3600) { fprintf(stderr, "Stopping\n"); break; }
    }

    mpz_clear(total); mpz_clear(fact_n_cube); mpz_clear(result); mpz_clear(temp);
    mpz_clear(fact_n);
    for (int i = 0; i < MAX_N; i++) {
        mpz_clear(S_arr[i]); mpz_clear(c_arr[i]); mpz_clear(w_arr[i]);
    }
    free(row_parts); free(sym_parts);
    return 0;
}
