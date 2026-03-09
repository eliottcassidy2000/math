/*
 * Compute A242095(n, k): n x n matrices over k symbols under S_n x S_n x S_k.
 * OpenMP-parallelized version.
 *
 * Same formula as a242095_gmp.c, parallelized over (s, u) pairs.
 *
 * Compile (macOS with Homebrew libomp):
 *   gcc -O3 -Xpreprocessor -fopenmp \
 *       -I/opt/homebrew/include -I/opt/homebrew/opt/libomp/include \
 *       -L/opt/homebrew/lib -L/opt/homebrew/opt/libomp/lib \
 *       -o a242095_omp a242095_omp.c -lgmp -lomp
 *
 * Usage: ./a242095_omp <max_n> <k> [num_threads]
 *   k=0: compute diagonal (k=n for each n) = A091058
 *
 * Author: opus-2026-03-08-S51
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <omp.h>

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
    int pk[64];
    int pm[64];
    int depth;
    mpz_t z_val;
    mpz_t pcount;
} part_t;

static long sigma_u(part_t *u, int L) {
    long s = 0;
    for (int i = 0; i < u->depth; i++) {
        if (L % u->pk[i] == 0)
            s += (long)u->pk[i] * u->pm[i];
    }
    return s;
}

static int pk_temp[64], pm_temp[64];

static void store_part(int depth, mpz_t fact, part_t *store, int *count) {
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
    mpz_divexact(store[idx].pcount, fact, store[idx].z_val);
}

static mpz_t enum_fact;

static void enumerate(int remaining, int max_part, int depth, part_t *store, int *count) {
    if (remaining == 0) { store_part(depth, enum_fact, store, count); return; }
    for (int p = (max_part < remaining ? max_part : remaining); p >= 1; p--) {
        int max_m = remaining / p;
        for (int m = 1; m <= max_m; m++) {
            pk_temp[depth] = p; pm_temp[depth] = m;
            enumerate(remaining - m * p, p - 1, depth + 1, store, count);
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <max_n> <k> [num_threads]\n", argv[0]);
        return 1;
    }

    int max_n = atoi(argv[1]);
    int fixed_k = atoi(argv[2]);
    int diagonal = (fixed_k == 0);
    int num_threads = (argc >= 4) ? atoi(argv[3]) : omp_get_max_threads();
    omp_set_num_threads(num_threads);

    if (diagonal)
        fprintf(stderr, "Computing A091058 (OpenMP, %d threads), n=0..%d\n\n", num_threads, max_n);
    else
        fprintf(stderr, "Computing A242095(n, %d) (OpenMP, %d threads), n=0..%d\n\n", fixed_k, num_threads, max_n);

    part_t *row_parts = calloc(MAX_PARTITIONS, sizeof(part_t));
    part_t *sym_parts = calloc(MAX_PARTITIONS, sizeof(part_t));
    mpz_init(enum_fact);

    mpz_t fact_n, fact_k, norm, result;
    mpz_init(fact_n); mpz_init(fact_k); mpz_init(norm); mpz_init(result);

    printf("0 1\n");
    fflush(stdout);

    for (int N = 1; N <= max_n; N++) {
        int K = diagonal ? N : fixed_k;

        struct timespec t_start, t_end;
        clock_gettime(CLOCK_MONOTONIC, &t_start);

        mpz_fac_ui(fact_n, N);
        mpz_fac_ui(fact_k, K);

        int num_row = 0, num_sym = 0;
        mpz_set(enum_fact, fact_n);
        enumerate(N, N, 0, row_parts, &num_row);
        mpz_set(enum_fact, fact_k);
        enumerate(K, K, 0, sym_parts, &num_sym);

        /* norm = (n!)^2 * k! */
        mpz_mul(norm, fact_n, fact_n);
        mpz_mul(norm, norm, fact_k);

        int total_pairs = num_row * num_sym;

        mpz_t *thread_totals = malloc(num_threads * sizeof(mpz_t));
        for (int t = 0; t < num_threads; t++)
            mpz_init_set_ui(thread_totals[t], 0);

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            mpz_t w_arr[MAX_N], S_arr[MAX_N], temp, falling;
            for (int i = 0; i < MAX_N; i++) {
                mpz_init(w_arr[i]);
                mpz_init(S_arr[i]);
            }
            mpz_init(temp);
            mpz_init(falling);

            #pragma omp for schedule(dynamic, 1)
            for (int pair = 0; pair < total_pairs; pair++) {
                int si = pair / num_sym;
                int ui = pair % num_sym;
                part_t *s = &row_parts[si];
                part_t *u = &sym_parts[ui];

                for (int j = 1; j <= N; j++) {
                    mpz_set_ui(w_arr[j], 1);
                    int w_is_zero = 0;
                    for (int ii = 0; ii < s->depth && !w_is_zero; ii++) {
                        int ci = s->pk[ii];
                        int mi = s->pm[ii];
                        int g = gcd_func(ci, j);
                        int L = lcm_func(ci, j);
                        long sig = sigma_u(u, L);
                        int exp_val = g * mi;
                        if (sig == 0 && exp_val > 0) {
                            mpz_set_ui(w_arr[j], 0);
                            w_is_zero = 1;
                        } else if (sig > 0 && exp_val > 0) {
                            mpz_set_ui(temp, sig);
                            mpz_pow_ui(temp, temp, exp_val);
                            mpz_mul(w_arr[j], w_arr[j], temp);
                        }
                    }
                }

                mpz_set_ui(S_arr[0], 1);
                for (int k = 1; k <= N; k++) {
                    mpz_set_ui(S_arr[k], 0);
                    mpz_set_ui(falling, 1);
                    for (int j = 1; j <= k; j++) {
                        if (j >= 2) mpz_mul_ui(falling, falling, k - j + 1);
                        mpz_mul(temp, w_arr[j], S_arr[k - j]);
                        mpz_mul(temp, temp, falling);
                        mpz_add(S_arr[k], S_arr[k], temp);
                    }
                }

                mpz_mul(temp, s->pcount, u->pcount);
                mpz_mul(temp, temp, S_arr[N]);
                mpz_add(thread_totals[tid], thread_totals[tid], temp);
            }

            for (int i = 0; i < MAX_N; i++) {
                mpz_clear(w_arr[i]);
                mpz_clear(S_arr[i]);
            }
            mpz_clear(temp);
            mpz_clear(falling);
        }

        mpz_set_ui(result, 0);
        for (int t = 0; t < num_threads; t++) {
            mpz_add(result, result, thread_totals[t]);
            mpz_clear(thread_totals[t]);
        }
        free(thread_totals);

        mpz_divexact(result, result, norm);

        clock_gettime(CLOCK_MONOTONIC, &t_end);
        double dt = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

        gmp_printf("%d %Zd\n", N, result);
        fprintf(stderr, "n=%d: %zu digits, %.3fs (%d row, %d sym, %d pairs)\n",
                N, mpz_sizeinbase(result, 10), dt, num_row, num_sym, total_pairs);
        fflush(stdout); fflush(stderr);

        for (int i = 0; i < num_row; i++) {
            mpz_clear(row_parts[i].z_val);
            mpz_clear(row_parts[i].pcount);
        }
        for (int i = 0; i < num_sym; i++) {
            mpz_clear(sym_parts[i].z_val);
            mpz_clear(sym_parts[i].pcount);
        }

        if (dt > 7200) { fprintf(stderr, "Stopping\n"); break; }
    }

    mpz_clear(fact_n); mpz_clear(fact_k); mpz_clear(norm); mpz_clear(result);
    mpz_clear(enum_fact);
    free(row_parts); free(sym_parts);
    return 0;
}
