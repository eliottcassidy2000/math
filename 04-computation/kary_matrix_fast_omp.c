/*
 * Compute k-ary n×n matrices under row/col perm (OpenMP-parallelized).
 * Column k of A246106. k=2: A002724, k=3: A052269, etc.
 *
 * a(n) = (1/n!) * sum over partitions q of n:
 *    permcount(q) * [x^n] exp(sum_{t=1}^{n} k^{K(q,t)} * x^t / t)
 *
 * K(q,t) = sum_i gcd(t, q_i) * m_i
 *
 * Parallelized over partitions.
 *
 * Compile (macOS):
 *   gcc -O3 -Xpreprocessor -fopenmp \
 *       -I/opt/homebrew/include -I/opt/homebrew/opt/libomp/include \
 *       -L/opt/homebrew/lib -L/opt/homebrew/opt/libomp/lib \
 *       -o kary_matrix_fast_omp kary_matrix_fast_omp.c -lgmp -lomp
 *
 * Usage: ./kary_matrix_fast_omp <k> <max_n> [num_threads]
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

#define MAX_N 512
#define MAX_PARTITIONS 200000

typedef struct {
    int pk[128];
    int pm[128];
    int depth;
    mpz_t z_val;
    mpz_t pcount;
} part_t;

static part_t *parts;
static int num_parts;
static int pk_temp[128], pm_temp[128];
static mpz_t fact_n;

static void store_part(int depth) {
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
    mpz_divexact(parts[idx].pcount, fact_n, parts[idx].z_val);
}

static void enumerate(int remaining, int max_part, int depth) {
    if (remaining == 0) { store_part(depth); return; }
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
        fprintf(stderr, "Usage: %s <k> <max_n> [num_threads]\n", argv[0]);
        return 1;
    }

    int K = atoi(argv[1]);
    int max_n = atoi(argv[2]);
    int num_threads = (argc >= 4) ? atoi(argv[3]) : omp_get_max_threads();
    omp_set_num_threads(num_threads);

    fprintf(stderr, "Computing %d-ary n×n matrices (OMP, %d threads), n=0..%d\n\n",
            K, num_threads, max_n);

    parts = calloc(MAX_PARTITIONS, sizeof(part_t));
    mpz_init(fact_n);

    printf("0 1\n");
    fflush(stdout);

    for (int N = 1; N <= max_n; N++) {
        struct timespec t_start, t_end;
        clock_gettime(CLOCK_MONOTONIC, &t_start);

        mpz_fac_ui(fact_n, N);
        num_parts = 0;
        enumerate(N, N, 0);

        mpz_t *thread_totals = malloc(num_threads * sizeof(mpz_t));
        for (int t = 0; t < num_threads; t++)
            mpz_init_set_ui(thread_totals[t], 0);

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            mpz_t S_arr[MAX_N], c_arr[MAX_N], temp, falling;
            for (int i = 0; i < MAX_N; i++) {
                mpz_init(S_arr[i]);
                mpz_init(c_arr[i]);
            }
            mpz_init(temp);
            mpz_init(falling);

            #pragma omp for schedule(dynamic, 1)
            for (int pi = 0; pi < num_parts; pi++) {
                part_t *p = &parts[pi];

                for (int t = 1; t <= N; t++) {
                    int Kval = 0;
                    for (int i = 0; i < p->depth; i++)
                        Kval += p->pm[i] * gcd_func(t, p->pk[i]);
                    mpz_set_ui(c_arr[t], K);
                    mpz_pow_ui(c_arr[t], c_arr[t], Kval);
                }

                mpz_set_ui(S_arr[0], 1);
                for (int k = 1; k <= N; k++) {
                    mpz_set_ui(S_arr[k], 0);
                    mpz_set_ui(falling, 1);
                    for (int j = 1; j <= k; j++) {
                        if (j >= 2) mpz_mul_ui(falling, falling, k - j + 1);
                        mpz_mul(temp, c_arr[j], S_arr[k - j]);
                        mpz_mul(temp, temp, falling);
                        mpz_add(S_arr[k], S_arr[k], temp);
                    }
                }

                mpz_mul(temp, p->pcount, S_arr[N]);
                mpz_add(thread_totals[tid], thread_totals[tid], temp);
            }

            for (int i = 0; i < MAX_N; i++) {
                mpz_clear(S_arr[i]);
                mpz_clear(c_arr[i]);
            }
            mpz_clear(temp);
            mpz_clear(falling);
        }

        mpz_t result;
        mpz_init_set_ui(result, 0);
        for (int t = 0; t < num_threads; t++) {
            mpz_add(result, result, thread_totals[t]);
            mpz_clear(thread_totals[t]);
        }
        free(thread_totals);

        mpz_t fact_n_sq;
        mpz_init(fact_n_sq);
        mpz_mul(fact_n_sq, fact_n, fact_n);
        mpz_divexact(result, result, fact_n_sq);

        clock_gettime(CLOCK_MONOTONIC, &t_end);
        double dt = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

        gmp_printf("%d %Zd\n", N, result);
        fprintf(stderr, "n=%d: %zu digits, %.3fs (%d partitions)\n",
                N, mpz_sizeinbase(result, 10), dt, num_parts);
        fflush(stdout); fflush(stderr);

        mpz_clear(result);
        mpz_clear(fact_n_sq);

        for (int i = 0; i < num_parts; i++) {
            mpz_clear(parts[i].z_val);
            mpz_clear(parts[i].pcount);
        }

        if (dt > 3600) { fprintf(stderr, "Stopping\n"); break; }
    }

    mpz_clear(fact_n);
    free(parts);
    return 0;
}
