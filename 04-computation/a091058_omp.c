/*
 * Compute A091058: n x n matrices over {1,...,n} under row/col/symbol perm.
 * OpenMP-parallelized version — ~4-8x faster than serial on multi-core.
 *
 * Same Burnside/GF approach as a091058_gmp.c but with the (s,u) pair loop
 * parallelized using OpenMP. Each thread has its own w_arr, S_arr, temp.
 *
 * Compile (macOS with Homebrew libomp):
 *   gcc -O3 -Xpreprocessor -fopenmp \
 *       -I/opt/homebrew/include -I/opt/homebrew/opt/libomp/include \
 *       -L/opt/homebrew/lib -L/opt/homebrew/opt/libomp/lib \
 *       -o a091058_omp a091058_omp.c -lgmp -lomp
 *
 * Usage: ./a091058_omp <max_n> [num_threads]
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

static part_t *parts;
static int num_parts;
static int pk_temp[64], pm_temp[64];
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
        fprintf(stderr, "Usage: %s <max_n> [num_threads]\n", argv[0]);
        return 1;
    }

    int max_n = atoi(argv[1]);
    int num_threads = (argc >= 3) ? atoi(argv[2]) : omp_get_max_threads();
    omp_set_num_threads(num_threads);

    fprintf(stderr, "Computing A091058 (OpenMP, %d threads), n=0..%d\n\n", num_threads, max_n);

    printf("0 1\n");
    fflush(stdout);

    parts = calloc(MAX_PARTITIONS, sizeof(part_t));
    mpz_init(fact_n);

    mpz_t fact_n_cube, result;
    mpz_init(fact_n_cube); mpz_init(result);

    for (int N = 1; N <= max_n; N++) {
        struct timespec t_start, t_end;
        clock_gettime(CLOCK_MONOTONIC, &t_start);

        mpz_fac_ui(fact_n, N);
        mpz_mul(fact_n_cube, fact_n, fact_n);
        mpz_mul(fact_n_cube, fact_n_cube, fact_n);

        num_parts = 0;
        enumerate(N, N, 0);

        int total_pairs = num_parts * num_parts;

        /* Thread-local partial sums, allocated once per N */
        mpz_t *thread_totals = malloc(num_threads * sizeof(mpz_t));
        for (int t = 0; t < num_threads; t++)
            mpz_init_set_ui(thread_totals[t], 0);

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            /* Thread-local GMP variables */
            mpz_t w_arr[MAX_N], S_arr[MAX_N], temp, falling;
            for (int i = 0; i < MAX_N; i++) {
                mpz_init(w_arr[i]);
                mpz_init(S_arr[i]);
            }
            mpz_init(temp);
            mpz_init(falling);

            #pragma omp for schedule(dynamic, 1)
            for (int pair = 0; pair < total_pairs; pair++) {
                int si = pair / num_parts;
                int ui = pair % num_parts;
                part_t *s = &parts[si];
                part_t *u = &parts[ui];

                /* Compute w_j for j = 1..N */
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

                /* S recurrence */
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

                /* Accumulate: permcount(s) * permcount(u) * S(N) */
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

        /* Reduce thread totals */
        mpz_set_ui(result, 0);
        for (int t = 0; t < num_threads; t++) {
            mpz_add(result, result, thread_totals[t]);
            mpz_clear(thread_totals[t]);
        }
        free(thread_totals);

        mpz_divexact(result, result, fact_n_cube);

        clock_gettime(CLOCK_MONOTONIC, &t_end);
        double dt = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

        gmp_printf("%d %Zd\n", N, result);
        fprintf(stderr, "n=%d: %zu digits, %.3fs (%d partitions, %d pairs)\n",
                N, mpz_sizeinbase(result, 10), dt, num_parts, total_pairs);
        fflush(stdout); fflush(stderr);

        for (int i = 0; i < num_parts; i++) {
            mpz_clear(parts[i].z_val);
            mpz_clear(parts[i].pcount);
        }

        if (dt > 7200) { fprintf(stderr, "Stopping (>2h per term)\n"); break; }
    }

    mpz_clear(fact_n_cube); mpz_clear(result);
    mpz_clear(fact_n);
    free(parts);
    return 0;
}
