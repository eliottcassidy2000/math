/*
 * Compute FOUR sequences simultaneously (OpenMP-parallelized):
 *   A002724: binary n×n matrices under row/col perm
 *   A006383: binary n×n matrices under row/col perm + col complement
 *   A122082: bicolored graphs invariant under color swap
 *   A007139: bicolored bipartite graphs = (A002724 + A122082)/2
 *
 * Same algorithm as matrix_quad_gmp.c with OpenMP parallelism over partitions.
 *
 * Compile (macOS):
 *   gcc -O3 -Xpreprocessor -fopenmp \
 *       -I/opt/homebrew/include -I/opt/homebrew/opt/libomp/include \
 *       -L/opt/homebrew/lib -L/opt/homebrew/opt/libomp/lib \
 *       -o matrix_quad_omp matrix_quad_omp.c -lgmp -lomp
 *
 * Usage: ./matrix_quad_omp <max_n> [num_threads]
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

static int v2(int x) {
    if (x == 0) return 0;
    int v = 0;
    while (x % 2 == 0) { x /= 2; v++; }
    return v;
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
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <max_n> [num_threads]\n", argv[0]);
        return 1;
    }

    int max_n = atoi(argv[1]);
    int num_threads = (argc >= 3) ? atoi(argv[2]) : omp_get_max_threads();
    omp_set_num_threads(num_threads);

    fprintf(stderr, "Computing A002724, A006383, A122082, A007139 (OpenMP, %d threads), n=0..%d\n\n",
            num_threads, max_n);

    parts = calloc(MAX_PARTITIONS, sizeof(part_t));
    mpz_init(fact_n);

    printf("0 1 1 1 1\n");
    fflush(stdout);

    for (int N = 1; N <= max_n; N++) {
        struct timespec t_start, t_end;
        clock_gettime(CLOCK_MONOTONIC, &t_start);

        mpz_fac_ui(fact_n, N);
        num_parts = 0;
        enumerate(N, N, 0);

        mpz_t *t724 = malloc(num_threads * sizeof(mpz_t));
        mpz_t *t383 = malloc(num_threads * sizeof(mpz_t));
        mpz_t *t122 = malloc(num_threads * sizeof(mpz_t));
        for (int t = 0; t < num_threads; t++) {
            mpz_init_set_ui(t724[t], 0);
            mpz_init_set_ui(t383[t], 0);
            mpz_init_set_ui(t122[t], 0);
        }

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            mpz_t S724[MAX_N], S383[MAX_N], c724[MAX_N], c383[MAX_N];
            mpz_t temp, falling;
            for (int i = 0; i < MAX_N; i++) {
                mpz_init(S724[i]); mpz_init(S383[i]);
                mpz_init(c724[i]); mpz_init(c383[i]);
            }
            mpz_init(temp); mpz_init(falling);

            #pragma omp for schedule(dynamic, 1)
            for (int pi = 0; pi < num_parts; pi++) {
                part_t *p = &parts[pi];
                int depth = p->depth;

                /* A122082: single partition formula */
                int num_cycles = 0;
                for (int i = 0; i < depth; i++) num_cycles += p->pm[i];
                int cycles[512];
                int idx = 0;
                for (int i = 0; i < depth; i++)
                    for (int j = 0; j < p->pm[i]; j++)
                        cycles[idx++] = p->pk[i];

                int edges = 0;
                for (int i = 0; i < num_cycles; i++)
                    for (int j = i + 1; j < num_cycles; j++)
                        edges += gcd_func(cycles[i], cycles[j]);
                for (int i = 0; i < num_cycles; i++)
                    edges += (cycles[i] + 1) / 2;

                mpz_set_ui(temp, 1);
                mpz_mul_2exp(temp, temp, edges);
                mpz_mul(temp, temp, p->pcount);
                mpz_add(t122[tid], t122[tid], temp);

                /* A002724 and A006383: GF approach */
                int min_v2_val = 100;
                for (int i = 0; i < depth; i++) {
                    int vv = v2(p->pk[i]);
                    if (vv < min_v2_val) min_v2_val = vv;
                }
                int e = 1 << min_v2_val;

                for (int t = 1; t <= N; t++) {
                    int K = 0;
                    for (int i = 0; i < depth; i++)
                        K += p->pm[i] * gcd_func(t, p->pk[i]);
                    mpz_set_ui(c724[t], 1);
                    mpz_mul_2exp(c724[t], c724[t], K);
                    int flag = (t % e == 0) ? 1 : 0;
                    mpz_set_ui(c383[t], 1);
                    mpz_mul_2exp(c383[t], c383[t], K - flag);
                }

                mpz_set_ui(S724[0], 1);
                mpz_set_ui(S383[0], 1);
                for (int k = 1; k <= N; k++) {
                    mpz_set_ui(S724[k], 0);
                    mpz_set_ui(S383[k], 0);
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
                }

                mpz_mul(temp, p->pcount, S724[N]);
                mpz_add(t724[tid], t724[tid], temp);

                mpz_mul(temp, p->pcount, S383[N]);
                mpz_add(t383[tid], t383[tid], temp);
            }

            for (int i = 0; i < MAX_N; i++) {
                mpz_clear(S724[i]); mpz_clear(S383[i]);
                mpz_clear(c724[i]); mpz_clear(c383[i]);
            }
            mpz_clear(temp); mpz_clear(falling);
        }

        /* Reduce */
        mpz_t total_724, total_383, total_122;
        mpz_init_set_ui(total_724, 0);
        mpz_init_set_ui(total_383, 0);
        mpz_init_set_ui(total_122, 0);
        for (int t = 0; t < num_threads; t++) {
            mpz_add(total_724, total_724, t724[t]);
            mpz_add(total_383, total_383, t383[t]);
            mpz_add(total_122, total_122, t122[t]);
            mpz_clear(t724[t]); mpz_clear(t383[t]); mpz_clear(t122[t]);
        }
        free(t724); free(t383); free(t122);

        mpz_t fact_n_sq, r724, r383, r122, r139;
        mpz_init(fact_n_sq); mpz_init(r724); mpz_init(r383);
        mpz_init(r122); mpz_init(r139);

        mpz_mul(fact_n_sq, fact_n, fact_n);
        mpz_divexact(r724, total_724, fact_n_sq);
        mpz_divexact(r383, total_383, fact_n_sq);
        mpz_divexact(r122, total_122, fact_n);
        mpz_add(r139, r724, r122);
        mpz_divexact_ui(r139, r139, 2);

        clock_gettime(CLOCK_MONOTONIC, &t_end);
        double dt = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

        gmp_printf("%d %Zd %Zd %Zd %Zd\n", N, r724, r383, r122, r139);
        fprintf(stderr, "n=%d: %.3fs (%d partitions)\n", N, dt, num_parts);
        fflush(stdout); fflush(stderr);

        mpz_clear(total_724); mpz_clear(total_383); mpz_clear(total_122);
        mpz_clear(fact_n_sq); mpz_clear(r724); mpz_clear(r383);
        mpz_clear(r122); mpz_clear(r139);

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
