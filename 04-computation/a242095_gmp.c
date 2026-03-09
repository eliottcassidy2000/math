/*
 * Compute A242095(n, k): n x n matrices over k symbols under S_n x S_n x S_k.
 *
 * a(n,k) = sum_{s part of n, u part of k, t part of n}
 *          fixA(s,t,u) / (z(s)*z(t)*z(u))
 *
 * fixA = prod_{i,j} (sum_{d|lcm(i,j)} d*u_d)^{gcd(i,j)*s_i*t_j}
 *
 * GF trick on t: for each (s,u), compute w_j and use exp GF.
 *
 * Diagonal (k=n): A091058
 * Column k=2: A091059, k=3: A091060, k=4: A091061, k=5: A091062
 * k=6..10: A246122..A246126
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o a242095_gmp a242095_gmp.c -lgmp
 *
 * Usage: ./a242095_gmp <max_n> <k>
 *   If k=0, computes diagonal (k=n for each n)
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
    int pk[64];
    int pm[64];
    int depth;
    mpz_t z_val;
    mpz_t pcount;
} part_t;

static part_t *row_parts;
static part_t *sym_parts;
static int num_row, num_sym;
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
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <max_n> <k>\n", argv[0]);
        fprintf(stderr, "  k=0: compute diagonal (k=n)\n");
        fprintf(stderr, "  k=2: A091059, k=3: A091060, ...\n");
        return 1;
    }

    int max_n = atoi(argv[1]);
    int fixed_k = atoi(argv[2]);
    int diagonal = (fixed_k == 0);

    if (diagonal)
        fprintf(stderr, "Computing A091058 (diagonal of A242095), n=0..%d\n\n", max_n);
    else
        fprintf(stderr, "Computing A242095(n, %d), n=0..%d\n\n", fixed_k, max_n);

    row_parts = calloc(MAX_PARTITIONS, sizeof(part_t));
    sym_parts = calloc(MAX_PARTITIONS, sizeof(part_t));
    mpz_init(enum_fact);

    mpz_t total, fact_n, fact_k, norm, result, temp;
    mpz_t S_arr[MAX_N], w_arr[MAX_N];
    mpz_init(total); mpz_init(fact_n); mpz_init(fact_k);
    mpz_init(norm); mpz_init(result); mpz_init(temp);
    for (int i = 0; i < MAX_N; i++) { mpz_init(S_arr[i]); mpz_init(w_arr[i]); }

    printf("0 1\n");
    fflush(stdout);

    for (int N = 1; N <= max_n; N++) {
        int K = diagonal ? N : fixed_k;

        struct timespec t_start, t_end;
        clock_gettime(CLOCK_MONOTONIC, &t_start);

        mpz_fac_ui(fact_n, N);
        mpz_fac_ui(fact_k, K);

        /* Enumerate row partitions of N */
        num_row = 0;
        mpz_set(enum_fact, fact_n);
        enumerate(N, N, 0, row_parts, &num_row);

        /* Enumerate symbol partitions of K */
        num_sym = 0;
        mpz_set(enum_fact, fact_k);
        enumerate(K, K, 0, sym_parts, &num_sym);

        /* norm = (n!)^2 * k! */
        mpz_mul(norm, fact_n, fact_n);
        mpz_mul(norm, norm, fact_k);

        mpz_set_ui(total, 0);

        for (int si = 0; si < num_row; si++) {
            part_t *s = &row_parts[si];
            for (int ui = 0; ui < num_sym; ui++) {
                part_t *u = &sym_parts[ui];

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

        mpz_divexact(result, total, norm);

        clock_gettime(CLOCK_MONOTONIC, &t_end);
        double dt = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

        gmp_printf("%d %Zd\n", N, result);
        fprintf(stderr, "n=%d: %zu digits, %.3fs (%d row parts, %d sym parts)\n",
                N, mpz_sizeinbase(result, 10), dt, num_row, num_sym);
        fflush(stdout); fflush(stderr);

        /* Cleanup */
        for (int i = 0; i < num_row; i++) {
            mpz_clear(row_parts[i].z_val);
            mpz_clear(row_parts[i].pcount);
        }
        for (int i = 0; i < num_sym; i++) {
            mpz_clear(sym_parts[i].z_val);
            mpz_clear(sym_parts[i].pcount);
        }

        if (dt > 3600) { fprintf(stderr, "Stopping\n"); break; }
    }

    mpz_clear(total); mpz_clear(fact_n); mpz_clear(fact_k);
    mpz_clear(norm); mpz_clear(result); mpz_clear(temp);
    mpz_clear(enum_fact);
    for (int i = 0; i < MAX_N; i++) { mpz_clear(S_arr[i]); mpz_clear(w_arr[i]); }
    free(row_parts); free(sym_parts);
    return 0;
}
