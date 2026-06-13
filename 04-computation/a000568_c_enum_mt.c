/*
 * a000568_c_enum_mt.c — A000568 mod p via multi-threaded partition enumeration.
 *
 * Like a000568_c_enum.c but uses pthreads for parallel enumeration.
 * Pure int64 arithmetic (no GMP) — each thread accumulates mod p.
 *
 * Compile:
 *   gcc -O3 -o a000568_c_enum_mt a000568_c_enum_mt.c -pthread
 *
 * Usage:
 *   ./a000568_c_enum_mt <n> <prime> [num_threads]
 *
 * Author: opus-2026-03-08-S48
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <stdatomic.h>

static int64_t MOD_P;
static int global_n;

static int64_t mod_pow(int64_t base, int64_t exp, int64_t mod) {
    int64_t result = 1;
    base %= mod;
    if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1)
            result = (__int128)result * base % mod;
        base = (__int128)base * base % mod;
        exp >>= 1;
    }
    return result;
}

static int64_t mod_inv(int64_t x, int64_t p) {
    return mod_pow(x, p - 2, p);
}

/* GCD table */
static int gcd_tab[502][502];
static void init_gcd(int n) {
    for (int i = 1; i <= n; i += 2)
        for (int j = i; j <= n; j += 2) {
            int a = i, b = j;
            while (b) { int t = b; b = a % b; a = t; }
            gcd_tab[i][j] = gcd_tab[j][i] = a;
        }
}

static int odd_parts[256];
static int num_odd;
static int64_t fact_mod[502];
static int64_t fact_inv_mod[502];

/* Precomputed: inv_k_pow[pi][m] = (k^{-1})^m mod p = k^{-m} mod p */
/* inv_k_pow[pi][0] = 1 */
static int64_t *inv_k_pow[256];
static int inv_k_max_m[256];

/* Work queue */
typedef struct { int pi; int m; } WorkItem;
static WorkItem *work_items;
static int num_work_items;
static atomic_int work_counter;

/* Per-thread data */
typedef struct { int k; int m; } PartEntry;

typedef struct {
    int64_t local_sum;
    PartEntry stack[256];
    int stack_size;
} ThreadData;

static void thread_enumerate(ThreadData *td, int remaining, int max_part_idx,
                              int64_t t_val, int64_t z_inv) {
    if (remaining == 0) {
        int64_t pow2 = mod_pow(2, t_val, MOD_P);
        int64_t contrib = (__int128)pow2 * z_inv % MOD_P;
        td->local_sum = (td->local_sum + contrib) % MOD_P;
        return;
    }

    for (int pi = max_part_idx; pi >= 0; pi--) {
        int k = odd_parts[pi];
        if (k > remaining) continue;
        int max_m = remaining / k;

        for (int m = 1; m <= max_m; m++) {
            int64_t dt_self = (int64_t)m * (m - 1) * k / 2 + (int64_t)m * (k - 1) / 2;
            int64_t dt_cross = 0;
            for (int s = 0; s < td->stack_size; s++)
                dt_cross += (int64_t)m * td->stack[s].m * gcd_tab[k][td->stack[s].k];

            /* z_inv update: z_inv *= k^{-m} * (m!)^{-1} */
            int64_t new_z_inv = (__int128)z_inv * inv_k_pow[pi][m] % MOD_P;
            new_z_inv = (__int128)new_z_inv * fact_inv_mod[m] % MOD_P;

            td->stack[td->stack_size].k = k;
            td->stack[td->stack_size].m = m;
            td->stack_size++;

            thread_enumerate(td, remaining - m * k, pi - 1,
                            t_val + dt_self + dt_cross, new_z_inv);

            td->stack_size--;
        }
    }
}

static void *thread_worker(void *arg) {
    ThreadData *td = (ThreadData *)arg;
    int n = global_n;

    while (1) {
        int idx = atomic_fetch_add(&work_counter, 1);
        if (idx >= num_work_items) break;

        int pi = work_items[idx].pi;
        int m = work_items[idx].m;
        int k = odd_parts[pi];

        int64_t dt_self = (int64_t)m * (m - 1) * k / 2 + (int64_t)m * (k - 1) / 2;
        int64_t z_inv = (__int128)inv_k_pow[pi][m] * fact_inv_mod[m] % MOD_P;

        td->stack[0].k = k;
        td->stack[0].m = m;
        td->stack_size = 1;

        thread_enumerate(td, n - m * k, pi - 1, dt_self, z_inv);
    }

    return NULL;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <n> <prime> [num_threads]\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    int64_t p = atoll(argv[2]);
    int num_threads = 1;
    if (argc >= 4) num_threads = atoi(argv[3]);

    MOD_P = p;
    global_n = n;

    if (n <= 1) {
        printf("%lld\n", (long long)(1 % p));
        return 0;
    }

    init_gcd(n);

    num_odd = 0;
    for (int k = 1; k <= n; k += 2)
        odd_parts[num_odd++] = k;

    /* Precompute factorials and inverses */
    fact_mod[0] = 1;
    for (int i = 1; i <= n; i++)
        fact_mod[i] = (__int128)fact_mod[i-1] * i % p;
    fact_inv_mod[n] = mod_inv(fact_mod[n], p);
    for (int i = n - 1; i >= 0; i--)
        fact_inv_mod[i] = (__int128)fact_inv_mod[i+1] * (i+1) % p;

    /* Precompute inv_k_pow[pi][m] = k^{-m} mod p */
    for (int pi = 0; pi < num_odd; pi++) {
        int k = odd_parts[pi];
        int max_m = n / k;
        inv_k_max_m[pi] = max_m;
        inv_k_pow[pi] = malloc((max_m + 1) * sizeof(int64_t));
        int64_t kinv = mod_inv(k, p);
        inv_k_pow[pi][0] = 1;
        for (int m = 1; m <= max_m; m++)
            inv_k_pow[pi][m] = (__int128)inv_k_pow[pi][m-1] * kinv % p;
    }

    /* Build work queue */
    num_work_items = 0;
    for (int pi = num_odd - 1; pi >= 0; pi--) {
        int k = odd_parts[pi];
        if (k > n) continue;
        num_work_items += n / k;
    }
    work_items = malloc(num_work_items * sizeof(WorkItem));
    int wi = 0;
    for (int pi = num_odd - 1; pi >= 0; pi--) {
        int k = odd_parts[pi];
        if (k > n) continue;
        for (int m = 1; m <= n / k; m++) {
            work_items[wi].pi = pi;
            work_items[wi].m = m;
            wi++;
        }
    }

    if (num_threads > num_work_items) num_threads = num_work_items;
    if (num_threads < 1) num_threads = 1;
    atomic_store(&work_counter, 0);

    ThreadData *threads = calloc(num_threads, sizeof(ThreadData));
    pthread_t *pthreads = calloc(num_threads, sizeof(pthread_t));

    for (int t = 0; t < num_threads; t++) {
        threads[t].local_sum = 0;
        threads[t].stack_size = 0;
    }

    for (int t = 0; t < num_threads; t++)
        pthread_create(&pthreads[t], NULL, thread_worker, &threads[t]);
    for (int t = 0; t < num_threads; t++)
        pthread_join(pthreads[t], NULL);

    int64_t result_sum = 0;
    for (int t = 0; t < num_threads; t++)
        result_sum = (result_sum + threads[t].local_sum) % p;

    printf("%lld\n", (long long)result_sum);

    /* Cleanup */
    for (int pi = 0; pi < num_odd; pi++)
        free(inv_k_pow[pi]);
    free(threads);
    free(pthreads);
    free(work_items);

    return 0;
}
