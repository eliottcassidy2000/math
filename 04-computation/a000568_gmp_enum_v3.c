/*
 * a000568_gmp_enum_v3.c — A000568: optimized GMP enum with precomputed factors.
 *
 * Optimizations over v1:
 *   1. Precomputed z-factor table: km_fac[k][m] = k^m * m! as GMP.
 *      Avoids calling mpz_ui_pow_ui at every leaf.
 *   2. Multi-threaded with atomic work-queue (from MT version).
 *   3. Avoids mpz_init/clear inside the hot loop (leaf processing).
 *
 * Compile:
 *   gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *       -o a000568_gmp_enum_v3 a000568_gmp_enum_v3.c -lgmp -pthread
 *
 * Author: opus-2026-03-08-S48
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <stdatomic.h>
#include <gmp.h>

#define MAX_N 1000
#define MAX_PARTS 500

/* Shared read-only data */
static int gcd_tab[MAX_N + 1][MAX_N + 1];
static int odd_parts[MAX_PARTS];
static int num_odd;
static mpz_t fact[MAX_N + 1];
static mpz_t LCD;
static int global_n;

/* Precomputed km_fac[odd_index][m] = k^m * m! */
static mpz_t *km_fac[MAX_PARTS]; /* dynamically allocated per odd part */
static int km_max_m[MAX_PARTS];

static int gcd_int(int a, int b) {
    while (b) { int t = b; b = a % b; a = t; }
    return a;
}

static void init_gcd(int n) {
    for (int i = 1; i <= n; i += 2)
        for (int j = i; j <= n; j += 2) {
            int g = gcd_int(i, j);
            gcd_tab[i][j] = gcd_tab[j][i] = g;
        }
}

/* Work queue */
typedef struct { int pi; int m; } WorkItem;
static WorkItem *work_items;
static int num_work_items;
static atomic_int work_counter;

/* Per-thread data */
typedef struct { int k_idx; int m; } PartEntry; /* k_idx into odd_parts[] */

typedef struct {
    mpz_t local_sum;
    long long local_count;
    PartEntry stack[MAX_PARTS];
    int stack_size;
    mpz_t tmp_z_val, tmp_contrib, tmp_lcd_over_z;
} ThreadData;

static void thread_enumerate(ThreadData *td, int remaining, int max_part_idx,
                              long long t_val)
{
    if (remaining == 0) {
        td->local_count++;

        /* Compute z_val from precomputed table */
        mpz_set_ui(td->tmp_z_val, 1);
        for (int i = 0; i < td->stack_size; i++) {
            int ki = td->stack[i].k_idx;
            int m = td->stack[i].m;
            mpz_mul(td->tmp_z_val, td->tmp_z_val, km_fac[ki][m]);
        }

        mpz_tdiv_q(td->tmp_lcd_over_z, LCD, td->tmp_z_val);
        mpz_mul_2exp(td->tmp_contrib, td->tmp_lcd_over_z, (unsigned long)t_val);
        mpz_add(td->local_sum, td->local_sum, td->tmp_contrib);
        return;
    }

    for (int pi = max_part_idx; pi >= 0; pi--) {
        int k = odd_parts[pi];
        if (k > remaining) continue;
        int max_m = remaining / k;

        for (int m = 1; m <= max_m; m++) {
            long long dt_self = (long long)m * (m - 1) * k / 2 + (long long)m * (k - 1) / 2;
            long long dt_cross = 0;
            for (int s = 0; s < td->stack_size; s++) {
                int sk = odd_parts[td->stack[s].k_idx];
                dt_cross += (long long)m * td->stack[s].m * gcd_tab[k][sk];
            }

            td->stack[td->stack_size].k_idx = pi;
            td->stack[td->stack_size].m = m;
            td->stack_size++;

            thread_enumerate(td, remaining - m * k, pi - 1, t_val + dt_self + dt_cross);

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

        long long dt_self = (long long)m * (m - 1) * k / 2 + (long long)m * (k - 1) / 2;

        td->stack[0].k_idx = pi;
        td->stack[0].m = m;
        td->stack_size = 1;

        thread_enumerate(td, n - m * k, pi - 1, dt_self);
    }

    return NULL;
}

void a000568(int n, int num_threads) {
    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    if (n <= 1) {
        printf("a(%d) = 1\n", n);
        return;
    }

    global_n = n;
    init_gcd(n);

    num_odd = 0;
    for (int k = 1; k <= n; k += 2)
        odd_parts[num_odd++] = k;

    /* Init factorials */
    for (int i = 0; i <= n; i++)
        mpz_init(fact[i]);
    mpz_set_ui(fact[0], 1);
    for (int i = 1; i <= n; i++)
        mpz_mul_ui(fact[i], fact[i - 1], i);

    /* Precompute km_fac[pi][m] = odd_parts[pi]^m * m! */
    for (int pi = 0; pi < num_odd; pi++) {
        int k = odd_parts[pi];
        int max_m = n / k;
        km_max_m[pi] = max_m;
        km_fac[pi] = malloc((max_m + 1) * sizeof(mpz_t));
        for (int m = 0; m <= max_m; m++) {
            mpz_init(km_fac[pi][m]);
        }
        mpz_set_ui(km_fac[pi][0], 1); /* k^0 * 0! = 1 */
        for (int m = 1; m <= max_m; m++) {
            /* km_fac[pi][m] = km_fac[pi][m-1] * k * m */
            /* Because k^m * m! = k^{m-1} * (m-1)! * k * m */
            mpz_mul_ui(km_fac[pi][m], km_fac[pi][m-1], k);
            mpz_mul_ui(km_fac[pi][m], km_fac[pi][m], m);
        }
    }

    /* Compute LCD */
    mpz_init_set_ui(LCD, 1);
    for (int pi = 0; pi < num_odd; pi++) {
        int max_m = km_max_m[pi];
        mpz_mul(LCD, LCD, km_fac[pi][max_m]);
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

    if (num_threads < 1) num_threads = 1;
    if (num_threads > num_work_items) num_threads = num_work_items;
    atomic_store(&work_counter, 0);

    ThreadData *threads = calloc(num_threads, sizeof(ThreadData));
    pthread_t *pthreads = calloc(num_threads, sizeof(pthread_t));

    for (int t = 0; t < num_threads; t++) {
        threads[t].local_count = 0;
        threads[t].stack_size = 0;
        mpz_init_set_ui(threads[t].local_sum, 0);
        mpz_init(threads[t].tmp_z_val);
        mpz_init(threads[t].tmp_contrib);
        mpz_init(threads[t].tmp_lcd_over_z);
    }

    for (int t = 0; t < num_threads; t++)
        pthread_create(&pthreads[t], NULL, thread_worker, &threads[t]);
    for (int t = 0; t < num_threads; t++)
        pthread_join(pthreads[t], NULL);

    /* Sum results */
    mpz_t total_sum;
    mpz_init_set_ui(total_sum, 0);
    long long total_parts = 0;
    for (int t = 0; t < num_threads; t++) {
        mpz_add(total_sum, total_sum, threads[t].local_sum);
        total_parts += threads[t].local_count;
    }

    mpz_t result;
    mpz_init(result);
    mpz_tdiv_q(result, total_sum, LCD);

    mpz_t remainder;
    mpz_init(remainder);
    mpz_tdiv_r(remainder, total_sum, LCD);
    if (mpz_sgn(remainder) != 0) {
        fprintf(stderr, "ERROR: not exactly divisible!\n");
    }
    mpz_clear(remainder);

    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;

    size_t ndigits = mpz_sizeinbase(result, 10);
    printf("a(%d) = ", n);
    mpz_out_str(stdout, 10, result);
    printf("\n");
    printf("(%zu digits, %.1fs, %lld partitions, %d threads)\n",
           ndigits, elapsed, total_parts, num_threads);

    /* Cleanup */
    mpz_clear(result);
    mpz_clear(total_sum);
    mpz_clear(LCD);
    for (int t = 0; t < num_threads; t++) {
        mpz_clear(threads[t].local_sum);
        mpz_clear(threads[t].tmp_z_val);
        mpz_clear(threads[t].tmp_contrib);
        mpz_clear(threads[t].tmp_lcd_over_z);
    }
    free(threads);
    free(pthreads);
    free(work_items);
    for (int pi = 0; pi < num_odd; pi++) {
        for (int m = 0; m <= km_max_m[pi]; m++)
            mpz_clear(km_fac[pi][m]);
        free(km_fac[pi]);
    }
    for (int i = 0; i <= n; i++)
        mpz_clear(fact[i]);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n> [num_threads]\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    int nt = 1;
    if (argc >= 3) nt = atoi(argv[2]);

    a000568(n, nt);
    fflush(stdout);
    return 0;
}
