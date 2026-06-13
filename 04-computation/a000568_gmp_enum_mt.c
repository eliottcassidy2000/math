/*
 * a000568_gmp_enum_mt.c — A000568 via multi-threaded partition enumeration with GMP.
 *
 * Uses a work-queue of top-level (k,m) choices distributed via atomic counter.
 * Each thread picks the next work item, processes the entire subtree,
 * and accumulates into its thread-local sum.
 *
 * Compile:
 *   gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *       -o a000568_gmp_enum_mt a000568_gmp_enum_mt.c -lgmp -pthread
 *
 * Usage:
 *   ./a000568_gmp_enum_mt <n> [num_threads]
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

/* Work queue: list of (pi, m) top-level choices */
typedef struct { int pi; int m; } WorkItem;
static WorkItem *work_items;
static int num_work_items;
static atomic_int work_counter;

/* Per-thread data */
typedef struct { int k; int m; } PartEntry;

typedef struct {
    int thread_id;
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

        /* Compute z_val */
        mpz_set_ui(td->tmp_z_val, 1);
        for (int i = 0; i < td->stack_size; i++) {
            int k = td->stack[i].k;
            int m = td->stack[i].m;
            mpz_t km;
            mpz_init(km);
            mpz_ui_pow_ui(km, k, m);
            mpz_mul(td->tmp_z_val, td->tmp_z_val, km);
            mpz_mul(td->tmp_z_val, td->tmp_z_val, fact[m]);
            mpz_clear(km);
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
                dt_cross += (long long)m * td->stack[s].m * gcd_tab[k][td->stack[s].k];
            }

            td->stack[td->stack_size].k = k;
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

        td->stack[0].k = k;
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

    /* Compute LCD */
    mpz_init_set_ui(LCD, 1);
    for (int k = 1; k <= n; k += 2) {
        int mk = n / k;
        mpz_t kpow;
        mpz_init(kpow);
        mpz_ui_pow_ui(kpow, k, mk);
        mpz_mul(LCD, LCD, kpow);
        mpz_mul(LCD, LCD, fact[mk]);
        mpz_clear(kpow);
    }

    /* Build work queue: all (pi, m) top-level choices */
    num_work_items = 0;
    for (int pi = num_odd - 1; pi >= 0; pi--) {
        int k = odd_parts[pi];
        if (k > n) continue;
        num_work_items += n / k;
    }
    work_items = malloc(num_work_items * sizeof(WorkItem));
    int wi = 0;
    /* Order: largest parts first (they tend to have smaller subtrees,
     * so mixing them gives better load balance) */
    for (int pi = num_odd - 1; pi >= 0; pi--) {
        int k = odd_parts[pi];
        if (k > n) continue;
        int max_m = n / k;
        for (int m = 1; m <= max_m; m++) {
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
        threads[t].thread_id = t;
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

    /* Divide by LCD */
    mpz_t result;
    mpz_init(result);
    mpz_tdiv_q(result, total_sum, LCD);

    /* Verify */
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
    if (argc >= 3)
        nt = atoi(argv[2]);

    a000568(n, nt);
    fflush(stdout);

    return 0;
}
