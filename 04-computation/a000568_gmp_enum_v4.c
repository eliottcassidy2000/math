/*
 * a000568_gmp_enum_v4.c — A000568: bucket accumulation by t-value.
 *
 * Key optimization: Instead of computing LCD/z_val * 2^t for each partition
 * and adding to a single huge accumulator, we bucket LCD/z_val by t-value.
 * Then at the end: total = sum_t bucket[t] * 2^t.
 *
 * This dramatically reduces the size of intermediate GMP numbers:
 * - bucket[t] values are much smaller than the full sum
 * - The final assembly does O(t_max) shifts instead of one per partition
 * - Each mpz_add operates on smaller numbers (no huge shifted values)
 *
 * Compile:
 *   gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *       -o a000568_gmp_enum_v4 a000568_gmp_enum_v4.c -lgmp -pthread
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
/* Max t-value: t <= C(n,2)/2 ~ n^2/4 for n=500 that's ~62500 */
#define MAX_T 200000

/* Shared read-only data */
static int gcd_tab[MAX_N + 1][MAX_N + 1];
static int odd_parts[MAX_PARTS];
static int num_odd;
static mpz_t fact[MAX_N + 1];
static mpz_t LCD;
static int global_n;

/* Precomputed km_fac[pi][m] = k^m * m! */
static mpz_t *km_fac[MAX_PARTS];
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
typedef struct { int k_idx; int m; } PartEntry;

typedef struct {
    long long local_count;
    PartEntry stack[MAX_PARTS];
    int stack_size;
    mpz_t tmp_z_val, tmp_lcd_over_z;
    /* Buckets indexed by t-value */
    mpz_t *buckets;
    int max_t_seen;
    int num_buckets;
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

        /* Add to bucket[t_val] instead of shifting and adding to total */
        int t = (int)t_val;
        if (t > td->max_t_seen) td->max_t_seen = t;
        mpz_add(td->buckets[t], td->buckets[t], td->tmp_lcd_over_z);
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

    /* Precompute km_fac */
    for (int pi = 0; pi < num_odd; pi++) {
        int k = odd_parts[pi];
        int max_m = n / k;
        km_max_m[pi] = max_m;
        km_fac[pi] = malloc((max_m + 1) * sizeof(mpz_t));
        for (int m = 0; m <= max_m; m++)
            mpz_init(km_fac[pi][m]);
        mpz_set_ui(km_fac[pi][0], 1);
        for (int m = 1; m <= max_m; m++) {
            mpz_mul_ui(km_fac[pi][m], km_fac[pi][m-1], k);
            mpz_mul_ui(km_fac[pi][m], km_fac[pi][m], m);
        }
    }

    /* Compute LCD */
    mpz_init_set_ui(LCD, 1);
    for (int pi = 0; pi < num_odd; pi++)
        mpz_mul(LCD, LCD, km_fac[pi][km_max_m[pi]]);

    /* Max possible t value: C(n,2) = n(n-1)/2 (partition 1^n) */
    int max_t = (long long)n * (n - 1) / 2 + 1;
    if (max_t > MAX_T) max_t = MAX_T;

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
        threads[t].max_t_seen = 0;
        threads[t].num_buckets = max_t + 1;
        mpz_init(threads[t].tmp_z_val);
        mpz_init(threads[t].tmp_lcd_over_z);
        threads[t].buckets = malloc((max_t + 1) * sizeof(mpz_t));
        for (int i = 0; i <= max_t; i++)
            mpz_init_set_ui(threads[t].buckets[i], 0);
    }

    for (int t = 0; t < num_threads; t++)
        pthread_create(&pthreads[t], NULL, thread_worker, &threads[t]);
    for (int t = 0; t < num_threads; t++)
        pthread_join(pthreads[t], NULL);

    /* Merge buckets across threads */
    int overall_max_t = 0;
    long long total_parts = 0;
    for (int t = 0; t < num_threads; t++) {
        if (threads[t].max_t_seen > overall_max_t)
            overall_max_t = threads[t].max_t_seen;
        total_parts += threads[t].local_count;
    }

    /* Merge into thread 0's buckets */
    for (int t = 1; t < num_threads; t++) {
        for (int i = 0; i <= overall_max_t; i++) {
            if (mpz_sgn(threads[t].buckets[i]) != 0)
                mpz_add(threads[0].buckets[i], threads[0].buckets[i], threads[t].buckets[i]);
        }
    }

    /* Final assembly: total_sum = sum_t bucket[t] * 2^t */
    mpz_t total_sum;
    mpz_init_set_ui(total_sum, 0);
    mpz_t shifted;
    mpz_init(shifted);

    for (int t = overall_max_t; t >= 0; t--) {
        if (mpz_sgn(threads[0].buckets[t]) != 0) {
            mpz_mul_2exp(shifted, threads[0].buckets[t], (unsigned long)t);
            mpz_add(total_sum, total_sum, shifted);
        }
    }
    mpz_clear(shifted);

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
    printf("(%zu digits, %.1fs, %lld partitions, %d threads, max_t=%d)\n",
           ndigits, elapsed, total_parts, num_threads, overall_max_t);

    /* Cleanup */
    mpz_clear(result);
    mpz_clear(total_sum);
    mpz_clear(LCD);
    for (int t = 0; t < num_threads; t++) {
        mpz_clear(threads[t].tmp_z_val);
        mpz_clear(threads[t].tmp_lcd_over_z);
        for (int i = 0; i <= max_t; i++)
            mpz_clear(threads[t].buckets[i]);
        free(threads[t].buckets);
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
