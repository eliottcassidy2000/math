/*
 * General k-uniform hypergraph enumeration on n unlabeled nodes.
 *
 * a(n,k) = (1/n!) * sum_{lambda ⊢ n} permcount(lambda) * 2^{c_k(lambda)}
 *
 * c_k(lambda) = (1/L) * sum_{t=0}^{L-1} [x^k] prod_a (1+x^{r_a/gcd(t,r_a)})^{m_a*gcd(t,r_a)}
 *
 * where L = lcm(lambda), sigma has cycle type lambda (compressed form (r_a, m_a)).
 *
 * Uses LCD scaling + bucket accumulation + multi-threading.
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o k_uniform_gmp_enum k_uniform_gmp_enum.c -lgmp -lpthread
 *
 * Usage: ./k_uniform_gmp_enum <k> <n> [threads]
 *
 * Sequences:
 *   k=2: A000088 (simple graphs)
 *   k=3: A000665 (3-uniform hypergraphs)
 *   k=4: A051240 (4-uniform hypergraphs)
 *   k=5+: higher uniform hypergraphs
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

static int global_n, global_k;
static int num_parts;
static int *parts_list;
static int *max_mult;
static mpz_t **km_fac;
static mpz_t LCD;
static int max_t_val;

typedef struct { int pi; int m; } WorkItem;
static WorkItem *work_items;
static int num_work_items;
static atomic_int work_counter;

typedef struct {
    int thread_id;
    mpz_t *buckets;
    long long count;
} ThreadData;

static int gcd_func(int a, int b) {
    while (b) { int t = b; b = a % b; a = t; }
    return a;
}

static long long lcm_func(long long a, long long b) {
    return a / gcd_func((int)a, (int)b) * b;
}

/* Binomial coefficient (for small values) */
static long long binomial(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    long long result = 1;
    for (int i = 0; i < k; i++) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

/*
 * Compute c_k(lambda) for a partition in compressed form.
 * pk[i] = part index, pm[i] = multiplicity, depth = number of distinct parts.
 *
 * c_k = (1/L) * sum_{t=0}^{L-1} f_k(t)
 * where f_k(t) = [x^k] prod_a (1 + x^{r_a/gcd(t,r_a)})^{m_a*gcd(t,r_a)}
 */
static long long compute_ck(int *pk, int *pm, int depth) {
    int k = global_k;

    /* Compute L = lcm of all parts used */
    long long L = 1;
    for (int i = 0; i < depth; i++) {
        L = lcm_func(L, parts_list[pk[i]]);
    }

    /* Sum over t = 0..L-1 */
    long long total = 0;

    for (long long t = 0; t < L; t++) {
        /* Compute [x^k] product of (1 + x^len)^cnt for each part */
        /* Since k is small, we use a fixed-size polynomial */
        long long poly[21];  /* support k up to 20 */
        memset(poly, 0, sizeof(poly));
        poly[0] = 1;

        for (int i = 0; i < depth; i++) {
            int r = parts_list[pk[i]], m = pm[i];
            int g = gcd_func((int)(t % r == 0 ? r : t % r), r);
            if (t == 0) g = r;
            else g = gcd_func((int)(t % r), r);
            /* Correct: gcd(t, r) where t can be larger than r */
            g = gcd_func((int)(t % r), r);
            if (t == 0) g = r;

            int length = r / g;
            int cnt = m * g;

            /* Multiply poly by (1 + x^length)^cnt, truncated to degree k */
            long long factor[21];
            memset(factor, 0, sizeof(factor));
            for (int j = 0; j <= cnt; j++) {
                int deg = j * length;
                if (deg > k) break;
                factor[deg] = binomial(cnt, j);
            }

            long long new_poly[21];
            memset(new_poly, 0, sizeof(new_poly));
            for (int a = 0; a <= k; a++) {
                if (poly[a] == 0) continue;
                for (int b = 0; b <= k - a; b++) {
                    if (factor[b] == 0) continue;
                    new_poly[a + b] += poly[a] * factor[b];
                }
            }
            memcpy(poly, new_poly, sizeof(poly));
        }

        total += poly[k];
    }

    /* c_k = total / L */
    if (total % L != 0) {
        fprintf(stderr, "ERROR: c_k not integer: %lld/%lld\n", total, L);
        return -1;
    }
    return total / L;
}

static void enumerate(
    ThreadData *td, int remaining, int max_pi,
    int *pk, int *pm, int depth, mpz_t lcd_over_z
) {
    if (remaining == 0) {
        td->count++;
        long long ck = compute_ck(pk, pm, depth);
        if (ck < 0 || ck > max_t_val) {
            fprintf(stderr, "ERROR: ck=%lld out of range [0,%d]\n", ck, max_t_val);
            return;
        }
        mpz_add(td->buckets[ck], td->buckets[ck], lcd_over_z);
        return;
    }

    for (int pi = max_pi; pi >= 0; pi--) {
        int k = parts_list[pi];
        if (k > remaining) continue;
        int mm = remaining / k;

        mpz_t cur_lcd;
        mpz_init(cur_lcd);
        for (int m = 1; m <= mm; m++) {
            mpz_divexact(cur_lcd, lcd_over_z, km_fac[pi][m]);
            if (m > 1) mpz_mul(cur_lcd, cur_lcd, km_fac[pi][m-1]);
            mpz_divexact(cur_lcd, lcd_over_z, km_fac[pi][m]);

            pk[depth] = pi;
            pm[depth] = m;
            enumerate(td, remaining - m * k, pi - 1, pk, pm, depth + 1, cur_lcd);
        }
        mpz_clear(cur_lcd);
    }
}

static void *thread_worker(void *arg) {
    ThreadData *td = (ThreadData *)arg;
    int *pk = (int *)calloc(global_n + 1, sizeof(int));
    int *pm = (int *)calloc(global_n + 1, sizeof(int));
    mpz_t lcd_over_z;
    mpz_init(lcd_over_z);

    while (1) {
        int idx = atomic_fetch_add(&work_counter, 1);
        if (idx >= num_work_items) break;

        int pi = work_items[idx].pi;
        int m = work_items[idx].m;

        mpz_divexact(lcd_over_z, LCD, km_fac[pi][m]);
        pk[0] = pi; pm[0] = m;

        int remaining = global_n - m * parts_list[pi];
        if (remaining == 0) {
            td->count++;
            long long ck = compute_ck(pk, pm, 1);
            if (ck >= 0 && ck <= max_t_val) {
                mpz_add(td->buckets[ck], td->buckets[ck], lcd_over_z);
            }
        } else {
            enumerate(td, remaining, pi - 1, pk, pm, 1, lcd_over_z);
        }
    }

    free(pk); free(pm);
    mpz_clear(lcd_over_z);
    return NULL;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <k> <n> [threads]\n", argv[0]);
        return 1;
    }

    global_k = atoi(argv[1]);
    global_n = atoi(argv[2]);
    int n = global_n;
    int k = global_k;
    int num_threads = (argc >= 4) ? atoi(argv[3]) : 1;

    if (k > 20) {
        fprintf(stderr, "k must be <= 20\n");
        return 1;
    }

    if (n < k) {
        printf("a(%d,%d) = 1\n", n, k);
        return 0;
    }

    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    /* Build parts list */
    num_parts = n;
    parts_list = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) parts_list[i] = i + 1;

    max_mult = (int *)malloc(num_parts * sizeof(int));
    for (int pi = 0; pi < num_parts; pi++)
        max_mult[pi] = n / parts_list[pi];

    /* Precompute km_fac[pi][m] = parts[pi]^m * m! */
    km_fac = (mpz_t **)malloc(num_parts * sizeof(mpz_t *));
    for (int pi = 0; pi < num_parts; pi++) {
        int mm = max_mult[pi];
        km_fac[pi] = (mpz_t *)malloc((mm + 1) * sizeof(mpz_t));
        mpz_init_set_ui(km_fac[pi][0], 1);
        for (int m = 1; m <= mm; m++) {
            mpz_init(km_fac[pi][m]);
            mpz_mul_ui(km_fac[pi][m], km_fac[pi][m-1], parts_list[pi]);
            mpz_mul_ui(km_fac[pi][m], km_fac[pi][m], m);
        }
    }

    mpz_init_set_ui(LCD, 1);
    for (int pi = 0; pi < num_parts; pi++)
        mpz_mul(LCD, LCD, km_fac[pi][max_mult[pi]]);

    /* Max c_k value: at most C(n, k) */
    long long cnk = 1;
    for (int i = 0; i < k; i++) cnk = cnk * (n - i) / (i + 1);
    max_t_val = (int)cnk + 1;
    if (max_t_val > 50000000) max_t_val = 50000000;

    /* Work items */
    work_items = (WorkItem *)malloc(n * num_parts * sizeof(WorkItem));
    num_work_items = 0;
    for (int pi = num_parts - 1; pi >= 0; pi--)
        for (int m = 1; m <= max_mult[pi]; m++) {
            work_items[num_work_items].pi = pi;
            work_items[num_work_items].m = m;
            num_work_items++;
        }
    atomic_store(&work_counter, 0);

    /* Threads */
    pthread_t *threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
    ThreadData *tdata = (ThreadData *)malloc(num_threads * sizeof(ThreadData));

    for (int t = 0; t < num_threads; t++) {
        tdata[t].thread_id = t;
        tdata[t].count = 0;
        tdata[t].buckets = (mpz_t *)malloc((max_t_val + 1) * sizeof(mpz_t));
        for (int i = 0; i <= max_t_val; i++) mpz_init(tdata[t].buckets[i]);
        pthread_create(&threads[t], NULL, thread_worker, &tdata[t]);
    }

    for (int t = 0; t < num_threads; t++)
        pthread_join(threads[t], NULL);

    /* Merge buckets */
    mpz_t *merged = (mpz_t *)malloc((max_t_val + 1) * sizeof(mpz_t));
    for (int i = 0; i <= max_t_val; i++) mpz_init(merged[i]);
    long long total_count = 0;

    for (int t = 0; t < num_threads; t++) {
        total_count += tdata[t].count;
        for (int i = 0; i <= max_t_val; i++)
            mpz_add(merged[i], merged[i], tdata[t].buckets[i]);
    }

    /* Assemble total = sum merged[t] * 2^t */
    mpz_t total, temp;
    mpz_init(total);
    mpz_init(temp);
    for (int t = max_t_val; t >= 0; t--) {
        if (mpz_sgn(merged[t]) != 0) {
            mpz_mul_2exp(temp, merged[t], t);
            mpz_add(total, total, temp);
        }
    }

    /* a(n,k) = total / LCD  (since total = sum LCD/z_lambda * 2^c_k) */
    mpz_divexact(total, total, LCD);

    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = (t_end.tv_sec - t_start.tv_sec) +
                     (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

    gmp_printf("a(%d,%d) = %Zd\n", n, k, total);
    printf("Time: %.3fs, Partitions: %lld, Threads: %d\n",
           elapsed, total_count, num_threads);

    /* Cleanup */
    mpz_clear(total); mpz_clear(temp); mpz_clear(LCD);
    for (int i = 0; i <= max_t_val; i++) mpz_clear(merged[i]);
    free(merged);
    for (int t = 0; t < num_threads; t++) {
        for (int i = 0; i <= max_t_val; i++) mpz_clear(tdata[t].buckets[i]);
        free(tdata[t].buckets);
    }
    free(tdata); free(threads); free(work_items);
    for (int pi = 0; pi < num_parts; pi++) {
        for (int m = 0; m <= max_mult[pi]; m++) mpz_clear(km_fac[pi][m]);
        free(km_fac[pi]);
    }
    free(km_fac); free(max_mult); free(parts_list);

    return 0;
}
