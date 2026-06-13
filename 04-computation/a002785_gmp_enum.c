/*
 * A002785: Self-complementary oriented graphs (self-converse tournaments).
 *
 * Formula (Andrew Howroyd, OEIS):
 *   a(n) = (1/n!) * sum_{p partition of floor(n/2), all odd parts}
 *          permcount(2*p) * 2^edges(p) * (n*2^#p if n odd, 1 if even)
 *
 * Compressed partition form with (part_k, mult_m):
 *   edges = sum_r m_r^2 * r + 2 * sum_{r<s} m_r * m_s * gcd(r,s)
 *   z_{2p} = prod (2r)^m_r * m_r!
 *   permcount(2p) = (2k)! / z_{2p}   where k = floor(n/2)
 *
 * Uses LCD scaling + bucket accumulation (same as a000568_gmp_enum_v4).
 * Multi-threaded with atomic work queue.
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o a002785_gmp_enum a002785_gmp_enum.c -lgmp -lpthread
 *
 * Usage: ./a002785_gmp_enum <n> [threads]
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

/* Globals */
static int global_n;
static int global_k;     /* floor(n/2) */
static int is_odd_n;
static int num_odd_parts;
static int *odd_parts;   /* array of odd parts 1,3,5,...,<=k */

/* Precomputed: km_fac[pi][m] = (2*k)^m * m!  for z_{2λ} */
static mpz_t **km_fac;   /* km_fac[pi][m] for pi=0..num_odd_parts-1, m=0..max_m */
static int *max_mult;    /* max multiplicity of each part */

/* LCD = product of all z_{2λ} denominators */
static mpz_t LCD;

/* GCD table */
static int **gcd_tab;

/* Bucket accumulation: buckets[t] += LCD/z_{2λ} * extra */
/* t ranges from 0 to max_t */
static int max_t_val;

/* Work queue for multi-threading */
typedef struct {
    int pi;  /* part index */
    int m;   /* multiplicity */
} WorkItem;

static WorkItem *work_items;
static int num_work_items;
static atomic_int work_counter;

/* Thread-local data */
typedef struct {
    int thread_id;
    mpz_t *buckets;
    long long count;
} ThreadData;

static int gcd_func(int a, int b) {
    while (b) { int t = b; b = a % b; a = t; }
    return a;
}

/*
 * Recursive partition enumeration.
 * parts_k[i], parts_m[i] for i=0..depth-1 are the chosen (part, mult) pairs.
 * lcd_over_z is LCD / z_{2λ_so_far} (maintained incrementally).
 * t_so_far is the edges value accumulated so far.
 */
static void enumerate(
    ThreadData *td,
    int remaining,
    int max_pi,
    int *parts_k, int *parts_m, int depth,
    long long t_so_far,
    int total_parts_so_far,
    mpz_t lcd_over_z  /* LCD / z_{2λ_partial} */
) {
    if (remaining == 0) {
        td->count++;

        /* t_val = edges + (num_parts if n is odd) */
        long long t_val = t_so_far;
        if (is_odd_n) {
            t_val += total_parts_so_far;
        }

        if (t_val < 0 || t_val > max_t_val) {
            fprintf(stderr, "ERROR: t_val=%lld out of range [0,%d]\n", t_val, max_t_val);
            return;
        }

        /* Extra factor for odd n: n * 2^#parts
         * But 2^#parts is absorbed into 2^t_val above.
         * We still need the factor of n. */
        if (is_odd_n) {
            mpz_t tmp;
            mpz_init(tmp);
            mpz_mul_ui(tmp, lcd_over_z, global_n);
            mpz_add(td->buckets[t_val], td->buckets[t_val], tmp);
            mpz_clear(tmp);
        } else {
            mpz_add(td->buckets[t_val], td->buckets[t_val], lcd_over_z);
        }
        return;
    }

    for (int pi = max_pi; pi >= 0; pi--) {
        int k = odd_parts[pi];
        if (k > remaining) continue;
        int mm = remaining / k;
        if (mm > max_mult[pi]) mm = max_mult[pi];

        mpz_t cur_lcd;
        mpz_init_set(cur_lcd, lcd_over_z);

        for (int m = 1; m <= mm; m++) {
            /* Divide lcd_over_z by (2k)^m * m! incrementally:
             * At step m: divide by (2k) * m */
            mpz_divexact(cur_lcd, cur_lcd, km_fac[pi][m]);
            if (m > 1) {
                mpz_mul(cur_lcd, cur_lcd, km_fac[pi][m-1]);
            }
            /* Actually, km_fac[pi][m] = (2k)^m * m!
             * We want LCD / (z_{prev} * (2k)^m * m!)
             * = (LCD / z_{prev}) / ((2k)^m * m!)
             * = lcd_over_z / km_fac[pi][m]
             *
             * But we're doing this incrementally. Let me just divide fresh.
             */
            mpz_divexact(cur_lcd, lcd_over_z, km_fac[pi][m]);

            /* Compute t contribution from this part */
            /* Self: m^2 * k */
            long long t_self = (long long)m * m * k;

            /* Cross with previous parts */
            long long t_cross = 0;
            for (int d = 0; d < depth; d++) {
                t_cross += 2LL * m * parts_m[d] * gcd_tab[pi][parts_k[d]];
                /* Note: parts_k[d] is the part INDEX, need actual gcd */
            }

            parts_k[depth] = pi;  /* store part index for gcd lookup */
            parts_m[depth] = m;

            enumerate(td, remaining - m * k, pi - 1,
                     parts_k, parts_m, depth + 1,
                     t_so_far + t_self + t_cross,
                     total_parts_so_far + m,
                     cur_lcd);
        }
        mpz_clear(cur_lcd);
    }
}

static void *thread_worker(void *arg) {
    ThreadData *td = (ThreadData *)arg;
    int k = global_k;

    /* Allocate local storage */
    int *parts_k = (int *)calloc(k + 1, sizeof(int));
    int *parts_m = (int *)calloc(k + 1, sizeof(int));
    mpz_t lcd_over_z;
    mpz_init(lcd_over_z);

    while (1) {
        int idx = atomic_fetch_add(&work_counter, 1);
        if (idx >= num_work_items) break;

        int pi = work_items[idx].pi;
        int m = work_items[idx].m;
        int part = odd_parts[pi];

        /* Start with LCD / km_fac[pi][m] */
        mpz_divexact(lcd_over_z, LCD, km_fac[pi][m]);

        /* t from this first part: m^2 * part */
        long long t_self = (long long)m * m * part;

        parts_k[0] = pi;
        parts_m[0] = m;

        int remaining = k - m * part;
        if (remaining == 0) {
            td->count++;
            long long t_val = t_self;
            int total_parts = m;
            if (is_odd_n) t_val += total_parts;

            if (t_val >= 0 && t_val <= max_t_val) {
                if (is_odd_n) {
                    mpz_t tmp;
                    mpz_init(tmp);
                    mpz_mul_ui(tmp, lcd_over_z, global_n);
                    mpz_add(td->buckets[t_val], td->buckets[t_val], tmp);
                    mpz_clear(tmp);
                } else {
                    mpz_add(td->buckets[t_val], td->buckets[t_val], lcd_over_z);
                }
            }
        } else {
            enumerate(td, remaining, pi - 1,
                     parts_k, parts_m, 1,
                     t_self, m, lcd_over_z);
        }
    }

    free(parts_k);
    free(parts_m);
    mpz_clear(lcd_over_z);
    return NULL;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n> [threads]\n", argv[0]);
        return 1;
    }

    global_n = atoi(argv[1]);
    int n = global_n;
    int num_threads = (argc >= 3) ? atoi(argv[2]) : 1;

    if (n <= 1) {
        printf("a(%d) = 1\n", n);
        return 0;
    }

    global_k = n / 2;
    int k = global_k;
    is_odd_n = n % 2;

    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    /* Build odd parts list */
    num_odd_parts = 0;
    odd_parts = (int *)malloc((k + 1) * sizeof(int));
    for (int i = 1; i <= k; i += 2) {
        odd_parts[num_odd_parts++] = i;
    }

    /* Precompute max multiplicities */
    max_mult = (int *)malloc(num_odd_parts * sizeof(int));
    for (int pi = 0; pi < num_odd_parts; pi++) {
        max_mult[pi] = k / odd_parts[pi];
    }

    /* Precompute GCD table: gcd_tab[pi1][pi2] = gcd(odd_parts[pi1], odd_parts[pi2]) */
    gcd_tab = (int **)malloc(num_odd_parts * sizeof(int *));
    for (int i = 0; i < num_odd_parts; i++) {
        gcd_tab[i] = (int *)malloc(num_odd_parts * sizeof(int));
        for (int j = 0; j < num_odd_parts; j++) {
            gcd_tab[i][j] = gcd_func(odd_parts[i], odd_parts[j]);
        }
    }

    /* Precompute km_fac[pi][m] = (2*odd_parts[pi])^m * m! */
    km_fac = (mpz_t **)malloc(num_odd_parts * sizeof(mpz_t *));
    for (int pi = 0; pi < num_odd_parts; pi++) {
        int mm = max_mult[pi];
        km_fac[pi] = (mpz_t *)malloc((mm + 1) * sizeof(mpz_t));
        mpz_init_set_ui(km_fac[pi][0], 1);
        for (int m = 1; m <= mm; m++) {
            mpz_init(km_fac[pi][m]);
            mpz_mul_ui(km_fac[pi][m], km_fac[pi][m-1], 2 * odd_parts[pi]);
            mpz_mul_ui(km_fac[pi][m], km_fac[pi][m], m);
        }
    }

    /* Compute LCD = product of all km_fac[pi][max_mult[pi]] */
    mpz_init_set_ui(LCD, 1);
    for (int pi = 0; pi < num_odd_parts; pi++) {
        mpz_mul(LCD, LCD, km_fac[pi][max_mult[pi]]);
    }

    /* Max t value: edges can be at most sum over all pairs = k*(k-1)/2 * max_gcd + k
     * Conservative bound: k^2 (since parts sum to k) plus k for odd_n extra */
    max_t_val = (long long)k * k + k + 1;

    /* Build work items: all (pi, m) top-level choices */
    work_items = (WorkItem *)malloc(k * num_odd_parts * sizeof(WorkItem));
    num_work_items = 0;
    for (int pi = num_odd_parts - 1; pi >= 0; pi--) {
        for (int m = 1; m <= max_mult[pi]; m++) {
            work_items[num_work_items].pi = pi;
            work_items[num_work_items].m = m;
            num_work_items++;
        }
    }
    atomic_store(&work_counter, 0);

    /* Create threads */
    pthread_t *threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
    ThreadData *tdata = (ThreadData *)malloc(num_threads * sizeof(ThreadData));

    for (int t = 0; t < num_threads; t++) {
        tdata[t].thread_id = t;
        tdata[t].count = 0;
        tdata[t].buckets = (mpz_t *)malloc((max_t_val + 1) * sizeof(mpz_t));
        for (int i = 0; i <= max_t_val; i++) {
            mpz_init(tdata[t].buckets[i]);
        }
        pthread_create(&threads[t], NULL, thread_worker, &tdata[t]);
    }

    for (int t = 0; t < num_threads; t++) {
        pthread_join(threads[t], NULL);
    }

    /* Merge buckets */
    mpz_t *merged = (mpz_t *)malloc((max_t_val + 1) * sizeof(mpz_t));
    for (int i = 0; i <= max_t_val; i++) {
        mpz_init(merged[i]);
    }
    long long total_count = 0;
    for (int t = 0; t < num_threads; t++) {
        total_count += tdata[t].count;
        for (int i = 0; i <= max_t_val; i++) {
            mpz_add(merged[i], merged[i], tdata[t].buckets[i]);
        }
    }

    /* Assemble: total = sum_t merged[t] * 2^t */
    /* Then result = total * (2k)! / (LCD * n!) */
    mpz_t total, temp;
    mpz_init(total);
    mpz_init(temp);
    for (int t = max_t_val; t >= 0; t--) {
        if (mpz_sgn(merged[t]) != 0) {
            mpz_mul_2exp(temp, merged[t], t);
            mpz_add(total, total, temp);
        }
    }

    /* total = sum of LCD/z_{2λ} * extra * 2^edges
     * We need: sum of permcount(2λ) * 2^edges * extra / n!
     * = sum of (2k)!/z_{2λ} * 2^edges * extra / n!
     * = (2k)! / (LCD * n!) * sum of LCD/z_{2λ} * 2^edges * extra
     * = (2k)! / (LCD * n!) * total
     */
    mpz_t fact_2k, fact_n;
    mpz_init(fact_2k);
    mpz_init(fact_n);
    mpz_fac_ui(fact_2k, 2 * k);
    mpz_fac_ui(fact_n, n);

    mpz_mul(total, total, fact_2k);
    mpz_t denom;
    mpz_init(denom);
    mpz_mul(denom, LCD, fact_n);
    mpz_divexact(total, total, denom);

    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = (t_end.tv_sec - t_start.tv_sec) +
                     (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

    gmp_printf("a(%d) = %Zd\n", n, total);
    printf("Time: %.3fs, Partitions: %lld, Threads: %d\n",
           elapsed, total_count, num_threads);

    /* Cleanup */
    mpz_clear(total);
    mpz_clear(temp);
    mpz_clear(fact_2k);
    mpz_clear(fact_n);
    mpz_clear(denom);
    mpz_clear(LCD);
    for (int i = 0; i <= max_t_val; i++) {
        mpz_clear(merged[i]);
    }
    free(merged);
    for (int t = 0; t < num_threads; t++) {
        for (int i = 0; i <= max_t_val; i++) {
            mpz_clear(tdata[t].buckets[i]);
        }
        free(tdata[t].buckets);
    }
    free(tdata);
    free(threads);
    free(work_items);
    for (int pi = 0; pi < num_odd_parts; pi++) {
        for (int m = 0; m <= max_mult[pi]; m++) {
            mpz_clear(km_fac[pi][m]);
        }
        free(km_fac[pi]);
        free(gcd_tab[pi]);
    }
    free(km_fac);
    free(gcd_tab);
    free(max_mult);
    free(odd_parts);

    return 0;
}
