/*
 * A000665: Number of 3-uniform hypergraphs on n unlabeled nodes.
 *
 * Formula (Andrew Howroyd, OEIS):
 *   a(n) = (1/n!) * sum_{p partition of n} permcount(p) * 2^edges(p)
 *
 * edges(p) in compressed form [(r_a, m_a)] with distinct parts:
 *   Single: sum_a m_a * ceil((r_a-1)*(r_a-2)/6)
 *   Pair-same: sum_a C(m_a,2) * r_a*(r_a-1)
 *   Pair-diff: sum_{a<b} m_a*m_b * [gcd(r_a,r_b)*(r_a+r_b-2+((r_a-r_b)/gcd)%2)/2]
 *   Triple-same: sum_a C(m_a,3) * r_a^2
 *   Triple-two: sum_{a!=b} C(m_a,2)*m_b * r_a*gcd(r_a,r_b)
 *   Triple-diff: sum_{a<b<c} m_a*m_b*m_c * r_a*r_b*r_c/lcm(r_a,r_b,r_c)
 *
 * Uses LCD scaling + bucket accumulation + multi-threading.
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o a000665_gmp_enum a000665_gmp_enum.c -lgmp -lpthread
 *
 * Usage: ./a000665_gmp_enum <n> [threads]
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

static int global_n;
static int num_parts;
static int *parts_list;   /* all parts 1..n */
static int *max_mult;
static mpz_t **km_fac;    /* km_fac[pi][m] = parts[pi]^m * m! */
static mpz_t LCD;
static int **gcd_tab;     /* gcd_tab[pi][pj] = gcd(parts[pi], parts[pj]) */
static int **lcm_tab;     /* lcm_tab[pi][pj] = lcm(parts[pi], parts[pj]) */
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

/*
 * Compute edges in compressed form for 3-uniform hypergraphs.
 * pk[i] = part index, pm[i] = multiplicity, depth = number of distinct parts.
 */
static long long compute_edges_3uniform(int *pk, int *pm, int depth) {
    long long result = 0;

    /* Single-cycle terms: m_a * ceil((r-1)*(r-2)/6) */
    for (int i = 0; i < depth; i++) {
        int r = parts_list[pk[i]], m = pm[i];
        result += (long long)m * (((long long)(r-1)*(r-2) + 5) / 6);
    }

    /* Pair terms */
    for (int i = 0; i < depth; i++) {
        int r = parts_list[pk[i]], mr = pm[i];
        /* Same-part pairs: C(mr,2) * r*(r-1) */
        result += (long long)mr * (mr - 1) / 2 * r * (r - 1);
        /* Cross-part pairs */
        for (int j = i + 1; j < depth; j++) {
            int s = parts_list[pk[j]], ms = pm[j];
            int g = gcd_tab[pk[i]][pk[j]];
            int pair_val = g * (r + s - 2 + ((r - s) / g) % 2) / 2;
            result += (long long)mr * ms * pair_val;
        }
    }

    /* Triple terms */
    for (int i = 0; i < depth; i++) {
        int r = parts_list[pk[i]], mr = pm[i];
        /* All three from same part: C(mr,3) * r^2 */
        result += (long long)mr * (mr - 1) * (mr - 2) / 6 * r * r;
        /* Two from part i, one from part j */
        long long c2 = (long long)mr * (mr - 1) / 2;
        for (int j = 0; j < depth; j++) {
            if (j == i) continue;
            int ms = pm[j];
            int g = gcd_tab[pk[i]][pk[j]];
            result += c2 * ms * r * g;
        }
        /* One each from three distinct parts (i < j < k) */
        for (int j = i + 1; j < depth; j++) {
            int s = parts_list[pk[j]], ms = pm[j];
            int lcm_ij = lcm_tab[pk[i]][pk[j]];
            for (int k = j + 1; k < depth; k++) {
                int t = parts_list[pk[k]], mt = pm[k];
                /* lcm(r,s,t) = lcm(lcm(r,s), t) */
                int lcm_ijk = lcm_ij / gcd_func(lcm_ij, t) * t;
                result += (long long)mr * ms * mt * ((long long)r * s * t / lcm_ijk);
            }
        }
    }

    return result;
}

static void enumerate(
    ThreadData *td, int remaining, int max_pi,
    int *pk, int *pm, int depth, mpz_t lcd_over_z
) {
    if (remaining == 0) {
        td->count++;
        long long edge_val = compute_edges_3uniform(pk, pm, depth);
        if (edge_val < 0 || edge_val > max_t_val) {
            fprintf(stderr, "ERROR: edge_val=%lld out of range [0,%d]\n", edge_val, max_t_val);
            return;
        }
        mpz_add(td->buckets[edge_val], td->buckets[edge_val], lcd_over_z);
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
            long long edge_val = compute_edges_3uniform(pk, pm, 1);
            if (edge_val >= 0 && edge_val <= max_t_val) {
                mpz_add(td->buckets[edge_val], td->buckets[edge_val], lcd_over_z);
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
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n> [threads]\n", argv[0]);
        return 1;
    }

    global_n = atoi(argv[1]);
    int n = global_n;
    int num_threads = (argc >= 3) ? atoi(argv[2]) : 1;

    if (n <= 2) {
        printf("a(%d) = 1\n", n);
        return 0;
    }

    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    /* Build parts list (all parts 1..n) */
    num_parts = n;
    parts_list = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) parts_list[i] = i + 1;

    max_mult = (int *)malloc(num_parts * sizeof(int));
    for (int pi = 0; pi < num_parts; pi++)
        max_mult[pi] = n / parts_list[pi];

    /* GCD and LCM tables */
    gcd_tab = (int **)malloc(num_parts * sizeof(int *));
    lcm_tab = (int **)malloc(num_parts * sizeof(int *));
    for (int i = 0; i < num_parts; i++) {
        gcd_tab[i] = (int *)malloc(num_parts * sizeof(int));
        lcm_tab[i] = (int *)malloc(num_parts * sizeof(int));
        for (int j = 0; j < num_parts; j++) {
            gcd_tab[i][j] = gcd_func(parts_list[i], parts_list[j]);
            lcm_tab[i][j] = parts_list[i] / gcd_tab[i][j] * parts_list[j];
        }
    }

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

    /* LCD = product of all km_fac[pi][max_mult[pi]] */
    mpz_init_set_ui(LCD, 1);
    for (int pi = 0; pi < num_parts; pi++)
        mpz_mul(LCD, LCD, km_fac[pi][max_mult[pi]]);

    /* Max t value: edges for 3-uniform is bounded by C(n,3) */
    /* Conservative: n^3/6 + n^2 */
    max_t_val = (long long)n * n * n / 6 + (long long)n * n + n + 1;
    if (max_t_val > 100000000) max_t_val = 100000000; /* safety cap */

    /* Build work items */
    work_items = (WorkItem *)malloc(n * num_parts * sizeof(WorkItem));
    num_work_items = 0;
    for (int pi = num_parts - 1; pi >= 0; pi--)
        for (int m = 1; m <= max_mult[pi]; m++) {
            work_items[num_work_items].pi = pi;
            work_items[num_work_items].m = m;
            num_work_items++;
        }
    atomic_store(&work_counter, 0);

    /* Create threads */
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

    /* Assemble: total = sum_t merged[t] * 2^t */
    mpz_t total, temp;
    mpz_init(total);
    mpz_init(temp);
    for (int t = max_t_val; t >= 0; t--) {
        if (mpz_sgn(merged[t]) != 0) {
            mpz_mul_2exp(temp, merged[t], t);
            mpz_add(total, total, temp);
        }
    }

    /* total = sum of (LCD/z_lambda) * 2^edges
     * a(n) = (1/n!) * sum_p permcount(p) * 2^edges(p)
     *       = (1/n!) * sum_p (n!/z_p) * 2^edges(p)
     *       = sum_p 2^edges(p) / z_p
     *       = (1/LCD) * sum_p (LCD/z_p) * 2^edges(p)
     *       = total / LCD
     * Wait... permcount(p) = n!/z_p, so
     *   a(n) = (1/n!) * sum_p (n!/z_p) * 2^edges = sum_p 2^edges/z_p
     *   = (1/LCD) * total
     * Hmm, but we have total = sum of (LCD/z_p) * 2^edges.
     * So a(n) = total / LCD.
     * No wait, the full formula is a(n) = (1/n!) * sum_p permcount(p) * 2^edges.
     * permcount(p) = s!/prod(p_i * k_i) where p is expanded.
     * In compressed form, z_lambda = prod(r^m * m!) for (r,m).
     * permcount = n!/z_lambda.
     * So a(n) = (1/n!) * sum (n!/z_lambda) * 2^edges = sum 2^edges/z_lambda.
     * With LCD scaling: total = sum (LCD/z_lambda) * 2^edges.
     * a(n) = total / LCD.
     */

    mpz_divexact(total, total, LCD);

    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = (t_end.tv_sec - t_start.tv_sec) +
                     (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

    gmp_printf("a(%d) = %Zd\n", n, total);
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
        free(km_fac[pi]); free(gcd_tab[pi]); free(lcm_tab[pi]);
    }
    free(km_fac); free(gcd_tab); free(lcm_tab); free(max_mult); free(parts_list);

    return 0;
}
