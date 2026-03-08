/*
 * Unified Burnside/Polya enumeration for graph-counting sequences.
 *
 * Supports:
 *   A000568 - Tournaments         (odd parts of n)
 *   A000273 - Digraphs            (all parts of n)
 *   A000595 - Binary relations    (all parts of n)
 *   A002785 - Self-comp tournaments (odd parts of n/2)
 *   A000171 - Self-comp graphs    (parts of n/4)
 *
 * All use the same Burnside/Polya pattern:
 *   a(n) = (1/n!) * sum_{p in Partitions} permcount(scale*p) * 2^edges(p) * extra
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o burnside_enum burnside_enum.c -lgmp -lpthread
 *
 * Usage: ./burnside_enum <sequence> <n> [threads]
 *   sequence: 568, 273, 595, 2785, 171
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

/* Sequence type */
typedef enum {
    SEQ_A000568,  /* Tournaments: odd parts of n */
    SEQ_A000273,  /* Digraphs: all parts of n */
    SEQ_A000595,  /* Binary relations: all parts of n */
    SEQ_A002785,  /* Self-comp tournaments: odd parts of n/2 */
    SEQ_A000171,  /* Self-comp graphs: parts of n/4 */
} SeqType;

/* Globals */
static SeqType seq_type;
static int global_n;
static int partition_target;  /* what we partition */
static int odd_parts_only;
static int num_parts;
static int *parts;            /* available parts */

static mpz_t **km_fac;       /* km_fac[pi][m] precomputed */
static int *max_mult;
static mpz_t LCD;
static int **gcd_tab;
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
 * Compute edges value for a partition with given parts and multiplicities.
 *
 * A000568: edges(v) = sum_{i<j} gcd(v_i,v_j) + (sum v_i - sum v_i)/2
 *   Actually: 2*sum_{i<j} gcd + sum v_i, then t = (edges - sum v_i)/2 ... no
 *   PARI: edges(v) = {sum(i=2, #v, sum(j=1, i-1, gcd(v[i], v[j]))) + sum(i=1, #v, v[i]\2)}
 *   So edges = sum_{i<j} gcd(v_i,v_j) + sum floor(v_i/2)
 *
 * A000273: edges(v) = 2*sum_{i<j} gcd(v_i,v_j) + sum(v_i - 1)
 *          = 2*sum_{i<j} gcd(v_i,v_j) + sum(v_i) - #v
 *
 * A000595: edges(v) = 2*sum_{i<j} gcd(v_i,v_j) + sum(v_i)
 *
 * A002785: edges(v) = 2*sum_{i<j} gcd(v_i,v_j) + sum(v_i)
 *   (but partition is of n/2)
 *
 * A000171: edges(v) = 4*sum_{i<j} gcd(v_i,v_j) + 2*sum(v_i)
 *
 * In compressed form with (part_k, mult_m):
 *   sum_{i<j} gcd(v_i,v_j) = sum_{r<s} m_r*m_s*gcd(r,s) + sum_r C(m_r,2)*r
 *     where the same-part pairs contribute gcd(r,r)=r
 *   sum v_i = sum_r m_r * r
 *   #v = sum m_r
 *   sum floor(v_i/2) = sum_r m_r * floor(r/2)
 */
static long long compute_edges(int *pk, int *pm, int depth) {
    long long gcd_sum = 0;  /* sum_{i<j} gcd(v_i,v_j) in full vector */
    long long sum_v = 0;
    long long num_v = 0;
    long long sum_half_v = 0;

    for (int i = 0; i < depth; i++) {
        int r = parts[pk[i]];
        int mr = pm[i];
        sum_v += (long long)mr * r;
        num_v += mr;
        sum_half_v += (long long)mr * (r / 2);

        /* Same-part pairs: C(mr,2) * gcd(r,r) = C(mr,2) * r */
        gcd_sum += (long long)mr * (mr - 1) / 2 * r;

        /* Cross-part pairs */
        for (int j = i + 1; j < depth; j++) {
            int s = parts[pk[j]];
            int ms = pm[j];
            gcd_sum += (long long)mr * ms * gcd_tab[pk[i]][pk[j]];
        }
    }

    long long edge_val;
    switch (seq_type) {
        case SEQ_A000568:
            /* edges = sum_{i<j} gcd + sum floor(v_i/2) */
            edge_val = gcd_sum + sum_half_v;
            break;
        case SEQ_A000273:
            /* edges = 2*sum_{i<j} gcd + sum(v_i - 1) = 2*gcd_sum + sum_v - num_v */
            edge_val = 2 * gcd_sum + sum_v - num_v;
            break;
        case SEQ_A000595:
            /* edges = 2*sum_{i<j} gcd + sum(v_i) */
            edge_val = 2 * gcd_sum + sum_v;
            break;
        case SEQ_A002785:
            /* edges = 2*sum_{i<j} gcd + sum(v_i) */
            edge_val = 2 * gcd_sum + sum_v;
            break;
        case SEQ_A000171:
            /* edges = 4*sum_{i<j} gcd + 2*sum(v_i) */
            edge_val = 4 * gcd_sum + 2 * sum_v;
            break;
        default:
            edge_val = 0;
    }
    return edge_val;
}

/* Compute the extra factor (for self-complementary sequences with odd n) */
static void compute_extra(mpz_t result, int *pk, int *pm, int depth) {
    mpz_set_ui(result, 1);

    if (seq_type == SEQ_A002785 && (global_n % 2 == 1)) {
        /* extra = n * 2^#parts */
        int total_parts = 0;
        for (int i = 0; i < depth; i++) total_parts += pm[i];
        mpz_set_ui(result, global_n);
        mpz_mul_2exp(result, result, total_parts);
    } else if (seq_type == SEQ_A000171 && (global_n % 2 == 1)) {
        /* extra = n * 2^#parts */
        int total_parts = 0;
        for (int i = 0; i < depth; i++) total_parts += pm[i];
        mpz_set_ui(result, global_n);
        mpz_mul_2exp(result, result, total_parts);
    }
}

/* Compute z_lambda (or z_{scale*lambda}) */
static void compute_z_scaled(mpz_t result, int *pk, int *pm, int depth) {
    int scale = 1;
    if (seq_type == SEQ_A002785) scale = 2;
    if (seq_type == SEQ_A000171) scale = 4;

    mpz_set_ui(result, 1);
    for (int i = 0; i < depth; i++) {
        int r = parts[pk[i]];
        int mr = pm[i];
        int sr = scale * r;
        /* z contribution: (scale*r)^mr * mr! */
        mpz_t tmp;
        mpz_init(tmp);
        mpz_ui_pow_ui(tmp, sr, mr);
        mpz_mul(result, result, tmp);
        mpz_fac_ui(tmp, mr);
        mpz_mul(result, result, tmp);
        mpz_clear(tmp);
    }
}

static void enumerate(
    ThreadData *td,
    int remaining,
    int max_pi,
    int *pk, int *pm, int depth,
    mpz_t lcd_over_z
) {
    if (remaining == 0) {
        td->count++;

        long long edge_val = compute_edges(pk, pm, depth);

        /* Extra factor */
        mpz_t extra;
        mpz_init(extra);
        compute_extra(extra, pk, pm, depth);

        if (edge_val < 0 || edge_val > max_t_val) {
            fprintf(stderr, "ERROR: edge_val=%lld out of range\n", edge_val);
            mpz_clear(extra);
            return;
        }

        mpz_t contrib;
        mpz_init(contrib);
        mpz_mul(contrib, lcd_over_z, extra);
        mpz_add(td->buckets[edge_val], td->buckets[edge_val], contrib);
        mpz_clear(contrib);
        mpz_clear(extra);
        return;
    }

    for (int pi = max_pi; pi >= 0; pi--) {
        int k = parts[pi];
        if (k > remaining) continue;
        int mm = remaining / k;

        mpz_t cur_lcd;
        mpz_init(cur_lcd);

        for (int m = 1; m <= mm; m++) {
            mpz_divexact(cur_lcd, lcd_over_z, km_fac[pi][m]);
            if (m > 1) mpz_mul(cur_lcd, cur_lcd, km_fac[pi][m-1]);

            /* Actually, cleaner: just divide fresh each time */
            mpz_divexact(cur_lcd, lcd_over_z, km_fac[pi][m]);

            pk[depth] = pi;
            pm[depth] = m;

            enumerate(td, remaining - m * k, pi - 1,
                     pk, pm, depth + 1, cur_lcd);
        }
        mpz_clear(cur_lcd);
    }
}

static void *thread_worker(void *arg) {
    ThreadData *td = (ThreadData *)arg;
    int *pk = (int *)calloc(partition_target + 1, sizeof(int));
    int *pm = (int *)calloc(partition_target + 1, sizeof(int));
    mpz_t lcd_over_z;
    mpz_init(lcd_over_z);

    while (1) {
        int idx = atomic_fetch_add(&work_counter, 1);
        if (idx >= num_work_items) break;

        int pi = work_items[idx].pi;
        int m = work_items[idx].m;
        int part = parts[pi];

        mpz_divexact(lcd_over_z, LCD, km_fac[pi][m]);

        pk[0] = pi;
        pm[0] = m;

        int remaining = partition_target - m * part;
        if (remaining == 0) {
            td->count++;
            long long edge_val = compute_edges(pk, pm, 1);
            if (edge_val >= 0 && edge_val <= max_t_val) {
                mpz_t extra;
                mpz_init(extra);
                compute_extra(extra, pk, pm, 1);
                mpz_t contrib;
                mpz_init(contrib);
                mpz_mul(contrib, lcd_over_z, extra);
                mpz_add(td->buckets[edge_val], td->buckets[edge_val], contrib);
                mpz_clear(contrib);
                mpz_clear(extra);
            }
        } else {
            enumerate(td, remaining, pi - 1, pk, pm, 1, lcd_over_z);
        }
    }

    free(pk);
    free(pm);
    mpz_clear(lcd_over_z);
    return NULL;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <seq: 568|273|595|2785|171> <n> [threads]\n", argv[0]);
        return 1;
    }

    int seq_id = atoi(argv[1]);
    global_n = atoi(argv[2]);
    int n = global_n;
    int num_threads = (argc >= 4) ? atoi(argv[3]) : 1;

    switch (seq_id) {
        case 568:  seq_type = SEQ_A000568; break;
        case 273:  seq_type = SEQ_A000273; break;
        case 595:  seq_type = SEQ_A000595; break;
        case 2785: seq_type = SEQ_A002785; break;
        case 171:  seq_type = SEQ_A000171; break;
        default:
            fprintf(stderr, "Unknown sequence: %d\n", seq_id);
            return 1;
    }

    /* Determine partition target and part types */
    switch (seq_type) {
        case SEQ_A000568:
            partition_target = n;
            odd_parts_only = 1;
            break;
        case SEQ_A000273:
        case SEQ_A000595:
            partition_target = n;
            odd_parts_only = 0;
            break;
        case SEQ_A002785:
            partition_target = n / 2;
            odd_parts_only = 1;
            break;
        case SEQ_A000171:
            partition_target = n / 4;
            odd_parts_only = 0;
            /* a(n) = 0 if n%4 >= 2 */
            if (n % 4 >= 2) {
                printf("a(%d) = 0\n", n);
                return 0;
            }
            break;
    }

    if (partition_target == 0) {
        printf("a(%d) = 1\n", n);
        return 0;
    }

    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    /* Build parts list */
    num_parts = 0;
    parts = (int *)malloc((partition_target + 1) * sizeof(int));
    if (odd_parts_only) {
        for (int i = 1; i <= partition_target; i += 2)
            parts[num_parts++] = i;
    } else {
        for (int i = 1; i <= partition_target; i++)
            parts[num_parts++] = i;
    }

    /* Max multiplicities */
    max_mult = (int *)malloc(num_parts * sizeof(int));
    for (int pi = 0; pi < num_parts; pi++)
        max_mult[pi] = partition_target / parts[pi];

    /* GCD table */
    gcd_tab = (int **)malloc(num_parts * sizeof(int *));
    for (int i = 0; i < num_parts; i++) {
        gcd_tab[i] = (int *)malloc(num_parts * sizeof(int));
        for (int j = 0; j < num_parts; j++)
            gcd_tab[i][j] = gcd_func(parts[i], parts[j]);
    }

    /* Precompute km_fac[pi][m] = (scale*parts[pi])^m * m! */
    int scale = 1;
    if (seq_type == SEQ_A002785) scale = 2;
    if (seq_type == SEQ_A000171) scale = 4;

    km_fac = (mpz_t **)malloc(num_parts * sizeof(mpz_t *));
    for (int pi = 0; pi < num_parts; pi++) {
        int mm = max_mult[pi];
        km_fac[pi] = (mpz_t *)malloc((mm + 1) * sizeof(mpz_t));
        mpz_init_set_ui(km_fac[pi][0], 1);
        for (int m = 1; m <= mm; m++) {
            mpz_init(km_fac[pi][m]);
            mpz_mul_ui(km_fac[pi][m], km_fac[pi][m-1], scale * parts[pi]);
            mpz_mul_ui(km_fac[pi][m], km_fac[pi][m], m);
        }
    }

    /* LCD */
    mpz_init_set_ui(LCD, 1);
    for (int pi = 0; pi < num_parts; pi++)
        mpz_mul(LCD, LCD, km_fac[pi][max_mult[pi]]);

    /* Max t value - conservative upper bound */
    /* For digraphs/relations: edges can be up to ~2*n^2 */
    long long pt = partition_target;
    switch (seq_type) {
        case SEQ_A000568:
            max_t_val = pt * (pt - 1) / 2 + 1;
            break;
        case SEQ_A000273:
            max_t_val = pt * (pt - 1) + pt + 1;
            break;
        case SEQ_A000595:
            max_t_val = pt * (pt - 1) + pt + 1;
            break;
        case SEQ_A002785:
            max_t_val = pt * pt + pt + 1;
            break;
        case SEQ_A000171:
            max_t_val = 2 * pt * (pt - 1) + 2 * pt + 1;
            break;
        default:
            max_t_val = pt * pt + 1;
    }

    /* Build work items */
    work_items = (WorkItem *)malloc(partition_target * num_parts * sizeof(WorkItem));
    num_work_items = 0;
    for (int pi = num_parts - 1; pi >= 0; pi--) {
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
        for (int i = 0; i <= max_t_val; i++)
            mpz_init(tdata[t].buckets[i]);
        pthread_create(&threads[t], NULL, thread_worker, &tdata[t]);
    }

    for (int t = 0; t < num_threads; t++)
        pthread_join(threads[t], NULL);

    /* Merge buckets */
    mpz_t *merged = (mpz_t *)malloc((max_t_val + 1) * sizeof(mpz_t));
    for (int i = 0; i <= max_t_val; i++)
        mpz_init(merged[i]);

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

    /* Multiply by permcount numerator: (sum of scaled partition)! / (LCD * n!) */
    /* For sequences partitioning n directly: sum of partition = n, permcount numer = n! */
    /* For A002785: sum of 2*partition = 2*k = n-1 or n, permcount numer = (2k)! */
    /* For A000171: sum of 4*partition = 4*(n/4) = n or n-1, permcount numer depends */
    mpz_t fact_sum, fact_n, denom;
    mpz_init(fact_sum);
    mpz_init(fact_n);
    mpz_init(denom);

    int sum_scaled = scale * partition_target;
    mpz_fac_ui(fact_sum, sum_scaled);
    mpz_fac_ui(fact_n, n);

    /* result = total * fact_sum / (LCD * fact_n) */
    mpz_mul(total, total, fact_sum);
    mpz_mul(denom, LCD, fact_n);
    mpz_divexact(total, total, denom);

    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = (t_end.tv_sec - t_start.tv_sec) +
                     (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

    gmp_printf("a(%d) = %Zd\n", n, total);
    printf("Time: %.3fs, Partitions: %lld, Threads: %d\n",
           elapsed, total_count, num_threads);

    /* Cleanup */
    mpz_clear(total); mpz_clear(temp); mpz_clear(fact_sum);
    mpz_clear(fact_n); mpz_clear(denom); mpz_clear(LCD);
    for (int i = 0; i <= max_t_val; i++) mpz_clear(merged[i]);
    free(merged);
    for (int t = 0; t < num_threads; t++) {
        for (int i = 0; i <= max_t_val; i++)
            mpz_clear(tdata[t].buckets[i]);
        free(tdata[t].buckets);
    }
    free(tdata); free(threads); free(work_items);
    for (int pi = 0; pi < num_parts; pi++) {
        for (int m = 0; m <= max_mult[pi]; m++)
            mpz_clear(km_fac[pi][m]);
        free(km_fac[pi]); free(gcd_tab[pi]);
    }
    free(km_fac); free(gcd_tab); free(max_mult); free(parts);

    return 0;
}
