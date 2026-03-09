/*
 * Unified Burnside/Polya enumeration for graph-counting sequences, v2.
 *
 * Supports 8 sequences:
 *   A000568 - Tournaments              (odd parts of n, base 2)
 *   A000273 - Digraphs                 (all parts of n, base 2)
 *   A000595 - Binary relations         (all parts of n, base 2)
 *   A000088 - Simple graphs            (all parts of n, base 2)
 *   A000666 - Symmetric relations      (all parts of n, base 2)
 *   A001174 - Oriented graphs          (all parts of n, base 3)
 *   A002785 - Self-comp tournaments    (odd parts of n/2, base 2)
 *   A000171 - Self-comp graphs         (parts of n/4, base 2)
 *
 * For base-2 sequences: uses bucket accumulation (group by t, shift at end).
 * For base-3 sequences: direct multiplication at each leaf.
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o burnside_enum_v2 burnside_enum_v2.c -lgmp -lpthread
 *
 * Usage: ./burnside_enum_v2 <sequence> <n> [threads]
 *   sequence: 568, 273, 595, 88, 666, 1174, 2785, 171
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

typedef enum {
    SEQ_A000568, SEQ_A000273, SEQ_A000595, SEQ_A000088,
    SEQ_A000666, SEQ_A001174, SEQ_A002785, SEQ_A000171,
    SEQ_A003086, SEQ_A005639, SEQ_A002499, SEQ_A002854,
} SeqType;

static SeqType seq_type;
static int global_n, partition_target, odd_parts_only;
static int base_val;  /* 2 or 3 */
static int num_parts;
static int *parts;
static mpz_t **km_fac;
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
    mpz_t *buckets;   /* used for base-2 */
    mpz_t direct_sum; /* used for base-3 */
    long long count;
} ThreadData;

static int gcd_func(int a, int b) {
    while (b) { int t = b; b = a % b; a = t; }
    return a;
}

/*
 * Compute edges in compressed form.
 *
 * All formulas decompose into:
 *   edges = coeff_cross * sum_{i<j_parts} gcd(r_i, r_j)
 *           + sum_r f(r, mr)
 *
 * where f depends on the sequence type.
 */
static long long compute_edges(int *pk, int *pm, int depth) {
    long long gcd_sum = 0;
    long long gcd_sum_even = 0;  /* For A005639: only pairs with ≥1 even cycle */
    long long sum_v = 0, num_v = 0;
    int has_odd_part = 0;

    for (int i = 0; i < depth; i++) {
        int r = parts[pk[i]], mr = pm[i];
        sum_v += (long long)mr * r;
        num_v += mr;
        if (r % 2 == 1) has_odd_part = 1;
        /* Same-part pairs: C(mr,2) * gcd(r,r) = C(mr,2) * r */
        long long same = (long long)mr * (mr - 1) / 2 * r;
        gcd_sum += same;
        if (r % 2 == 0) gcd_sum_even += same;
        for (int j = i + 1; j < depth; j++) {
            long long cross = (long long)mr * pm[j] * gcd_tab[pk[i]][pk[j]];
            gcd_sum += cross;
            if (parts[pk[i]] % 2 == 0 || parts[pk[j]] % 2 == 0)
                gcd_sum_even += cross;
        }
    }

    long long self_term = 0;
    for (int i = 0; i < depth; i++) {
        int r = parts[pk[i]], mr = pm[i];
        switch (seq_type) {
            case SEQ_A000568:
                /* sum floor(v_i/2) → mr * floor(r/2) */
                self_term += (long long)mr * (r / 2);
                break;
            case SEQ_A000088:
                /* sum floor(v_i/2) → same as A000568 */
                self_term += (long long)mr * (r / 2);
                break;
            case SEQ_A000273:
                /* sum(v_i - 1) = sum_v - #v → handled below */
                break;
            case SEQ_A000595:
                /* sum(v_i) = sum_v → handled below */
                break;
            case SEQ_A000666:
                /* sum(floor(v_i/2) + 1) = sum floor(v_i/2) + #v */
                self_term += (long long)mr * (r / 2 + 1);
                break;
            case SEQ_A001174:
                /* sum floor((v_i-1)/2) */
                self_term += (long long)mr * ((r - 1) / 2);
                break;
            case SEQ_A002785:
                /* sum(v_i) → handled with cross coeff */
                break;
            case SEQ_A000171:
                /* 2*sum(v_i) → handled with cross coeff */
                break;
            case SEQ_A003086:
                /* same as A000273: sum(v_i-1) handled below */
                break;
            case SEQ_A005639:
                /* For even r: ((r-2)/4)*2 + 1; for odd r: 0 */
                if (r % 2 == 0)
                    self_term += (long long)mr * ((r - 2) / 4 * 2 + 1);
                break;
            case SEQ_A002499:
                /* v\2 + (if v even: (v-2)\4*2+1, else 0) */
                if (r % 2 == 0)
                    self_term += (long long)mr * (r / 2 + (r - 2) / 4 * 2 + 1);
                else
                    self_term += (long long)mr * (r / 2);  /* = (r-1)/2 */
                break;
            case SEQ_A002854:
                /* floor(r/2) - 1 per cycle */
                self_term += (long long)mr * (r / 2 - 1);
                break;
        }
    }

    long long edge_val;
    switch (seq_type) {
        case SEQ_A000568:
        case SEQ_A000088:
            edge_val = gcd_sum + self_term;
            break;
        case SEQ_A000273:
            edge_val = 2 * gcd_sum + sum_v - num_v;
            break;
        case SEQ_A000595:
            edge_val = 2 * gcd_sum + sum_v;
            break;
        case SEQ_A000666:
            edge_val = gcd_sum + self_term;
            break;
        case SEQ_A001174:
            edge_val = gcd_sum + self_term;
            break;
        case SEQ_A002785:
            edge_val = 2 * gcd_sum + sum_v;
            break;
        case SEQ_A000171:
            edge_val = 4 * gcd_sum + 2 * sum_v;
            break;
        case SEQ_A003086:
            /* 4*sum gcd(p_i,p_j) + sum(2*p_i - 1) with unscaled parts */
            edge_val = 4 * gcd_sum + 2 * sum_v - num_v;
            break;
        case SEQ_A005639:
            /* Only cross pairs with ≥1 even cycle, plus self-term for even cycles */
            edge_val = gcd_sum_even + self_term;
            break;
        case SEQ_A002499:
            /* gcd*2 for even pairs, gcd*1 for odd pairs = gcd_sum + gcd_sum_even */
            edge_val = gcd_sum + gcd_sum_even + self_term;
            break;
        case SEQ_A002854:
            /* gcd_sum + sum(floor(r/2)-1)*m + [has odd part] */
            edge_val = gcd_sum + self_term + has_odd_part;
            break;
        default:
            edge_val = 0;
    }
    return edge_val;
}

static void compute_extra(mpz_t result, int *pm, int depth) {
    mpz_set_ui(result, 1);
    if ((seq_type == SEQ_A002785 || seq_type == SEQ_A000171) && (global_n % 2 == 1)) {
        int total_parts = 0;
        for (int i = 0; i < depth; i++) total_parts += pm[i];
        mpz_set_ui(result, global_n);
        mpz_mul_2exp(result, result, total_parts);  /* n * 2^(#parts) */
    }
    if (seq_type == SEQ_A003086 && (global_n % 2 == 1)) {
        int total_parts = 0;
        for (int i = 0; i < depth; i++) total_parts += pm[i];
        mpz_set_ui(result, global_n);
        mpz_mul_2exp(result, result, 2 * total_parts);  /* n * 4^(#parts) */
    }
}

/* Precomputed powers of 3 for base-3 sequences */
static mpz_t *pow3_table;
static int pow3_max;

static void enumerate(
    ThreadData *td, int remaining, int max_pi,
    int *pk, int *pm, int depth, mpz_t lcd_over_z
) {
    if (remaining == 0) {
        td->count++;
        long long edge_val = compute_edges(pk, pm, depth);
        if (edge_val < 0 || edge_val > max_t_val) {
            fprintf(stderr, "ERROR: edge_val=%lld out of range [0,%d]\n", edge_val, max_t_val);
            return;
        }

        mpz_t extra, contrib;
        mpz_init(extra);
        mpz_init(contrib);
        compute_extra(extra, pm, depth);
        mpz_mul(contrib, lcd_over_z, extra);

        if (base_val == 2) {
            mpz_add(td->buckets[edge_val], td->buckets[edge_val], contrib);
        } else {
            /* base 3: multiply by 3^edge_val directly */
            mpz_t term;
            mpz_init(term);
            if (edge_val <= pow3_max) {
                mpz_mul(term, contrib, pow3_table[edge_val]);
            } else {
                mpz_t p3;
                mpz_init(p3);
                mpz_ui_pow_ui(p3, 3, edge_val);
                mpz_mul(term, contrib, p3);
                mpz_clear(p3);
            }
            mpz_add(td->direct_sum, td->direct_sum, term);
            mpz_clear(term);
        }

        mpz_clear(extra);
        mpz_clear(contrib);
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
    int *pk = (int *)calloc(partition_target + 1, sizeof(int));
    int *pm = (int *)calloc(partition_target + 1, sizeof(int));
    mpz_t lcd_over_z;
    mpz_init(lcd_over_z);

    while (1) {
        int idx = atomic_fetch_add(&work_counter, 1);
        if (idx >= num_work_items) break;

        int pi = work_items[idx].pi;
        int m = work_items[idx].m;

        mpz_divexact(lcd_over_z, LCD, km_fac[pi][m]);
        pk[0] = pi; pm[0] = m;

        int remaining = partition_target - m * parts[pi];
        if (remaining == 0) {
            td->count++;
            long long edge_val = compute_edges(pk, pm, 1);
            if (edge_val >= 0 && edge_val <= max_t_val) {
                mpz_t extra, contrib;
                mpz_init(extra); mpz_init(contrib);
                compute_extra(extra, pm, 1);
                mpz_mul(contrib, lcd_over_z, extra);
                if (base_val == 2) {
                    mpz_add(td->buckets[edge_val], td->buckets[edge_val], contrib);
                } else {
                    mpz_t term; mpz_init(term);
                    if (edge_val <= pow3_max) {
                        mpz_mul(term, contrib, pow3_table[edge_val]);
                    } else {
                        mpz_t p3; mpz_init(p3);
                        mpz_ui_pow_ui(p3, 3, edge_val);
                        mpz_mul(term, contrib, p3);
                        mpz_clear(p3);
                    }
                    mpz_add(td->direct_sum, td->direct_sum, term);
                    mpz_clear(term);
                }
                mpz_clear(extra); mpz_clear(contrib);
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
        fprintf(stderr, "Usage: %s <seq: 568|273|595|88|666|1174|2785|171|3086|5639|2499|2854> <n> [threads]\n", argv[0]);
        return 1;
    }

    int seq_id = atoi(argv[1]);
    global_n = atoi(argv[2]);
    int n = global_n;
    int num_threads = (argc >= 4) ? atoi(argv[3]) : 1;

    base_val = 2;
    int scale = 1;

    switch (seq_id) {
        case 568:  seq_type = SEQ_A000568; partition_target = n; odd_parts_only = 1; break;
        case 273:  seq_type = SEQ_A000273; partition_target = n; odd_parts_only = 0; break;
        case 595:  seq_type = SEQ_A000595; partition_target = n; odd_parts_only = 0; break;
        case 88:   seq_type = SEQ_A000088; partition_target = n; odd_parts_only = 0; break;
        case 666:  seq_type = SEQ_A000666; partition_target = n; odd_parts_only = 0; break;
        case 1174: seq_type = SEQ_A001174; partition_target = n; odd_parts_only = 0; base_val = 3; break;
        case 2785: seq_type = SEQ_A002785; partition_target = n/2; odd_parts_only = 1; scale = 2; break;
        case 171:  seq_type = SEQ_A000171; partition_target = n/4; odd_parts_only = 0; scale = 4;
            if (n % 4 >= 2) { printf("a(%d) = 0\n", n); return 0; }
            break;
        case 3086: seq_type = SEQ_A003086; partition_target = n/2; odd_parts_only = 0; scale = 2; break;
        case 5639: seq_type = SEQ_A005639; partition_target = n; odd_parts_only = 0; base_val = 3; break;
        case 2499: seq_type = SEQ_A002499; partition_target = n; odd_parts_only = 0; break;
        case 2854: seq_type = SEQ_A002854; partition_target = n; odd_parts_only = 0; break;
        default:
            fprintf(stderr, "Unknown sequence: %d\n", seq_id);
            return 1;
    }

    if (partition_target == 0) {
        /* a(0) = 1 for most sequences */
        printf("a(%d) = 1\n", n);
        return 0;
    }

    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    /* Build parts list */
    num_parts = 0;
    parts = (int *)malloc((partition_target + 1) * sizeof(int));
    if (odd_parts_only) {
        for (int i = 1; i <= partition_target; i += 2) parts[num_parts++] = i;
    } else {
        for (int i = 1; i <= partition_target; i++) parts[num_parts++] = i;
    }

    max_mult = (int *)malloc(num_parts * sizeof(int));
    for (int pi = 0; pi < num_parts; pi++)
        max_mult[pi] = partition_target / parts[pi];

    gcd_tab = (int **)malloc(num_parts * sizeof(int *));
    for (int i = 0; i < num_parts; i++) {
        gcd_tab[i] = (int *)malloc(num_parts * sizeof(int));
        for (int j = 0; j < num_parts; j++)
            gcd_tab[i][j] = gcd_func(parts[i], parts[j]);
    }

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

    mpz_init_set_ui(LCD, 1);
    for (int pi = 0; pi < num_parts; pi++)
        mpz_mul(LCD, LCD, km_fac[pi][max_mult[pi]]);

    /* Max t value */
    long long pt = partition_target;
    switch (seq_type) {
        case SEQ_A000568: case SEQ_A000088:
            max_t_val = pt * (pt - 1) / 2 + pt / 2 + 1; break;
        case SEQ_A000273:
            max_t_val = pt * (pt - 1) + pt + 1; break;
        case SEQ_A000595:
            max_t_val = pt * (pt - 1) + pt + 1; break;
        case SEQ_A000666:
            max_t_val = pt * (pt - 1) / 2 + pt + 1; break;
        case SEQ_A001174:
            max_t_val = pt * (pt - 1) / 2 + pt + 1; break;
        case SEQ_A002785:
            max_t_val = pt * pt + pt + 1; break;
        case SEQ_A000171:
            max_t_val = 2 * pt * (pt - 1) + 2 * pt + 1; break;
        case SEQ_A003086:
            /* 4*sum gcd + 2*sum_v - num_v ≤ 4*C(pt,2)*pt + 2*pt - 1 */
            max_t_val = 2 * pt * pt + pt + 1; break;
        case SEQ_A005639:
            /* gcd_sum_even + self_term ≤ pt*(pt-1)/2 + pt */
            max_t_val = pt * (pt - 1) / 2 + pt + 1; break;
        case SEQ_A002499:
            /* gcd_sum + gcd_sum_even + self ≤ 2*C(pt,2)*pt + 2*pt */
            max_t_val = pt * pt + pt + 1; break;
        case SEQ_A002854:
            /* gcd_sum + self + 1 ≤ C(pt,2)*pt + pt/2 + 1 */
            max_t_val = pt * (pt - 1) / 2 + pt / 2 + 2; break;
        default:
            max_t_val = pt * pt + 1;
    }

    /* Precompute powers of 3 if needed */
    pow3_max = 0;
    pow3_table = NULL;
    if (base_val == 3) {
        pow3_max = max_t_val;
        pow3_table = (mpz_t *)malloc((pow3_max + 1) * sizeof(mpz_t));
        mpz_init_set_ui(pow3_table[0], 1);
        for (int i = 1; i <= pow3_max; i++) {
            mpz_init(pow3_table[i]);
            mpz_mul_ui(pow3_table[i], pow3_table[i-1], 3);
        }
    }

    /* Work items */
    work_items = (WorkItem *)malloc(partition_target * num_parts * sizeof(WorkItem));
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
        mpz_init(tdata[t].direct_sum);
        if (base_val == 2) {
            tdata[t].buckets = (mpz_t *)malloc((max_t_val + 1) * sizeof(mpz_t));
            for (int i = 0; i <= max_t_val; i++) mpz_init(tdata[t].buckets[i]);
        } else {
            tdata[t].buckets = NULL;
        }
        pthread_create(&threads[t], NULL, thread_worker, &tdata[t]);
    }

    for (int t = 0; t < num_threads; t++)
        pthread_join(threads[t], NULL);

    /* Merge and assemble */
    mpz_t total, temp;
    mpz_init(total);
    mpz_init(temp);
    long long total_count = 0;

    if (base_val == 2) {
        mpz_t *merged = (mpz_t *)malloc((max_t_val + 1) * sizeof(mpz_t));
        for (int i = 0; i <= max_t_val; i++) mpz_init(merged[i]);

        for (int t = 0; t < num_threads; t++) {
            total_count += tdata[t].count;
            for (int i = 0; i <= max_t_val; i++)
                mpz_add(merged[i], merged[i], tdata[t].buckets[i]);
        }

        for (int t = max_t_val; t >= 0; t--) {
            if (mpz_sgn(merged[t]) != 0) {
                mpz_mul_2exp(temp, merged[t], t);
                mpz_add(total, total, temp);
            }
        }

        for (int i = 0; i <= max_t_val; i++) mpz_clear(merged[i]);
        free(merged);
    } else {
        /* base 3: direct sum */
        for (int t = 0; t < num_threads; t++) {
            total_count += tdata[t].count;
            mpz_add(total, total, tdata[t].direct_sum);
        }
    }

    /* Multiply by permcount numerator and divide */
    int sum_scaled = scale * partition_target;
    mpz_t fact_sum, fact_n, denom;
    mpz_init(fact_sum); mpz_init(fact_n); mpz_init(denom);
    mpz_fac_ui(fact_sum, sum_scaled);
    mpz_fac_ui(fact_n, n);
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
    mpz_clear(total); mpz_clear(temp);
    mpz_clear(fact_sum); mpz_clear(fact_n); mpz_clear(denom); mpz_clear(LCD);
    for (int t = 0; t < num_threads; t++) {
        mpz_clear(tdata[t].direct_sum);
        if (tdata[t].buckets) {
            for (int i = 0; i <= max_t_val; i++) mpz_clear(tdata[t].buckets[i]);
            free(tdata[t].buckets);
        }
    }
    free(tdata); free(threads); free(work_items);
    if (pow3_table) {
        for (int i = 0; i <= pow3_max; i++) mpz_clear(pow3_table[i]);
        free(pow3_table);
    }
    for (int pi = 0; pi < num_parts; pi++) {
        for (int m = 0; m <= max_mult[pi]; m++) mpz_clear(km_fac[pi][m]);
        free(km_fac[pi]); free(gcd_tab[pi]);
    }
    free(km_fac); free(gcd_tab); free(max_mult); free(parts);

    return 0;
}
