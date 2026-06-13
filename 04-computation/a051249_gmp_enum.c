/*
 * A051249: Number of 5-uniform hypergraphs on n unlabeled nodes.
 *
 * Uses closed-form c_5 computation via divisor-signature approach.
 * Key: [x^5] prod_a (1+x^{l_a})^{c_a} is nonzero only when l_a <= 5,
 * so only divisors d with r/d in {1,2,3,4,5} matter.
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o a051249_gmp_enum a051249_gmp_enum.c -lgmp -lpthread
 *
 * Usage: ./a051249_gmp_enum <n> [threads]
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

#define MAX_DEPTH 200
#define K_VAL 5

static int global_n;
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

static int mobius_func(int n) {
    if (n == 1) return 1;
    int factors = 0, m = n;
    for (int p = 2; p * p <= m; p++) {
        if (m % p == 0) {
            m /= p; factors++;
            if (m % p == 0) return 0;
        }
    }
    if (m > 1) factors++;
    return (factors & 1) ? -1 : 1;
}

static int get_divisors(int n, int *divs, int max_divs) {
    int cnt = 0;
    for (int i = 1; i * i <= n; i++) {
        if (n % i == 0) {
            if (cnt < max_divs) divs[cnt++] = i;
            if (i != n / i && cnt < max_divs) divs[cnt++] = n / i;
        }
    }
    for (int i = 0; i < cnt - 1; i++)
        for (int j = i + 1; j < cnt; j++)
            if (divs[i] > divs[j]) { int t = divs[i]; divs[i] = divs[j]; divs[j] = t; }
    return cnt;
}

static long long count_gcd_pattern(int *rs, int *ds, int p, long long L) {
    int qs[5];
    int div_lists[5][64];
    int div_counts[5];
    for (int i = 0; i < p; i++) {
        qs[i] = rs[i] / ds[i];
        div_counts[i] = get_divisors(qs[i], div_lists[i], 64);
    }

    long long total = 0;

    /* Unrolled for p = 1..5 */
    if (p == 1) {
        for (int i0 = 0; i0 < div_counts[0]; i0++) {
            int mu0 = mobius_func(div_lists[0][i0]);
            if (!mu0) continue;
            long long d0e0 = (long long)ds[0] * div_lists[0][i0];
            if (L % d0e0) continue;
            total += mu0 * (L / d0e0);
        }
    } else if (p == 2) {
        for (int i0 = 0; i0 < div_counts[0]; i0++) {
            int mu0 = mobius_func(div_lists[0][i0]);
            if (!mu0) continue;
            long long d0e0 = (long long)ds[0] * div_lists[0][i0];
            for (int i1 = 0; i1 < div_counts[1]; i1++) {
                int mu1 = mobius_func(div_lists[1][i1]);
                if (!mu1) continue;
                long long lc = lcm_func(d0e0, (long long)ds[1] * div_lists[1][i1]);
                if (lc > L || L % lc) continue;
                total += mu0 * mu1 * (L / lc);
            }
        }
    } else if (p == 3) {
        for (int i0 = 0; i0 < div_counts[0]; i0++) {
            int mu0 = mobius_func(div_lists[0][i0]);
            if (!mu0) continue;
            long long d0e0 = (long long)ds[0] * div_lists[0][i0];
            for (int i1 = 0; i1 < div_counts[1]; i1++) {
                int mu1 = mobius_func(div_lists[1][i1]);
                if (!mu1) continue;
                long long lc01 = lcm_func(d0e0, (long long)ds[1] * div_lists[1][i1]);
                if (lc01 > L) continue;
                for (int i2 = 0; i2 < div_counts[2]; i2++) {
                    int mu2 = mobius_func(div_lists[2][i2]);
                    if (!mu2) continue;
                    long long lc = lcm_func(lc01, (long long)ds[2] * div_lists[2][i2]);
                    if (lc > L || L % lc) continue;
                    total += mu0 * mu1 * mu2 * (L / lc);
                }
            }
        }
    } else if (p == 4) {
        for (int i0 = 0; i0 < div_counts[0]; i0++) {
            int mu0 = mobius_func(div_lists[0][i0]);
            if (!mu0) continue;
            long long d0e0 = (long long)ds[0] * div_lists[0][i0];
            for (int i1 = 0; i1 < div_counts[1]; i1++) {
                int mu1 = mobius_func(div_lists[1][i1]);
                if (!mu1) continue;
                long long lc01 = lcm_func(d0e0, (long long)ds[1] * div_lists[1][i1]);
                if (lc01 > L) continue;
                for (int i2 = 0; i2 < div_counts[2]; i2++) {
                    int mu2 = mobius_func(div_lists[2][i2]);
                    if (!mu2) continue;
                    long long lc012 = lcm_func(lc01, (long long)ds[2] * div_lists[2][i2]);
                    if (lc012 > L) continue;
                    for (int i3 = 0; i3 < div_counts[3]; i3++) {
                        int mu3 = mobius_func(div_lists[3][i3]);
                        if (!mu3) continue;
                        long long lc = lcm_func(lc012, (long long)ds[3] * div_lists[3][i3]);
                        if (lc > L || L % lc) continue;
                        total += mu0 * mu1 * mu2 * mu3 * (L / lc);
                    }
                }
            }
        }
    } else { /* p == 5 */
        for (int i0 = 0; i0 < div_counts[0]; i0++) {
            int mu0 = mobius_func(div_lists[0][i0]);
            if (!mu0) continue;
            long long d0e0 = (long long)ds[0] * div_lists[0][i0];
            for (int i1 = 0; i1 < div_counts[1]; i1++) {
                int mu1 = mobius_func(div_lists[1][i1]);
                if (!mu1) continue;
                long long lc01 = lcm_func(d0e0, (long long)ds[1] * div_lists[1][i1]);
                if (lc01 > L) continue;
                for (int i2 = 0; i2 < div_counts[2]; i2++) {
                    int mu2 = mobius_func(div_lists[2][i2]);
                    if (!mu2) continue;
                    long long lc012 = lcm_func(lc01, (long long)ds[2] * div_lists[2][i2]);
                    if (lc012 > L) continue;
                    for (int i3 = 0; i3 < div_counts[3]; i3++) {
                        int mu3 = mobius_func(div_lists[3][i3]);
                        if (!mu3) continue;
                        long long lc0123 = lcm_func(lc012, (long long)ds[3] * div_lists[3][i3]);
                        if (lc0123 > L) continue;
                        for (int i4 = 0; i4 < div_counts[4]; i4++) {
                            int mu4 = mobius_func(div_lists[4][i4]);
                            if (!mu4) continue;
                            long long lc = lcm_func(lc0123, (long long)ds[4] * div_lists[4][i4]);
                            if (lc > L || L % lc) continue;
                            total += (long long)mu0 * mu1 * mu2 * mu3 * mu4 * (L / lc);
                        }
                    }
                }
            }
        }
    }

    return total;
}

/*
 * Enumerate all compositions of K_VAL into p positive parts (sorted decreasing)
 * and their distinct permutations, then assign to combinations of p parts from depth.
 *
 * Generic approach: iterate over all (j_0,...,j_{p-1}) with sum = K and j_i > 0.
 * For each, iterate over all C(depth, p) choices of active parts.
 * For each, iterate over all (d_val) choices per active part.
 */

/* Precomputed compositions of 5 into p parts */
/* p=1: {5}
 * p=2: {4,1},{3,2},{2,3},{1,4}  -- all ordered
 * p=3: {3,1,1},{1,3,1},{1,1,3},{2,2,1},{2,1,2},{1,2,2}
 * p=4: {2,1,1,1},{1,2,1,1},{1,1,2,1},{1,1,1,2}
 * p=5: {1,1,1,1,1}
 *
 * Actually for efficiency, enumerate sorted compositions and then
 * generate distinct permutations. But in C this is complex.
 * Since K=5 and max p=5, the total number of ordered compositions is small.
 * compositions(5, p) for p=1..5: 1, 4, 6, 4, 1 = 16 total.
 * Just enumerate them directly.
 */

static int all_comps[][5] = {
    /* p=1 */
    {5, 0, 0, 0, 0},
    /* p=2 */
    {4, 1, 0, 0, 0}, {3, 2, 0, 0, 0}, {2, 3, 0, 0, 0}, {1, 4, 0, 0, 0},
    /* p=3 */
    {3, 1, 1, 0, 0}, {1, 3, 1, 0, 0}, {1, 1, 3, 0, 0},
    {2, 2, 1, 0, 0}, {2, 1, 2, 0, 0}, {1, 2, 2, 0, 0},
    /* p=4 */
    {2, 1, 1, 1, 0}, {1, 2, 1, 1, 0}, {1, 1, 2, 1, 0}, {1, 1, 1, 2, 0},
    /* p=5 */
    {1, 1, 1, 1, 1}
};
static int comp_sizes[] = {
    1, 1, 1, 1, 1,
    2, 2, 2, 2, 2, 2,
    3, 3, 3, 3, 3, 3,
    4, 4, 4, 4,
    5
};
/* Map: which compositions have p parts */
static int comp_p[] = {
    1,
    2, 2, 2, 2,
    3, 3, 3, 3, 3, 3,
    4, 4, 4, 4,
    5
};
#define NUM_COMPS 16

/* Choose p indices from depth (combinations) */
static void choose_indices(int *indices, int depth, int p, int start, int cur,
                          int *pk, int *pm, long long L,
                          int *js, long long *total_ptr) {
    if (cur == p) {
        /* For each active index, try all valid (d_val, l) options */
        /* Recursive enumeration of options */
        int nopts[5];
        int opt_d_arr[5][5]; /* max 5 options per part (l=1..5) */
        long long opt_bv_arr[5][5];
        int rs_arr[5], ds_arr[5];

        for (int i = 0; i < p; i++) {
            int pi = pk[indices[i]];
            int r = parts_list[pi];
            int m = pm[indices[i]];
            int j = js[i];
            nopts[i] = 0;
            for (int l = 1; l <= K_VAL; l++) {
                if (r % l != 0) continue;
                if (j % l != 0) continue;
                int d_val = r / l;
                int nc = j / l;
                int avail = m * d_val;
                if (nc > avail) continue;
                opt_d_arr[i][nopts[i]] = d_val;
                opt_bv_arr[i][nopts[i]] = binomial(avail, nc);
                nopts[i]++;
            }
            rs_arr[i] = r;
            if (nopts[i] == 0) return;
        }

        /* Enumerate all option combinations */
        /* Unroll for p=1..5 */
        if (p == 1) {
            for (int o0 = 0; o0 < nopts[0]; o0++) {
                ds_arr[0] = opt_d_arr[0][o0];
                long long cnt = count_gcd_pattern(rs_arr, ds_arr, 1, L);
                *total_ptr += cnt * opt_bv_arr[0][o0];
            }
        } else if (p == 2) {
            for (int o0 = 0; o0 < nopts[0]; o0++) {
                ds_arr[0] = opt_d_arr[0][o0];
                for (int o1 = 0; o1 < nopts[1]; o1++) {
                    ds_arr[1] = opt_d_arr[1][o1];
                    long long cnt = count_gcd_pattern(rs_arr, ds_arr, 2, L);
                    *total_ptr += cnt * opt_bv_arr[0][o0] * opt_bv_arr[1][o1];
                }
            }
        } else if (p == 3) {
            for (int o0 = 0; o0 < nopts[0]; o0++) {
                ds_arr[0] = opt_d_arr[0][o0];
                for (int o1 = 0; o1 < nopts[1]; o1++) {
                    ds_arr[1] = opt_d_arr[1][o1];
                    for (int o2 = 0; o2 < nopts[2]; o2++) {
                        ds_arr[2] = opt_d_arr[2][o2];
                        long long cnt = count_gcd_pattern(rs_arr, ds_arr, 3, L);
                        *total_ptr += cnt * opt_bv_arr[0][o0] * opt_bv_arr[1][o1] * opt_bv_arr[2][o2];
                    }
                }
            }
        } else if (p == 4) {
            for (int o0 = 0; o0 < nopts[0]; o0++) {
                ds_arr[0] = opt_d_arr[0][o0];
                for (int o1 = 0; o1 < nopts[1]; o1++) {
                    ds_arr[1] = opt_d_arr[1][o1];
                    for (int o2 = 0; o2 < nopts[2]; o2++) {
                        ds_arr[2] = opt_d_arr[2][o2];
                        for (int o3 = 0; o3 < nopts[3]; o3++) {
                            ds_arr[3] = opt_d_arr[3][o3];
                            long long cnt = count_gcd_pattern(rs_arr, ds_arr, 4, L);
                            *total_ptr += cnt * opt_bv_arr[0][o0] * opt_bv_arr[1][o1]
                                         * opt_bv_arr[2][o2] * opt_bv_arr[3][o3];
                        }
                    }
                }
            }
        } else { /* p == 5 */
            for (int o0 = 0; o0 < nopts[0]; o0++) {
                ds_arr[0] = opt_d_arr[0][o0];
                for (int o1 = 0; o1 < nopts[1]; o1++) {
                    ds_arr[1] = opt_d_arr[1][o1];
                    for (int o2 = 0; o2 < nopts[2]; o2++) {
                        ds_arr[2] = opt_d_arr[2][o2];
                        for (int o3 = 0; o3 < nopts[3]; o3++) {
                            ds_arr[3] = opt_d_arr[3][o3];
                            for (int o4 = 0; o4 < nopts[4]; o4++) {
                                ds_arr[4] = opt_d_arr[4][o4];
                                long long cnt = count_gcd_pattern(rs_arr, ds_arr, 5, L);
                                *total_ptr += cnt * opt_bv_arr[0][o0] * opt_bv_arr[1][o1]
                                             * opt_bv_arr[2][o2] * opt_bv_arr[3][o3] * opt_bv_arr[4][o4];
                            }
                        }
                    }
                }
            }
        }
        return;
    }

    for (int i = start; i < depth; i++) {
        indices[cur] = i;
        choose_indices(indices, depth, p, i + 1, cur + 1, pk, pm, L, js, total_ptr);
    }
}

static long long compute_ck(int *pk, int *pm, int depth) {
    long long L = 1;
    for (int i = 0; i < depth; i++)
        L = lcm_func(L, parts_list[pk[i]]);

    long long total = 0;
    int indices[5];

    for (int ci = 0; ci < NUM_COMPS; ci++) {
        int p = comp_p[ci];
        if (p > depth) continue;
        choose_indices(indices, depth, p, 0, 0, pk, pm, L, all_comps[ci], &total);
    }

    if (total % L != 0) {
        fprintf(stderr, "ERROR: ck not integer: %lld/%lld\n", total, L);
        return -1;
    }
    return total / L;
}


static void enumerate(
    ThreadData *td, int remaining, int max_pi,
    int *pk, int *pm, int depth_val, mpz_t lcd_over_z
) {
    if (remaining == 0) {
        td->count++;
        long long ck = compute_ck(pk, pm, depth_val);
        if (ck < 0 || ck > max_t_val) {
            fprintf(stderr, "ERROR: ck=%lld out of range\n", ck);
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

            pk[depth_val] = pi;
            pm[depth_val] = m;
            enumerate(td, remaining - m * k, pi - 1, pk, pm, depth_val + 1, cur_lcd);
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
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n> [threads]\n", argv[0]);
        return 1;
    }

    global_n = atoi(argv[1]);
    int n = global_n;
    int num_threads = (argc >= 3) ? atoi(argv[2]) : 1;

    if (n < K_VAL) {
        printf("a(%d) = 1\n", n);
        return 0;
    }

    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    num_parts = n;
    parts_list = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) parts_list[i] = i + 1;

    max_mult = (int *)malloc(num_parts * sizeof(int));
    for (int pi = 0; pi < num_parts; pi++)
        max_mult[pi] = n / parts_list[pi];

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

    /* Max c_5 value: C(n, 5) */
    long long cnk = (long long)n * (n-1) * (n-2) * (n-3) * (n-4) / 120;
    max_t_val = (int)(cnk > 50000000 ? 50000000 : cnk) + 1;

    work_items = (WorkItem *)malloc(n * num_parts * sizeof(WorkItem));
    num_work_items = 0;
    for (int pi = num_parts - 1; pi >= 0; pi--)
        for (int m = 1; m <= max_mult[pi]; m++) {
            work_items[num_work_items].pi = pi;
            work_items[num_work_items].m = m;
            num_work_items++;
        }
    atomic_store(&work_counter, 0);

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

    mpz_t *merged = (mpz_t *)malloc((max_t_val + 1) * sizeof(mpz_t));
    for (int i = 0; i <= max_t_val; i++) mpz_init(merged[i]);
    long long total_count = 0;

    for (int t = 0; t < num_threads; t++) {
        total_count += tdata[t].count;
        for (int i = 0; i <= max_t_val; i++)
            mpz_add(merged[i], merged[i], tdata[t].buckets[i]);
    }

    mpz_t total, temp;
    mpz_init(total);
    mpz_init(temp);
    for (int t = max_t_val; t >= 0; t--) {
        if (mpz_sgn(merged[t]) != 0) {
            mpz_mul_2exp(temp, merged[t], t);
            mpz_add(total, total, temp);
        }
    }

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
        free(km_fac[pi]);
    }
    free(km_fac); free(max_mult); free(parts_list);

    return 0;
}
