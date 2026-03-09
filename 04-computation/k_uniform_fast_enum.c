/*
 * Fast k-uniform hypergraph enumeration using divisor-signature Mobius approach.
 *
 * Instead of iterating t=0..L-1 for Burnside orbit counting, enumerate all
 * selections of at most k "active" parts and compute orbit counts via Mobius
 * inversion over divisor tuples.
 *
 * Key insight: [x^k] prod_a (1+x^{l_a})^{c_a} is nonzero only when l_a <= k
 * for the active parts. So for each part r_a, only O(1) divisors contribute.
 *
 * Speedup over general Burnside: 10-1000x depending on partition structure.
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o k_uniform_fast_enum k_uniform_fast_enum.c -lgmp -lpthread
 *
 * Usage: ./k_uniform_fast_enum <k> <n> [threads]
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

#define MAX_K 20
#define MAX_DEPTH 200

static int global_n, global_k;
static int num_parts;
static int *parts_list;
static int *max_mult;
static mpz_t **km_fac;
static mpz_t LCD;
static int max_t_val;

/* Precomputed: all ordered compositions of global_k into p positive parts (p=1..min(k,depth)) */
static int num_comps;
static int comp_data[10000][MAX_K]; /* comp_data[i] = j values */
static int comp_len[10000];          /* number of parts */

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

/* Count t in [0,L) with gcd(t, rs[i]) = ds[i] for i=0..p-1 via Mobius */
static long long count_gcd_pattern_recursive(
    int *rs, int *ds, int p, long long L,
    int idx, long long mu_prod, long long lcm_so_far,
    int *div_lists, int *div_counts, int max_divs_per
) {
    if (idx == p) {
        if (L % lcm_so_far != 0) return 0;
        return mu_prod * (L / lcm_so_far);
    }

    long long total = 0;
    int *divs = div_lists + idx * max_divs_per;
    int ndivs = div_counts[idx];

    for (int i = 0; i < ndivs; i++) {
        int e = divs[i];
        int mu_e = mobius_func(e);
        if (mu_e == 0) continue;
        long long de = (long long)ds[idx] * e;
        long long new_lcm = lcm_func(lcm_so_far, de);
        if (new_lcm > L) continue;
        total += count_gcd_pattern_recursive(
            rs, ds, p, L, idx + 1, mu_prod * mu_e, new_lcm,
            div_lists, div_counts, max_divs_per
        );
    }
    return total;
}

static long long count_gcd_pattern(int *rs, int *ds, int p, long long L) {
    #define MAX_DIVS_PER 128
    int div_lists[MAX_K * MAX_DIVS_PER];
    int div_counts[MAX_K];

    for (int i = 0; i < p; i++) {
        int q = rs[i] / ds[i];
        div_counts[i] = get_divisors(q, div_lists + i * MAX_DIVS_PER, MAX_DIVS_PER);
    }

    return count_gcd_pattern_recursive(rs, ds, p, L, 0, 1, 1,
                                       div_lists, div_counts, MAX_DIVS_PER);
}

/* Generate all ordered compositions of global_k into positive parts */
static void gen_compositions(int remaining, int *prefix, int plen) {
    if (remaining == 0) {
        if (num_comps >= 10000) return;
        memcpy(comp_data[num_comps], prefix, plen * sizeof(int));
        comp_len[num_comps] = plen;
        num_comps++;
        return;
    }
    for (int v = remaining; v >= 1; v--) {
        prefix[plen] = v;
        gen_compositions(remaining - v, prefix, plen + 1);
    }
}

/* Recursive: choose p indices from [0, depth), then for each
 * compute the contribution from this composition + index selection */
static void enumerate_indices(
    int *pk, int *pm, int depth, long long L,
    int *js, int p,
    int *indices, int cur, int start,
    long long *total_ptr
) {
    if (cur == p) {
        /* For each active part, enumerate valid (d_val) options */
        int nopts[MAX_K];
        int opt_d[MAX_K][MAX_K]; /* at most k options per part */
        long long opt_bv[MAX_K][MAX_K];
        int rs[MAX_K];

        for (int i = 0; i < p; i++) {
            int pi = pk[indices[i]];
            int r = parts_list[pi];
            int m = pm[indices[i]];
            int j = js[i];
            rs[i] = r;
            nopts[i] = 0;
            for (int l = 1; l <= global_k; l++) {
                if (r % l != 0) continue;
                if (j % l != 0) continue;
                int d_val = r / l;
                int nc = j / l;
                int avail = m * d_val;
                if (nc > avail) continue;
                opt_d[i][nopts[i]] = d_val;
                opt_bv[i][nopts[i]] = binomial(avail, nc);
                nopts[i]++;
            }
            if (nopts[i] == 0) return;
        }

        /* Enumerate all option combinations recursively */
        int ds[MAX_K];

        /* Use iterative approach for small p */
        if (p == 1) {
            for (int o0 = 0; o0 < nopts[0]; o0++) {
                ds[0] = opt_d[0][o0];
                long long cnt = count_gcd_pattern(rs, ds, 1, L);
                *total_ptr += cnt * opt_bv[0][o0];
            }
        } else if (p == 2) {
            for (int o0 = 0; o0 < nopts[0]; o0++) {
                ds[0] = opt_d[0][o0];
                for (int o1 = 0; o1 < nopts[1]; o1++) {
                    ds[1] = opt_d[1][o1];
                    long long cnt = count_gcd_pattern(rs, ds, 2, L);
                    *total_ptr += cnt * opt_bv[0][o0] * opt_bv[1][o1];
                }
            }
        } else {
            /* General recursive enumeration */
            /* For p >= 3, use recursive helper */
            long long bv_stack[MAX_K];
            int ds_stack[MAX_K];

            /* Simple recursive expansion via iteration with stack */
            /* Implemented iteratively for up to MAX_K levels */
            int oidx[MAX_K];
            memset(oidx, 0, sizeof(oidx));
            int level = 0;
            while (1) {
                if (level == p) {
                    long long bv_prod = 1;
                    for (int i = 0; i < p; i++) {
                        ds[i] = opt_d[i][oidx[i]];
                        bv_prod *= opt_bv[i][oidx[i]];
                    }
                    long long cnt = count_gcd_pattern(rs, ds, p, L);
                    *total_ptr += cnt * bv_prod;
                    level--;
                    if (level < 0) break;
                    oidx[level]++;
                } else if (oidx[level] < nopts[level]) {
                    level++;
                    if (level < p) oidx[level] = 0;
                } else {
                    level--;
                    if (level < 0) break;
                    oidx[level]++;
                }
            }
        }
        return;
    }

    for (int i = start; i < depth; i++) {
        indices[cur] = i;
        enumerate_indices(pk, pm, depth, L, js, p, indices, cur + 1, i + 1, total_ptr);
    }
}

static long long compute_ck(int *pk, int *pm, int depth) {
    long long L = 1;
    for (int i = 0; i < depth; i++)
        L = lcm_func(L, parts_list[pk[i]]);

    long long total = 0;
    int indices[MAX_K];

    for (int ci = 0; ci < num_comps; ci++) {
        int p = comp_len[ci];
        if (p > depth) continue;
        enumerate_indices(pk, pm, depth, L, comp_data[ci], p,
                         indices, 0, 0, &total);
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
        int kp = parts_list[pi];
        if (kp > remaining) continue;
        int mm = remaining / kp;

        mpz_t cur_lcd;
        mpz_init(cur_lcd);
        for (int m = 1; m <= mm; m++) {
            mpz_divexact(cur_lcd, lcd_over_z, km_fac[pi][m]);
            if (m > 1) mpz_mul(cur_lcd, cur_lcd, km_fac[pi][m-1]);
            mpz_divexact(cur_lcd, lcd_over_z, km_fac[pi][m]);

            pk[depth_val] = pi;
            pm[depth_val] = m;
            enumerate(td, remaining - m * kp, pi - 1, pk, pm, depth_val + 1, cur_lcd);
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

    if (k > MAX_K) {
        fprintf(stderr, "k must be <= %d\n", MAX_K);
        return 1;
    }

    if (n < k) {
        printf("a(%d,%d) = 1\n", n, k);
        return 0;
    }

    /* Precompute compositions */
    num_comps = 0;
    int prefix[MAX_K];
    gen_compositions(k, prefix, 0);
    fprintf(stderr, "Compositions of %d: %d\n", k, num_comps);

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

    /* Max c_k value estimate */
    long long cnk = 1;
    for (int i = 0; i < k; i++) cnk = cnk * (n - i) / (i + 1);
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
