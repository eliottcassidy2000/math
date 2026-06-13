/*
 * A051240: Number of 4-uniform hypergraphs on n unlabeled nodes.
 *
 * Uses closed-form c_4 computation via divisor-signature approach:
 * instead of iterating t=0..L-1, enumerate active part selections and
 * compute orbit counts via Mobius inversion.
 *
 * Key: [x^4] prod_a (1+x^{l_a})^{c_a} is nonzero only when l_a <= 4,
 * so for each part r_a, only divisors d with r_a/d in {1,2,3,4} matter.
 *
 * Compile: gcc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *          -o a051240_gmp_enum a051240_gmp_enum.c -lgmp -lpthread
 *
 * Usage: ./a051240_gmp_enum <n> [threads]
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

#define MAX_PARTS 200
#define MAX_DEPTH 200
#define K_VAL 4

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

/* Mobius function */
static int mobius_func(int n) {
    if (n == 1) return 1;
    int factors = 0;
    int m = n;
    for (int p = 2; p * p <= m; p++) {
        if (m % p == 0) {
            m /= p;
            factors++;
            if (m % p == 0) return 0;
        }
    }
    if (m > 1) factors++;
    return (factors & 1) ? -1 : 1;
}

/* Divisors of n, stored in divs array, return count */
static int get_divisors(int n, int *divs, int max_divs) {
    int cnt = 0;
    for (int i = 1; i * i <= n; i++) {
        if (n % i == 0) {
            if (cnt < max_divs) divs[cnt++] = i;
            if (i != n / i && cnt < max_divs) divs[cnt++] = n / i;
        }
    }
    /* Sort */
    for (int i = 0; i < cnt - 1; i++)
        for (int j = i + 1; j < cnt; j++)
            if (divs[i] > divs[j]) { int t = divs[i]; divs[i] = divs[j]; divs[j] = t; }
    return cnt;
}

/*
 * Count t in [0,L) with gcd(t, rs[i]) = ds[i] for i=0..p-1.
 * Uses Mobius inversion.
 */
static long long count_gcd_pattern(int *rs, int *ds, int p, long long L) {
    /* q[i] = r[i] / d[i] */
    int qs[4];
    int div_lists[4][64];
    int div_counts[4];

    for (int i = 0; i < p; i++) {
        qs[i] = rs[i] / ds[i];
        div_counts[i] = get_divisors(qs[i], div_lists[i], 64);
    }

    long long total = 0;

    /* Iterate over divisor tuples using nested loops (p <= 4) */
    if (p == 1) {
        for (int i0 = 0; i0 < div_counts[0]; i0++) {
            int e0 = div_lists[0][i0];
            int mu0 = mobius_func(e0);
            if (mu0 == 0) continue;
            long long d0e0 = (long long)ds[0] * e0;
            if (L % d0e0 != 0) continue;
            total += mu0 * (L / d0e0);
        }
    } else if (p == 2) {
        for (int i0 = 0; i0 < div_counts[0]; i0++) {
            int e0 = div_lists[0][i0];
            int mu0 = mobius_func(e0);
            if (mu0 == 0) continue;
            long long d0e0 = (long long)ds[0] * e0;
            for (int i1 = 0; i1 < div_counts[1]; i1++) {
                int e1 = div_lists[1][i1];
                int mu1 = mobius_func(e1);
                if (mu1 == 0) continue;
                long long d1e1 = (long long)ds[1] * e1;
                long long lc = lcm_func(d0e0, d1e1);
                if (lc > L || L % lc != 0) continue;
                total += mu0 * mu1 * (L / lc);
            }
        }
    } else if (p == 3) {
        for (int i0 = 0; i0 < div_counts[0]; i0++) {
            int e0 = div_lists[0][i0];
            int mu0 = mobius_func(e0);
            if (mu0 == 0) continue;
            long long d0e0 = (long long)ds[0] * e0;
            for (int i1 = 0; i1 < div_counts[1]; i1++) {
                int e1 = div_lists[1][i1];
                int mu1 = mobius_func(e1);
                if (mu1 == 0) continue;
                long long lc01 = lcm_func(d0e0, (long long)ds[1] * e1);
                if (lc01 > L) continue;
                for (int i2 = 0; i2 < div_counts[2]; i2++) {
                    int e2 = div_lists[2][i2];
                    int mu2 = mobius_func(e2);
                    if (mu2 == 0) continue;
                    long long lc = lcm_func(lc01, (long long)ds[2] * e2);
                    if (lc > L || L % lc != 0) continue;
                    total += mu0 * mu1 * mu2 * (L / lc);
                }
            }
        }
    } else { /* p == 4 */
        for (int i0 = 0; i0 < div_counts[0]; i0++) {
            int e0 = div_lists[0][i0];
            int mu0 = mobius_func(e0);
            if (mu0 == 0) continue;
            long long d0e0 = (long long)ds[0] * e0;
            for (int i1 = 0; i1 < div_counts[1]; i1++) {
                int e1 = div_lists[1][i1];
                int mu1 = mobius_func(e1);
                if (mu1 == 0) continue;
                long long lc01 = lcm_func(d0e0, (long long)ds[1] * e1);
                if (lc01 > L) continue;
                for (int i2 = 0; i2 < div_counts[2]; i2++) {
                    int e2 = div_lists[2][i2];
                    int mu2 = mobius_func(e2);
                    if (mu2 == 0) continue;
                    long long lc012 = lcm_func(lc01, (long long)ds[2] * e2);
                    if (lc012 > L) continue;
                    for (int i3 = 0; i3 < div_counts[3]; i3++) {
                        int e3 = div_lists[3][i3];
                        int mu3 = mobius_func(e3);
                        if (mu3 == 0) continue;
                        long long lc = lcm_func(lc012, (long long)ds[3] * e3);
                        if (lc > L || L % lc != 0) continue;
                        total += mu0 * mu1 * mu2 * mu3 * (L / lc);
                    }
                }
            }
        }
    }

    return total;
}

/*
 * For a compressed partition pk[0..depth-1] (part indices) and pm[0..depth-1] (multiplicities),
 * compute c_4 = number of orbits on 4-subsets.
 *
 * Enumerates all selections of at most 4 "active" parts and valid (l, j) assignments.
 */
static long long compute_c4(int *pk, int *pm, int depth) {
    long long L = 1;
    for (int i = 0; i < depth; i++)
        L = lcm_func(L, parts_list[pk[i]]);

    long long total = 0;

    /* For each part, compute valid (d_val, l) options where l in {1,2,3,4} and l|r */
    /* options[i][j] = (d_val, l) */
    int nopts[MAX_DEPTH];
    int opt_d[MAX_DEPTH][4];
    int opt_l[MAX_DEPTH][4];

    for (int i = 0; i < depth; i++) {
        int r = parts_list[pk[i]];
        nopts[i] = 0;
        for (int l = 1; l <= 4; l++) {
            if (r % l == 0) {
                opt_d[i][nopts[i]] = r / l;
                opt_l[i][nopts[i]] = l;
                nopts[i]++;
            }
        }
    }

    /* Compositions of 4 into positive parts, with parts sorted decreasingly.
     * For k=4: (4), (3,1), (2,2), (2,1,1), (1,1,1,1)
     * For each composition, generate all distinct permutations.
     */
    /* p=1: (4) */
    for (int a = 0; a < depth; a++) {
        for (int oa = 0; oa < nopts[a]; oa++) {
            if (4 % opt_l[a][oa] != 0) continue;
            int nc = 4 / opt_l[a][oa];
            int avail = pm[a] * opt_d[a][oa];
            if (nc > avail) continue;
            long long bv = binomial(avail, nc);
            int rs[1] = {parts_list[pk[a]]};
            int ds[1] = {opt_d[a][oa]};
            long long cnt = count_gcd_pattern(rs, ds, 1, L);
            total += cnt * bv;
        }
    }

    /* p=2 */
    /* Compositions: (3,1), (2,2), (1,3) — but (1,3) is just (3,1) with parts swapped */
    /* (3,1) and all permutations: (3,1), (1,3) */
    for (int a = 0; a < depth; a++) {
        for (int b = a + 1; b < depth; b++) {
            /* For each pair (a,b), try all (ja, jb) with ja+jb=4, ja>0, jb>0 */
            int j_pairs[3][2] = {{3,1}, {2,2}, {1,3}};
            for (int jp = 0; jp < 3; jp++) {
                int ja = j_pairs[jp][0], jb = j_pairs[jp][1];
                for (int oa = 0; oa < nopts[a]; oa++) {
                    if (ja % opt_l[a][oa] != 0) continue;
                    int nca = ja / opt_l[a][oa];
                    int avail_a = pm[a] * opt_d[a][oa];
                    if (nca > avail_a) continue;
                    long long bva = binomial(avail_a, nca);

                    for (int ob = 0; ob < nopts[b]; ob++) {
                        if (jb % opt_l[b][ob] != 0) continue;
                        int ncb = jb / opt_l[b][ob];
                        int avail_b = pm[b] * opt_d[b][ob];
                        if (ncb > avail_b) continue;
                        long long bvb = binomial(avail_b, ncb);

                        int rs[2] = {parts_list[pk[a]], parts_list[pk[b]]};
                        int ds[2] = {opt_d[a][oa], opt_d[b][ob]};
                        long long cnt = count_gcd_pattern(rs, ds, 2, L);
                        total += cnt * bva * bvb;
                    }
                }
            }
        }
    }

    /* p=3: compositions (2,1,1) and all distinct permutations */
    for (int a = 0; a < depth; a++) {
        for (int b = a + 1; b < depth; b++) {
            for (int c = b + 1; c < depth; c++) {
                /* Permutations of (2,1,1): (2,1,1),(1,2,1),(1,1,2) */
                int perms[3][3] = {{2,1,1}, {1,2,1}, {1,1,2}};
                for (int pp = 0; pp < 3; pp++) {
                    int js[3];
                    int parts_idx[3] = {a, b, c};
                    js[0] = perms[pp][0];
                    js[1] = perms[pp][1];
                    js[2] = perms[pp][2];

                    int valid = 1;
                    long long bv_prod = 1;
                    int rs[3], ds[3];

                    /* For each active part, find best option */
                    /* Actually need to enumerate all option combinations */
                    /* But since each part has at most 4 options, 4^3 = 64 max */
                    for (int o0 = 0; o0 < nopts[parts_idx[0]] && valid; o0++) {
                        if (js[0] % opt_l[parts_idx[0]][o0] != 0) continue;
                        int nc0 = js[0] / opt_l[parts_idx[0]][o0];
                        int av0 = pm[parts_idx[0]] * opt_d[parts_idx[0]][o0];
                        if (nc0 > av0) continue;
                        long long bv0 = binomial(av0, nc0);
                        rs[0] = parts_list[pk[parts_idx[0]]];
                        ds[0] = opt_d[parts_idx[0]][o0];

                        for (int o1 = 0; o1 < nopts[parts_idx[1]]; o1++) {
                            if (js[1] % opt_l[parts_idx[1]][o1] != 0) continue;
                            int nc1 = js[1] / opt_l[parts_idx[1]][o1];
                            int av1 = pm[parts_idx[1]] * opt_d[parts_idx[1]][o1];
                            if (nc1 > av1) continue;
                            long long bv1 = binomial(av1, nc1);
                            rs[1] = parts_list[pk[parts_idx[1]]];
                            ds[1] = opt_d[parts_idx[1]][o1];

                            for (int o2 = 0; o2 < nopts[parts_idx[2]]; o2++) {
                                if (js[2] % opt_l[parts_idx[2]][o2] != 0) continue;
                                int nc2 = js[2] / opt_l[parts_idx[2]][o2];
                                int av2 = pm[parts_idx[2]] * opt_d[parts_idx[2]][o2];
                                if (nc2 > av2) continue;
                                long long bv2 = binomial(av2, nc2);
                                rs[2] = parts_list[pk[parts_idx[2]]];
                                ds[2] = opt_d[parts_idx[2]][o2];

                                long long cnt = count_gcd_pattern(rs, ds, 3, L);
                                total += cnt * bv0 * bv1 * bv2;
                            }
                        }
                    }
                }
            }
        }
    }

    /* p=4: composition (1,1,1,1) — only one permutation */
    if (depth >= 4) {
        for (int a = 0; a < depth; a++) {
            for (int b = a + 1; b < depth; b++) {
                for (int c = b + 1; c < depth; c++) {
                    for (int d_idx = c + 1; d_idx < depth; d_idx++) {
                        int pi[4] = {a, b, c, d_idx};
                        /* j = 1 for all four parts */
                        /* Need l | 1, so l must be 1, so d_val = r */
                        /* Check: l=1 means gcd(t,r) = r, i.e., r|t */
                        /* C(m*r, 1) = m*r */

                        int rs[4], ds[4];
                        long long bv = 1;
                        int valid = 1;
                        for (int i = 0; i < 4; i++) {
                            int r = parts_list[pk[pi[i]]];
                            int m = pm[pi[i]];
                            rs[i] = r;
                            ds[i] = r;  /* l=1, d=r */
                            bv *= m * r;  /* C(m*r, 1) */
                        }

                        long long cnt = count_gcd_pattern(rs, ds, 4, L);
                        total += cnt * bv;
                    }
                }
            }
        }
    }

    if (total % L != 0) {
        fprintf(stderr, "ERROR: c4 not integer: %lld/%lld\n", total, L);
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
        long long ck = compute_c4(pk, pm, depth_val);
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
            long long ck = compute_c4(pk, pm, 1);
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

    if (n < 4) {
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

    /* Max c_4 value: C(n, 4) */
    long long cnk = (long long)n * (n-1) * (n-2) * (n-3) / 24;
    max_t_val = (int)cnk + 1;
    if (max_t_val > 50000000) max_t_val = 50000000;

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
