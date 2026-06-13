/*
 * a000568_c_batch.c — Batched CRT: one DP pass, vectorized over primes.
 *
 * Each DP weight is an array of int64_t, one per prime.
 * All arithmetic is modular and uses machine words.
 *
 * Compile: gcc -O3 -o a000568_c_batch a000568_c_batch.c -lm
 *
 * Usage: ./a000568_c_batch N
 * Outputs: a(N) as a decimal number.
 *
 * Author: opus-2026-03-07-S46f
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define MAX_DIVS 40
#define MAX_PRIMES 200

/* ---- Prime generation ---- */

static int is_prime(int64_t n) {
    if (n < 2) return 0;
    if (n < 4) return 1;
    if (n % 2 == 0) return 0;
    for (int64_t i = 3; i * i <= n; i += 2)
        if (n % i == 0) return 0;
    return 1;
}

static int gen_primes_above(int64_t start, int count, int64_t *out) {
    int found = 0;
    int64_t n = start + 1;
    if (n % 2 == 0) n++;
    while (found < count) {
        if (is_prime(n)) out[found++] = n;
        n += 2;
    }
    return found;
}

/* ---- Modular arithmetic ---- */

static int64_t mod_pow(int64_t base, int64_t exp, int64_t mod) {
    int64_t result = 1;
    base %= mod;
    if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1)
            result = ((__int128)result * base) % mod;
        base = ((__int128)base * base) % mod;
        exp >>= 1;
    }
    return result;
}

/* ---- Euler totient ---- */

static int euler_totient(int n) {
    int result = n, temp = n;
    for (int p = 2; p * p <= temp; p++) {
        if (temp % p == 0) {
            while (temp % p == 0) temp /= p;
            result -= result / p;
        }
    }
    if (temp > 1) result -= result / temp;
    return result;
}

/* ---- Hash table with vector weights ---- */

typedef struct {
    int total;
    int state[MAX_DIVS];
    int ndivs;
} StateKey;

typedef struct {
    StateKey key;
    int64_t *weight;  /* array of num_primes int64_t values */
    int occupied;
} VEntry;

typedef struct {
    VEntry *entries;
    int capacity;
    int size;
    int num_primes;
    int64_t *primes;
} VHashMap;

static uint64_t hash_state(const StateKey *k) {
    uint64_t h = (uint64_t)k->total * 2654435761ULL;
    for (int i = 0; i < k->ndivs; i++) {
        h ^= (uint64_t)k->state[i] * (2654435761ULL + 2 * i);
        h = (h << 13) | (h >> 51);
    }
    return h;
}

static int keys_equal(const StateKey *a, const StateKey *b) {
    if (a->total != b->total || a->ndivs != b->ndivs) return 0;
    return memcmp(a->state, b->state, a->ndivs * sizeof(int)) == 0;
}

static VHashMap *vhm_create(int capacity, int num_primes, int64_t *primes) {
    VHashMap *hm = malloc(sizeof(VHashMap));
    int cap = 1;
    while (cap < capacity * 2) cap <<= 1;
    hm->capacity = cap;
    hm->size = 0;
    hm->num_primes = num_primes;
    hm->primes = primes;
    hm->entries = calloc(cap, sizeof(VEntry));
    return hm;
}

static void vhm_destroy(VHashMap *hm) {
    for (int i = 0; i < hm->capacity; i++) {
        if (hm->entries[i].occupied)
            free(hm->entries[i].weight);
    }
    free(hm->entries);
    free(hm);
}

static void vhm_add(VHashMap *hm, const StateKey *key, const int64_t *weight) {
    uint64_t mask = (uint64_t)(hm->capacity - 1);
    uint64_t h = hash_state(key) & mask;
    int NP = hm->num_primes;
    int64_t *P = hm->primes;

    while (1) {
        VEntry *e = &hm->entries[h];
        if (!e->occupied) {
            e->key = *key;
            e->weight = malloc(NP * sizeof(int64_t));
            for (int j = 0; j < NP; j++) {
                e->weight[j] = weight[j] % P[j];
                if (e->weight[j] < 0) e->weight[j] += P[j];
            }
            e->occupied = 1;
            hm->size++;
            return;
        }
        if (keys_equal(&e->key, key)) {
            for (int j = 0; j < NP; j++) {
                e->weight[j] = (e->weight[j] + weight[j]) % P[j];
                if (e->weight[j] < 0) e->weight[j] += P[j];
            }
            return;
        }
        h = (h + 1) & mask;
    }
}

/* ---- Main computation ---- */

/* Big number print (for CRT result) using Python-style arbitrary precision */
/* We'll use a simpler approach: output residues and let Python do CRT */

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s N\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    if (n <= 1) {
        printf("1\n");
        return 0;
    }

    /* Estimate bits needed */
    double est_bits = n * (n - 1) / 2.0 - n * log2(fmax(n, 2)) + n + 100;
    if (est_bits < 100) est_bits = 100;

    int prime_bits = 30;
    int num_primes = (int)(est_bits / prime_bits) + 5;
    if (num_primes > MAX_PRIMES) num_primes = MAX_PRIMES;

    int64_t primes[MAX_PRIMES];
    int64_t start = n > (1 << 30) ? n : (1 << 30);
    gen_primes_above(start, num_primes, primes);

    fprintf(stderr, "n=%d: %d primes (%d-bit), est %.0f bits\n",
            n, num_primes, prime_bits, est_bits);

    /* Precompute */
    int phi[n + 1];
    for (int i = 1; i <= n; i++) phi[i] = euler_totient(i);

    int odd_parts[(n + 1) / 2 + 1];
    int num_odd = 0;
    for (int k = 1; k <= n; k += 2) odd_parts[num_odd++] = k;

    int divisors[n + 1][32];
    int ndivs_arr[n + 1];
    memset(ndivs_arr, 0, sizeof(ndivs_arr));
    for (int i = 0; i < num_odd; i++) {
        int k = odd_parts[i];
        for (int d = 1; d <= k; d++)
            if (k % d == 0 && d % 2 == 1)
                divisors[k][ndivs_arr[k]++] = d;
    }

    int max_part_for_div[n + 1];
    memset(max_part_for_div, 0, sizeof(max_part_for_div));
    for (int i = 0; i < num_odd; i++) {
        int k = odd_parts[i];
        for (int j = 0; j < ndivs_arr[k]; j++)
            max_part_for_div[divisors[k][j]] = k;
    }

    int active[MAX_DIVS];
    int num_active = 0;
    for (int d = 1; d <= n; d += 2)
        if (max_part_for_div[d] >= odd_parts[0])
            active[num_active++] = d;

    int div_index[n + 1];
    memset(div_index, -1, sizeof(div_index));
    for (int i = 0; i < num_active; i++)
        div_index[active[i]] = i;

    /* Factorials mod each prime */
    int64_t fact_mod[n + 1][MAX_PRIMES];
    for (int j = 0; j < num_primes; j++) fact_mod[0][j] = 1;
    for (int i = 1; i <= n; i++)
        for (int j = 0; j < num_primes; j++)
            fact_mod[i][j] = ((__int128)fact_mod[i-1][j] * i) % primes[j];

    /* Init DP */
    int ht_cap = n * n * n + 1024;
    VHashMap *dp = vhm_create(ht_cap, num_primes, primes);

    StateKey init_key;
    init_key.total = 0;
    init_key.ndivs = num_active;
    memset(init_key.state, 0, sizeof(init_key.state));

    int64_t ones[MAX_PRIMES];
    for (int j = 0; j < num_primes; j++) ones[j] = 1;
    vhm_add(dp, &init_key, ones);

    /* Precompute 2^{2^i} mod each prime for fast pow2 */
    int64_t two_pow2i[64][MAX_PRIMES];
    for (int j = 0; j < num_primes; j++) two_pow2i[0][j] = 2 % primes[j];
    for (int i = 1; i < 64; i++)
        for (int j = 0; j < num_primes; j++)
            two_pow2i[i][j] = ((__int128)two_pow2i[i-1][j] * two_pow2i[i-1][j]) % primes[j];

    /* Fast vectorized pow(2, exp, primes[j]) */
    int64_t tmp_pow2[MAX_PRIMES];

    /* Temp arrays */
    int64_t weight_tmp[MAX_PRIMES];

    for (int ki = 0; ki < num_odd; ki++) {
        int k = odd_parts[ki];
        if (k > n) break;

        int div_idx_k[32], phi_k[32];
        int num_dk = 0;
        for (int j = 0; j < ndivs_arr[k]; j++) {
            int d = divisors[k][j];
            if (div_index[d] >= 0) {
                div_idx_k[num_dk] = div_index[d];
                phi_k[num_dk] = phi[d];
                num_dk++;
            }
        }

        int max_m = n / k;

        /* Precompute k^{-1} mod p, then (k^m)^{-1} mod p */
        int64_t k_inv[MAX_PRIMES];
        for (int j = 0; j < num_primes; j++)
            k_inv[j] = mod_pow(k, primes[j] - 2, primes[j]);

        int64_t k_pow_inv[max_m + 1][num_primes];
        for (int j = 0; j < num_primes; j++) k_pow_inv[0][j] = 1;
        for (int m = 1; m <= max_m; m++)
            for (int j = 0; j < num_primes; j++)
                k_pow_inv[m][j] = ((__int128)k_pow_inv[m-1][j] * k_inv[j]) % primes[j];

        /* (m!)^{-1} mod p */
        int64_t fact_inv[max_m + 1][num_primes];
        for (int m = 0; m <= max_m; m++)
            for (int j = 0; j < num_primes; j++)
                fact_inv[m][j] = mod_pow(fact_mod[m][j], primes[j] - 2, primes[j]);

        /* Precompute self_t powers and combined coeff */
        int64_t coeff[max_m + 1][num_primes];
        for (int m = 0; m <= max_m; m++) {
            int64_t st = (int64_t)m * (m - 1) * k / 2 + (int64_t)m * (k - 1) / 2;
            /* 2^st mod each prime */
            for (int j = 0; j < num_primes; j++) {
                int64_t r = 1;
                int64_t e = st;
                int bit = 0;
                while (e > 0) {
                    if (e & 1) r = ((__int128)r * two_pow2i[bit][j]) % primes[j];
                    e >>= 1;
                    bit++;
                }
                coeff[m][j] = ((__int128)((__int128)r * k_pow_inv[m][j] % primes[j]) * fact_inv[m][j]) % primes[j];
            }
        }

        VHashMap *new_dp = vhm_create(ht_cap, num_primes, primes);

        for (int idx = 0; idx < dp->capacity; idx++) {
            VEntry *e = &dp->entries[idx];
            if (!e->occupied) continue;

            int total = e->key.total;
            int max_m_here = (n - total) / k;

            int64_t cross = 0;
            for (int i = 0; i < num_dk; i++)
                cross += (int64_t)phi_k[i] * e->key.state[div_idx_k[i]];

            /* 2^cross mod each prime */
            int64_t two_c[MAX_PRIMES];
            {
                int64_t ec = cross;
                for (int j = 0; j < num_primes; j++) {
                    int64_t r = 1;
                    int64_t ee = ec;
                    int bit = 0;
                    while (ee > 0) {
                        if (ee & 1) r = ((__int128)r * two_pow2i[bit][j]) % primes[j];
                        ee >>= 1;
                        bit++;
                    }
                    two_c[j] = r;
                }
            }

            /* (2^cross)^m via repeated multiply */
            int64_t two_c_pow[MAX_PRIMES];
            for (int j = 0; j < num_primes; j++) two_c_pow[j] = 1;

            for (int m = 0; m <= max_m_here; m++) {
                int new_total = total + m * k;

                /* factor = two_c_pow * coeff[m] */
                for (int j = 0; j < num_primes; j++)
                    weight_tmp[j] = ((__int128)((__int128)e->weight[j] * two_c_pow[j] % primes[j]) * coeff[m][j]) % primes[j];

                StateKey new_key;
                new_key.total = new_total;
                new_key.ndivs = num_active;
                memcpy(new_key.state, e->key.state, num_active * sizeof(int));
                for (int i = 0; i < num_dk; i++)
                    new_key.state[div_idx_k[i]] += m;

                vhm_add(new_dp, &new_key, weight_tmp);

                /* Update for next m */
                for (int j = 0; j < num_primes; j++)
                    two_c_pow[j] = ((__int128)two_c_pow[j] * two_c[j]) % primes[j];
            }
        }

        vhm_destroy(dp);
        dp = new_dp;

        /* Project state */
        if (ki + 1 < num_odd) {
            int next_k = odd_parts[ki + 1];
            int new_active[MAX_DIVS];
            int new_num_active = 0;
            for (int d = 1; d <= n; d += 2)
                if (max_part_for_div[d] >= next_k)
                    new_active[new_num_active++] = d;

            if (new_num_active < num_active) {
                int new_div_index[n + 1];
                memset(new_div_index, -1, sizeof(new_div_index));
                for (int i = 0; i < new_num_active; i++)
                    new_div_index[new_active[i]] = i;

                int keep_indices[MAX_DIVS];
                int nkeep = 0;
                for (int i = 0; i < new_num_active; i++)
                    if (div_index[new_active[i]] >= 0)
                        keep_indices[nkeep++] = div_index[new_active[i]];

                VHashMap *projected = vhm_create(ht_cap, num_primes, primes);
                for (int idx = 0; idx < dp->capacity; idx++) {
                    VEntry *e = &dp->entries[idx];
                    if (!e->occupied) continue;

                    StateKey new_key;
                    new_key.total = e->key.total;
                    new_key.ndivs = new_num_active;
                    memset(new_key.state, 0, sizeof(new_key.state));
                    for (int i = 0; i < nkeep; i++)
                        new_key.state[i] = e->key.state[keep_indices[i]];

                    vhm_add(projected, &new_key, e->weight);
                }

                vhm_destroy(dp);
                dp = projected;
                num_active = new_num_active;
                memcpy(active, new_active, new_num_active * sizeof(int));
                memset(div_index, -1, sizeof(div_index));
                for (int i = 0; i < num_active; i++)
                    div_index[active[i]] = i;
            }
        }

        fprintf(stderr, "  k=%3d: %d states\n", k, dp->size);
    }

    /* Sum states with total = n — output residues for Python CRT */
    int64_t result[MAX_PRIMES];
    memset(result, 0, sizeof(result));
    for (int idx = 0; idx < dp->capacity; idx++) {
        VEntry *e = &dp->entries[idx];
        if (e->occupied && e->key.total == n) {
            for (int j = 0; j < num_primes; j++)
                result[j] = (result[j] + e->weight[j]) % primes[j];
        }
    }

    /* Output: primes and residues for Python CRT reconstruction */
    printf("%d\n", num_primes);
    for (int j = 0; j < num_primes; j++)
        printf("%lld %lld\n", (long long)primes[j], (long long)result[j]);

    vhm_destroy(dp);
    return 0;
}
