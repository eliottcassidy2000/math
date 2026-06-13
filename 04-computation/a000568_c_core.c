/*
 * a000568_c_core.c — C implementation of the A000568 mod-p DP.
 *
 * Compile: gcc -O3 -shared -fPIC -o a000568_c_core.so a000568_c_core.c
 *
 * The DP is implemented using hash tables for the state space.
 * States are (total, S_d1, S_d2, ..., S_dk) tuples with integer weights mod p.
 *
 * Author: opus-2026-03-07-S46f
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* Hash table for DP states */
#define MAX_DIVS 40
#define MAX_STATES 500000

typedef struct {
    int total;
    int state[MAX_DIVS];
    int ndivs;
} StateKey;

typedef struct Entry {
    StateKey key;
    int64_t weight;  /* mod p */
    int occupied;
} Entry;

typedef struct {
    Entry *entries;
    int capacity;
    int size;
} HashMap;

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

static HashMap *hm_create(int capacity) {
    HashMap *hm = malloc(sizeof(HashMap));
    /* Round up to power of 2 for fast modulo */
    int cap = 1;
    while (cap < capacity * 2) cap <<= 1;  /* 2x load factor headroom */
    hm->capacity = cap;
    hm->size = 0;
    hm->entries = calloc(cap, sizeof(Entry));
    return hm;
}

static void hm_destroy(HashMap *hm) {
    free(hm->entries);
    free(hm);
}

static void hm_add(HashMap *hm, const StateKey *key, int64_t weight, int64_t p) {
    uint64_t mask = (uint64_t)(hm->capacity - 1);
    uint64_t h = hash_state(key) & mask;
    while (1) {
        Entry *e = &hm->entries[h];
        if (!e->occupied) {
            e->key = *key;
            e->weight = weight % p;
            if (e->weight < 0) e->weight += p;
            e->occupied = 1;
            hm->size++;
            return;
        }
        if (keys_equal(&e->key, key)) {
            e->weight = (e->weight + weight) % p;
            if (e->weight < 0) e->weight += p;
            return;
        }
        h = (h + 1) & mask;
    }
}

static int64_t mod_pow(int64_t base, int64_t exp, int64_t mod) {
    int64_t result = 1;
    base %= mod;
    if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1) {
            result = (__int128)result * base % mod;
        }
        base = (__int128)base * base % mod;
        exp >>= 1;
    }
    return result;
}

static int64_t mod_inv(int64_t x, int64_t p) {
    return mod_pow(x, p - 2, p);
}

static int64_t mod_mul(int64_t a, int64_t b, int64_t p) {
    return (__int128)a * b % p;
}

/* Euler totient */
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

/*
 * Compute A000568(n) mod p.
 * p must be prime and > n.
 * Returns the result.
 */
int64_t a000568_mod_p(int n, int64_t p) {
    if (n <= 1) return 1;

    /* Precompute phi */
    int phi[n + 1];
    for (int i = 1; i <= n; i++) phi[i] = euler_totient(i);

    /* Odd parts */
    int odd_parts[(n + 1) / 2 + 1];
    int num_odd = 0;
    for (int k = 1; k <= n; k += 2) odd_parts[num_odd++] = k;

    /* Divisors of each odd part */
    int divisors[n + 1][32];
    int ndivs[n + 1];
    memset(ndivs, 0, sizeof(ndivs));
    for (int i = 0; i < num_odd; i++) {
        int k = odd_parts[i];
        for (int d = 1; d <= k; d++) {
            if (k % d == 0 && d % 2 == 1) {
                divisors[k][ndivs[k]++] = d;
            }
        }
    }

    /* Max part for each divisor */
    int max_part_for_div[n + 1];
    memset(max_part_for_div, 0, sizeof(max_part_for_div));
    for (int i = 0; i < num_odd; i++) {
        int k = odd_parts[i];
        for (int j = 0; j < ndivs[k]; j++) {
            max_part_for_div[divisors[k][j]] = k;
        }
    }

    /* Active divisors at each stage */
    /* Start: all divisors that divide some part >= odd_parts[0] */
    int active[MAX_DIVS];
    int num_active = 0;
    for (int d = 1; d <= n; d += 2) {
        if (max_part_for_div[d] >= odd_parts[0]) {
            active[num_active++] = d;
        }
    }

    int div_index[n + 1];
    memset(div_index, -1, sizeof(div_index));
    for (int i = 0; i < num_active; i++) {
        div_index[active[i]] = i;
    }

    /* Init DP */
    HashMap *dp = hm_create(n * n * n + 1024);
    StateKey init_key;
    init_key.total = 0;
    init_key.ndivs = num_active;
    memset(init_key.state, 0, sizeof(init_key.state));
    hm_add(dp, &init_key, 1, p);

    /* Precompute factorials mod p */
    int64_t fact[n + 1];
    fact[0] = 1;
    for (int i = 1; i <= n; i++) fact[i] = mod_mul(fact[i - 1], i, p);

    for (int ki = 0; ki < num_odd; ki++) {
        int k = odd_parts[ki];
        if (k > n) break;

        int div_idx_k[32], phi_k[32];
        int num_dk = 0;
        for (int j = 0; j < ndivs[k]; j++) {
            int d = divisors[k][j];
            if (div_index[d] >= 0) {
                div_idx_k[num_dk] = div_index[d];
                phi_k[num_dk] = phi[d];
                num_dk++;
            }
        }

        int max_m = n / k;

        /* Precompute k^m mod p and (m!)^{-1} mod p */
        int64_t k_pow[max_m + 1], fact_inv[max_m + 1];
        k_pow[0] = 1;
        for (int m = 1; m <= max_m; m++) k_pow[m] = mod_mul(k_pow[m - 1], k, p);
        for (int m = 0; m <= max_m; m++) fact_inv[m] = mod_inv(fact[m], p);

        HashMap *new_dp = hm_create(n * n * n + 1024);

        for (int idx = 0; idx < dp->capacity; idx++) {
            Entry *e = &dp->entries[idx];
            if (!e->occupied) continue;

            int total = e->key.total;
            int64_t weight = e->weight;
            int max_m_here = (n - total) / k;

            for (int m = 0; m <= max_m_here; m++) {
                int new_total = total + m * k;

                /* Compute cross */
                int64_t cross = 0;
                for (int i = 0; i < num_dk; i++) {
                    cross += (int64_t)phi_k[i] * e->key.state[div_idx_k[i]];
                }

                int64_t delta_t = (int64_t)m * cross +
                    (int64_t)m * (m - 1) * k / 2 +
                    (int64_t)m * (k - 1) / 2;

                int64_t pow2 = mod_pow(2, delta_t, p);
                int64_t inv_km = mod_inv(k_pow[m], p);
                int64_t factor = mod_mul(mod_mul(pow2, inv_km, p), fact_inv[m], p);

                StateKey new_key;
                new_key.total = new_total;
                new_key.ndivs = num_active;
                memcpy(new_key.state, e->key.state, num_active * sizeof(int));
                for (int i = 0; i < num_dk; i++) {
                    new_key.state[div_idx_k[i]] += m;
                }

                int64_t val = mod_mul(weight, factor, p);
                hm_add(new_dp, &new_key, val, p);
            }
        }

        hm_destroy(dp);
        dp = new_dp;

        /* Project state */
        if (ki + 1 < num_odd) {
            int next_k = odd_parts[ki + 1];
            int new_active[MAX_DIVS];
            int new_num_active = 0;
            for (int d = 1; d <= n; d += 2) {
                if (max_part_for_div[d] >= next_k) {
                    new_active[new_num_active++] = d;
                }
            }

            if (new_num_active < num_active) {
                int new_div_index[n + 1];
                memset(new_div_index, -1, sizeof(new_div_index));
                for (int i = 0; i < new_num_active; i++) {
                    new_div_index[new_active[i]] = i;
                }

                int keep_indices[MAX_DIVS];
                int nkeep = 0;
                for (int i = 0; i < new_num_active; i++) {
                    if (div_index[new_active[i]] >= 0) {
                        keep_indices[nkeep++] = div_index[new_active[i]];
                    }
                }

                HashMap *projected = hm_create(n * n * n + 1024);
                for (int idx = 0; idx < dp->capacity; idx++) {
                    Entry *e = &dp->entries[idx];
                    if (!e->occupied) continue;

                    StateKey new_key;
                    new_key.total = e->key.total;
                    new_key.ndivs = new_num_active;
                    for (int i = 0; i < nkeep; i++) {
                        new_key.state[i] = e->key.state[keep_indices[i]];
                    }

                    hm_add(projected, &new_key, e->weight, p);
                }

                hm_destroy(dp);
                dp = projected;

                num_active = new_num_active;
                memcpy(active, new_active, new_num_active * sizeof(int));
                memset(div_index, -1, sizeof(div_index));
                for (int i = 0; i < num_active; i++) {
                    div_index[active[i]] = i;
                }
            }
        }
    }

    /* Sum states with total = n */
    int64_t result = 0;
    for (int idx = 0; idx < dp->capacity; idx++) {
        Entry *e = &dp->entries[idx];
        if (e->occupied && e->key.total == n) {
            result = (result + e->weight) % p;
        }
    }

    hm_destroy(dp);
    return result;
}

/* Test harness */
int main(int argc, char *argv[]) {
    int n = 12;
    int64_t p = 1000000007;

    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) p = atoll(argv[2]);

    int64_t result = a000568_mod_p(n, p);
    printf("a(%d) mod %lld = %lld\n", n, (long long)p, (long long)result);

    return 0;
}
