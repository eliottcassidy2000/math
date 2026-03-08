/*
 * a000568_c_enum.c — A000568 via enumeration with CRT.
 *
 * Enumerates all partitions of n into odd parts, computes contribution
 * 2^{t(λ)} / z_λ mod p for each, sums them.
 *
 * Compile: gcc -O3 -o a000568_c_enum a000568_c_enum.c -lm
 * Usage: ./a000568_c_enum <n> <prime>
 *
 * Author: opus-2026-03-07-S46f
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int64_t MOD_P;

static int64_t mod_pow(int64_t base, int64_t exp, int64_t mod) {
    int64_t result = 1;
    base %= mod;
    if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1)
            result = (__int128)result * base % mod;
        base = (__int128)base * base % mod;
        exp >>= 1;
    }
    return result;
}

static int64_t mod_inv(int64_t x, int64_t p) {
    return mod_pow(x, p - 2, p);
}

/* GCD table */
static int gcd_tab[502][502];
static void init_gcd(int n) {
    for (int i = 1; i <= n; i += 2)
        for (int j = i; j <= n; j += 2) {
            int a = i, b = j;
            while (b) { int t = b; b = a % b; a = t; }
            gcd_tab[i][j] = gcd_tab[j][i] = a;
        }
}

/* Precomputed values */
static int n_val;
static int odd_parts[256];
static int num_odd;
static int64_t fact_mod[502];      /* i! mod p */
static int64_t fact_inv_mod[502];  /* (i!)^{-1} mod p */
static int64_t result_sum;

/* Stack for partition enumeration */
typedef struct {
    int k;
    int m;
} PartEntry;
static PartEntry stack[256];
static int stack_size = 0;

static void enumerate(int remaining, int max_part_idx, int64_t t_val, int64_t z_inv) {
    if (remaining == 0) {
        /* Contribution: 2^{t_val} * z_inv mod p */
        int64_t pow2 = mod_pow(2, t_val, MOD_P);
        int64_t contrib = (__int128)pow2 * z_inv % MOD_P;
        result_sum = (result_sum + contrib) % MOD_P;
        return;
    }

    for (int pi = max_part_idx; pi >= 0; pi--) {
        int k = odd_parts[pi];
        if (k > remaining) continue;
        int max_m = remaining / k;

        /* For each multiplicity m of part k */
        int64_t k_inv_pow = 1;  /* (k^{-1})^m mod p */
        int64_t k_inv = mod_inv(k, MOD_P);

        for (int m = 1; m <= max_m; m++) {
            /* Incremental update */
            k_inv_pow = (__int128)k_inv_pow * k_inv % MOD_P;

            /* Self contribution to t */
            int64_t dt_self = (int64_t)m * (m - 1) * k / 2 + (int64_t)m * (k - 1) / 2;

            /* Cross contribution to t */
            int64_t dt_cross = 0;
            for (int s = 0; s < stack_size; s++) {
                dt_cross += (int64_t)m * stack[s].m * gcd_tab[k][stack[s].k];
            }

            int64_t new_t = t_val + dt_self + dt_cross;

            /* z_inv update: z_inv *= k^{-m} * (m!)^{-1} */
            int64_t new_z_inv = (__int128)z_inv * k_inv_pow % MOD_P;
            new_z_inv = (__int128)new_z_inv * fact_inv_mod[m] % MOD_P;

            stack[stack_size].k = k;
            stack[stack_size].m = m;
            stack_size++;

            enumerate(remaining - m * k, pi - 1, new_t, new_z_inv);

            stack_size--;
        }
    }
}

int64_t a000568_mod_p_enum(int n, int64_t p) {
    MOD_P = p;
    n_val = n;
    if (n <= 1) return 1 % p;

    init_gcd(n);

    num_odd = 0;
    for (int k = 1; k <= n; k += 2)
        odd_parts[num_odd++] = k;

    /* Precompute factorials and their inverses */
    fact_mod[0] = 1;
    for (int i = 1; i <= n; i++)
        fact_mod[i] = (__int128)fact_mod[i-1] * i % p;
    fact_inv_mod[n] = mod_inv(fact_mod[n], p);
    for (int i = n - 1; i >= 0; i--)
        fact_inv_mod[i] = (__int128)fact_inv_mod[i+1] * (i+1) % p;

    result_sum = 0;
    stack_size = 0;
    enumerate(n, num_odd - 1, 0, 1);

    return result_sum;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <n> <prime>\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    int64_t p = atoll(argv[2]);

    int64_t result = a000568_mod_p_enum(n, p);
    printf("%lld\n", (long long)result);

    return 0;
}
