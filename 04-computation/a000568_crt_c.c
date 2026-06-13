/*
 * a000568_crt_c.c — opus-2026-04-04-S22
 * Fast A000568 via CRT with native 64-bit modular arithmetic.
 * Compile: gcc -O3 -o a000568_crt_c a000568_crt_c.c -lgmp -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>

typedef unsigned long long u64;
typedef long long i64;

static inline u64 mulmod(u64 a, u64 b, u64 m) {
    return (unsigned __int128)a * b % m;
}
static u64 powmod(u64 base, u64 exp, u64 mod) {
    u64 r = 1; base %= mod;
    while (exp > 0) {
        if (exp & 1) r = mulmod(r, base, mod);
        base = mulmod(base, base, mod);
        exp >>= 1;
    }
    return r;
}
static int gcd_int(int a, int b) {
    while (b) { int t = b; b = a % b; a = t; } return a;
}

/* Global state for enumeration */
static int odd_parts[500], num_odd;
static int gcd_tab_g[500][500];
static int placed_idx[100], placed_mult[100], depth;
static u64 total_acc, p_mod, nfact_mod_g;
static u64 *pow2_tab;
static int max_edges_g;

static void enum_partitions(int remaining, int max_pi) {
    if (remaining == 0) {
        i64 edges = 0;
        u64 z_mod = 1;
        for (int i = 0; i < depth; i++) {
            int ki = odd_parts[placed_idx[i]], mi = placed_mult[i];
            edges += (i64)mi * (ki - 1) / 2;
            edges += (i64)mi * (mi - 1) / 2 * ki;
            z_mod = mulmod(z_mod, powmod(ki, mi, p_mod), p_mod);
            for (int f = 1; f <= mi; f++)
                z_mod = mulmod(z_mod, f, p_mod);
            for (int j = i + 1; j < depth; j++) {
                int mj = placed_mult[j];
                edges += (i64)mi * mj * gcd_tab_g[placed_idx[i]][placed_idx[j]];
            }
        }
        if (edges > max_edges_g) edges = max_edges_g;
        if (z_mod == 0) return;
        u64 z_inv = powmod(z_mod, p_mod - 2, p_mod);
        u64 c = mulmod(mulmod(nfact_mod_g, z_inv, p_mod), pow2_tab[edges], p_mod);
        total_acc = (total_acc + c) % p_mod;
        return;
    }
    for (int pi = max_pi; pi >= 0; pi--) {
        int k = odd_parts[pi];
        if (k > remaining) continue;
        for (int m = 1; m * k <= remaining; m++) {
            placed_idx[depth] = pi;
            placed_mult[depth] = m;
            depth++;
            enum_partitions(remaining - m * k, pi - 1);
            depth--;
        }
    }
}

static u64 burnside_mod_p(int n, u64 p) {
    if (n <= 2) return 1 % p;
    num_odd = 0;
    for (int k = 1; k <= n; k += 2) odd_parts[num_odd++] = k;
    for (int i = 0; i < num_odd; i++)
        for (int j = 0; j < num_odd; j++)
            gcd_tab_g[i][j] = gcd_int(odd_parts[i], odd_parts[j]);

    u64 nf = 1;
    for (int k = 2; k <= n; k++) nf = mulmod(nf, k, p);
    u64 nf_inv = powmod(nf, p - 2, p);

    int me = n * (n - 1) / 2;
    pow2_tab = (u64*)malloc((me + 1) * sizeof(u64));
    pow2_tab[0] = 1;
    for (int i = 1; i <= me; i++) pow2_tab[i] = mulmod(pow2_tab[i-1], 2, p);

    p_mod = p; nfact_mod_g = nf; max_edges_g = me;
    depth = 0; total_acc = 0;
    enum_partitions(n, num_odd - 1);

    free(pow2_tab);
    return mulmod(total_acc, nf_inv, p);
}

static int is_prime_u64(u64 n) {
    if (n < 2) return 0;
    if (n < 4) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    for (u64 i = 5; i * i <= n; i += 6)
        if (n % i == 0 || n % (i+2) == 0) return 0;
    return 1;
}
static u64 next_prime_u64(u64 n) {
    n |= 1;
    while (!is_prime_u64(n)) n += 2;
    return n;
}

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <n>\n", argv[0]); return 1; }
    int n = atoi(argv[1]);

    double el2 = n * (n-1) / 2.0;
    for (int k = 2; k <= n; k++) el2 -= log2(k);
    int ebits = (int)el2 + 20;
    int np = ebits / 62 + 2;

    printf("n=%d: ~%d bits, %d primes\n", n, ebits, np);

    u64 *primes = malloc(np * sizeof(u64));
    u64 *resids = malloc(np * sizeof(u64));
    u64 p = (1ULL << 62) + 1;
    for (int i = 0; i < np; i++) { p = next_prime_u64(p); primes[i] = p; p += 2; }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int i = 0; i < np; i++) resids[i] = burnside_mod_p(n, primes[i]);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double sec = (t1.tv_sec-t0.tv_sec) + (t1.tv_nsec-t0.tv_nsec)/1e9;
    printf("Mod computation: %.3fs\n", sec);

    /* CRT with GMP */
    mpz_t result, modulus, temp;
    mpz_init_set_ui(result, resids[0]);
    mpz_init_set_ui(modulus, primes[0]);
    mpz_init(temp);
    for (int i = 1; i < np; i++) {
        u64 r = resids[i], pi = primes[i];
        u64 rm = mpz_fdiv_ui(result, pi);
        u64 d = (r + pi - rm) % pi;
        u64 mi = powmod(mpz_fdiv_ui(modulus, pi), pi - 2, pi);
        u64 tv = mulmod(d, mi, pi);
        mpz_mul_ui(temp, modulus, tv);
        mpz_add(result, result, temp);
        mpz_mul_ui(modulus, modulus, pi);
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    sec = (t1.tv_sec-t0.tv_sec) + (t1.tv_nsec-t0.tv_nsec)/1e9;

    printf("a(%d) = ", n);
    mpz_out_str(stdout, 10, result);
    printf("\nTime: %.3fs, Partitions per prime: ~%d\n", sec, num_odd > 0 ? num_odd : 0);

    mpz_clear(result); mpz_clear(modulus); mpz_clear(temp);
    free(primes); free(resids);
    return 0;
}
