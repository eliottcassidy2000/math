/*
 * a000568_crt_v2.c — opus-2026-04-04-S23
 *
 * Fast A000568 via CRT. Fixed version: no nested functions.
 * Each prime computed independently with native 64-bit arithmetic.
 *
 * Compile: gcc -O3 -o a000568_crt_v2 a000568_crt_v2.c -lgmp -lm
 * Usage:   ./a000568_crt_v2 <n>
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
        base = mulmod(base, base, mod); exp >>= 1;
    }
    return r;
}

static int gcd_int(int a, int b) {
    while (b) { int t = b; b = a % b; a = t; } return a;
}

/* ---- Partition enumeration state (global for recursion) ---- */
static int odd_parts[500], num_odd;
static int gcd_tab[500][500];
static int pl_idx[100], pl_mul[100], depth;
static u64 acc, mod_p, nfact_p;
static u64 *pow2_p;
static int max_e;

static void enum_parts(int rem, int mpi) {
    if (rem == 0) {
        i64 edges = 0;
        u64 z = 1;
        for (int i = 0; i < depth; i++) {
            int k = odd_parts[pl_idx[i]], m = pl_mul[i];
            edges += (i64)m * (k - 1) / 2;
            edges += (i64)m * (m - 1) / 2 * k;
            z = mulmod(z, powmod(k, m, mod_p), mod_p);
            for (int f = 1; f <= m; f++) z = mulmod(z, f, mod_p);
            for (int j = i + 1; j < depth; j++)
                edges += (i64)m * pl_mul[j] * gcd_tab[pl_idx[i]][pl_idx[j]];
        }
        if (edges > max_e) edges = max_e;
        if (z == 0) return;
        u64 zi = powmod(z, mod_p - 2, mod_p);
        u64 c = mulmod(mulmod(nfact_p, zi, mod_p), pow2_p[edges], mod_p);
        acc = (acc + c) % mod_p;
        return;
    }
    for (int pi = mpi; pi >= 0; pi--) {
        int k = odd_parts[pi];
        if (k > rem) continue;
        for (int m = 1; m * k <= rem; m++) {
            pl_idx[depth] = pi; pl_mul[depth] = m; depth++;
            enum_parts(rem - m * k, pi - 1);
            depth--;
        }
    }
}

static u64 burnside_mod(int n, u64 p) {
    if (n <= 2) return 1 % p;
    num_odd = 0;
    for (int k = 1; k <= n; k += 2) odd_parts[num_odd++] = k;
    for (int i = 0; i < num_odd; i++)
        for (int j = 0; j < num_odd; j++)
            gcd_tab[i][j] = gcd_int(odd_parts[i], odd_parts[j]);
    u64 nf = 1;
    for (int k = 2; k <= n; k++) nf = mulmod(nf, k, p);
    int me = n * (n - 1) / 2;
    u64 *p2 = malloc((me + 1) * sizeof(u64));
    p2[0] = 1;
    for (int i = 1; i <= me; i++) p2[i] = mulmod(p2[i-1], 2, p);

    mod_p = p; nfact_p = nf; pow2_p = p2; max_e = me;
    depth = 0; acc = 0;
    enum_parts(n, num_odd - 1);
    free(p2);
    return mulmod(acc, powmod(nf, p - 2, p), p);
}

/* ---- Prime generation ---- */
static int is_prime64(u64 n) {
    if (n < 2) return 0; if (n < 4) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    for (u64 i = 5; i * i <= n; i += 6)
        if (n % i == 0 || n % (i+2) == 0) return 0;
    return 1;
}
static u64 next_prime64(u64 n) {
    n |= 1; while (!is_prime64(n)) n += 2; return n;
}

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <n>\n", argv[0]); return 1; }
    int n = atoi(argv[1]);

    double el2 = n * (n - 1) / 2.0;
    for (int k = 2; k <= n; k++) el2 -= log2(k);
    int ebits = (int)el2 + 20;
    int np = ebits / 62 + 2;
    if (np < 2) np = 2;

    u64 *primes = malloc(np * sizeof(u64));
    u64 *resids = malloc(np * sizeof(u64));
    u64 p = (1ULL << 62) + 1;
    for (int i = 0; i < np; i++) { p = next_prime64(p); primes[i] = p; p += 2; }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (int i = 0; i < np; i++)
        resids[i] = burnside_mod(n, primes[i]);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double sec = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

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
    double total = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    printf("n=%d, %d primes, mod=%.3fs, total=%.3fs\n", n, np, sec, total);
    printf("a(%d) = ", n);
    mpz_out_str(stdout, 10, result);
    printf("\n");

    mpz_clear(result); mpz_clear(modulus); mpz_clear(temp);
    free(primes); free(resids);
    return 0;
}
