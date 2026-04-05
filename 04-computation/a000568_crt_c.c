/*
 * a000568_crt_c.c — opus-2026-04-04-S22
 *
 * Fast A000568 via CRT: compute a(n) mod p for many primes using
 * native 64-bit modular arithmetic, reconstruct via CRT with GMP.
 *
 * BINARY DECISION PERSPECTIVE: Each tournament = C(n,2) binary choices.
 * Burnside compresses into partition sum. CRT compresses bignum arithmetic.
 * Each prime runs independently → perfect parallelism.
 *
 * Compile: gcc -O3 -o a000568_crt_c a000568_crt_c.c -lgmp -lm
 * Usage: ./a000568_crt_c <n>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>

typedef unsigned long long u64;
typedef long long i64;

/* Modular arithmetic */
static inline u64 mulmod(u64 a, u64 b, u64 m) {
    return (unsigned __int128)a * b % m;
}

static u64 powmod(u64 base, u64 exp, u64 mod) {
    u64 result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = mulmod(result, base, mod);
        base = mulmod(base, base, mod);
        exp >>= 1;
    }
    return result;
}

/* GCD */
static int gcd_int(int a, int b) {
    while (b) { int t = b; b = a % b; a = t; }
    return a;
}

/* Odd-part partitions of n: enumerate and compute Burnside sum mod p */
static int odd_parts[500];   /* available odd parts */
static int num_odd;
static int gcd_tab[500][500]; /* precomputed GCDs */
static int n_global;

/* The recursive partition enumerator */
static u64 burnside_mod_p(int n, u64 p) {
    /* Precompute */
    num_odd = 0;
    for (int k = 1; k <= n; k += 2)
        odd_parts[num_odd++] = k;

    for (int i = 0; i < num_odd; i++)
        for (int j = 0; j < num_odd; j++)
            gcd_tab[i][j] = gcd_int(odd_parts[i], odd_parts[j]);

    /* n! mod p */
    u64 nfact = 1;
    for (int k = 2; k <= n; k++)
        nfact = mulmod(nfact, k, p);
    u64 nfact_inv = powmod(nfact, p - 2, p);

    /* pow2[i] = 2^i mod p, precomputed */
    int max_edges = n * (n - 1) / 2;
    u64 *pow2 = (u64*)malloc((max_edges + 1) * sizeof(u64));
    pow2[0] = 1;
    for (int i = 1; i <= max_edges; i++)
        pow2[i] = mulmod(pow2[i-1], 2, p);

    /* Recursive enumeration of odd-part partitions */
    /* Track: placed parts as (part_index, multiplicity) pairs */
    int placed_idx[100], placed_mult[100];
    int depth = 0;
    u64 total = 0;

    /* Stack-based enumeration to avoid recursion overhead */
    /* For each part index pi (descending), try multiplicity 0, 1, 2, ... */

    /* Simple recursive version first */
    void enumerate(int remaining, int max_pi) {
        if (remaining == 0) {
            /* Compute e(λ) = edges and z(λ) mod p */
            i64 edges = 0;
            u64 z_mod = 1;

            for (int i = 0; i < depth; i++) {
                int ki = odd_parts[placed_idx[i]];
                int mi = placed_mult[i];

                /* Self-cycle: mi * (ki-1) / 2 */
                edges += (i64)mi * (ki - 1) / 2;
                /* Same-size pairs: C(mi,2) * ki */
                edges += (i64)mi * (mi - 1) / 2 * ki;
                /* z contribution: ki^mi * mi! */
                z_mod = mulmod(z_mod, powmod(ki, mi, p), p);
                for (int f = 1; f <= mi; f++)
                    z_mod = mulmod(z_mod, f, p);

                /* Cross pairs */
                for (int j = i + 1; j < depth; j++) {
                    int kj = odd_parts[placed_idx[j]];
                    int mj = placed_mult[j];
                    edges += (i64)mi * mj * gcd_tab[placed_idx[i]][placed_idx[j]];
                }
            }

            if (edges > max_edges) edges = max_edges; /* safety */
            if (z_mod == 0) return; /* p | z, skip */

            u64 z_inv = powmod(z_mod, p - 2, p);
            u64 contrib = mulmod(mulmod(nfact, z_inv, p), pow2[edges], p);
            total = (total + contrib) % p;
            return;
        }

        for (int pi = max_pi; pi >= 0; pi--) {
            int k = odd_parts[pi];
            if (k > remaining) continue;

            int max_m = remaining / k;
            for (int m = 1; m <= max_m; m++) {
                placed_idx[depth] = pi;
                placed_mult[depth] = m;
                depth++;
                enumerate(remaining - m * k, pi - 1);
                depth--;
            }
        }
    }

    enumerate(n, num_odd - 1);
    free(pow2);

    /* Result: total = n! * a(n) mod p */
    return mulmod(total, nfact_inv, p);
}

/* Generate primes near 2^62 */
static int is_prime(u64 n) {
    if (n < 2) return 0;
    if (n < 4) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    for (u64 i = 5; i * i <= n; i += 6)
        if (n % i == 0 || n % (i+2) == 0) return 0;
    return 1;
}

static u64 next_prime(u64 n) {
    if (n < 2) return 2;
    n |= 1; /* make odd */
    while (!is_prime(n)) n += 2;
    return n;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n>\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    n_global = n;

    /* Estimate bits */
    double est_log2 = n * (n-1) / 2.0;
    for (int k = 2; k <= n; k++)
        est_log2 -= log2(k);
    int est_bits = (int)est_log2 + 20; /* safety margin */

    /* Choose primes */
    int prime_bits = 62; /* use large primes for fewer CRT steps */
    int num_primes = est_bits / prime_bits + 2;

    printf("n=%d: est %d bits, using %d primes of ~%d bits\n",
           n, est_bits, num_primes, prime_bits);

    u64 *primes = (u64*)malloc(num_primes * sizeof(u64));
    u64 *residues = (u64*)malloc(num_primes * sizeof(u64));

    /* Generate primes near 2^62 */
    u64 p = (1ULL << 62) + 1;
    for (int i = 0; i < num_primes; i++) {
        p = next_prime(p);
        primes[i] = p;
        p += 2;
    }

    /* Compute residues */
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (int i = 0; i < num_primes; i++) {
        residues[i] = burnside_mod_p(n, primes[i]);
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    printf("Modular computation: %.3fs\n", elapsed);

    /* CRT reconstruction using GMP */
    mpz_t result, modulus, temp, diff;
    mpz_init_set_ui(result, residues[0]);
    mpz_init_set_ui(modulus, primes[0]);
    mpz_init(temp);
    mpz_init(diff);

    for (int i = 1; i < num_primes; i++) {
        /* result + modulus * t where t ≡ (r - result) * modulus^{-1} (mod p_i) */
        u64 r = residues[i];
        u64 pi = primes[i];

        u64 result_mod_pi = mpz_fdiv_ui(result, pi);
        u64 d = (r + pi - result_mod_pi) % pi;
        u64 mod_inv = powmod(mpz_fdiv_ui(modulus, pi), pi - 2, pi);
        u64 t_val = mulmod(d, mod_inv, pi);

        mpz_mul_ui(temp, modulus, t_val);
        mpz_add(result, result, temp);
        mpz_mul_ui(modulus, modulus, pi);
    }

    /* Print result */
    printf("a(%d) = ", n);
    mpz_out_str(stdout, 10, result);
    printf("\n");

    clock_gettime(CLOCK_MONOTONIC, &t1);
    elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    printf("Total time: %.3fs\n", elapsed);

    /* Cleanup */
    mpz_clear(result);
    mpz_clear(modulus);
    mpz_clear(temp);
    mpz_clear(diff);
    free(primes);
    free(residues);

    return 0;
}
