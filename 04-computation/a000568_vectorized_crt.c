/*
 * a000568_vectorized_crt.c — opus-2026-04-04-S23
 *
 * KEY INNOVATION: Enumerate partitions ONCE, accumulate residues
 * for ALL primes simultaneously. This amortizes the partition
 * enumeration cost across all primes.
 *
 * For each partition λ:
 *   For each prime p_i:
 *     bucket_i[e(λ)] += (n!/z_λ) mod p_i
 *
 * Then merge buckets and do CRT reconstruction.
 *
 * Cost: O(partitions × num_primes) for the inner loop
 * vs O(partitions × num_primes) for separate runs.
 * BUT: the inner loop is a single mulmod per prime (vs full recursion).
 *
 * Compile: gcc -O3 -o a000568_vec_crt a000568_vectorized_crt.c -lgmp -lm
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
static u64 powmod(u64 b, u64 e, u64 m) {
    u64 r = 1; b %= m;
    while (e > 0) { if (e&1) r = mulmod(r,b,m); b = mulmod(b,b,m); e >>= 1; }
    return r;
}
static int gcd_int(int a, int b) {
    while (b) { int t=b; b=a%b; a=t; } return a;
}

/* Global state */
static int odd_parts[500], num_odd;
static int gcd_tab[500][500];
static int pl_idx[100], pl_mul[100], depth;

/* Multi-prime state */
static int np_global;       /* number of primes */
static u64 *primes_g;       /* array of primes */
static u64 *nfact_g;        /* n! mod p_i */
static u64 **pow2_g;        /* pow2_g[i][j] = 2^j mod p_i */
static u64 **buckets_g;     /* buckets_g[i][e] = partial sum mod p_i */
static int max_e_g;

static void enum_vec(int rem, int mpi) {
    if (rem == 0) {
        /* Compute edges and z for this partition */
        i64 edges = 0;

        /* Compute z as product of (k^m * m!) — need this mod each prime */
        /* We'll compute z mod p_i for each prime */
        for (int i = 0; i < depth; i++) {
            int k = odd_parts[pl_idx[i]], m = pl_mul[i];
            edges += (i64)m * (k - 1) / 2;
            edges += (i64)m * (m - 1) / 2 * k;
            for (int j = i + 1; j < depth; j++)
                edges += (i64)m * pl_mul[j] * gcd_tab[pl_idx[i]][pl_idx[j]];
        }
        if (edges > max_e_g) edges = max_e_g;

        /* For each prime: compute n!/z mod p and add to bucket[edges] */
        for (int pi = 0; pi < np_global; pi++) {
            u64 p = primes_g[pi];
            u64 z = 1;
            for (int i = 0; i < depth; i++) {
                int k = odd_parts[pl_idx[i]], m = pl_mul[i];
                z = mulmod(z, powmod(k, m, p), p);
                for (int f = 1; f <= m; f++) z = mulmod(z, f, p);
            }
            if (z == 0) continue;
            u64 zi = powmod(z, p - 2, p);
            u64 c = mulmod(nfact_g[pi], zi, p);
            buckets_g[pi][edges] = (buckets_g[pi][edges] + c) % p;
        }
        return;
    }
    for (int pi = mpi; pi >= 0; pi--) {
        int k = odd_parts[pi];
        if (k > rem) continue;
        for (int m = 1; m * k <= rem; m++) {
            pl_idx[depth] = pi; pl_mul[depth] = m; depth++;
            enum_vec(rem - m * k, pi - 1);
            depth--;
        }
    }
}

/* Prime generation */
static int is_prime64(u64 n) {
    if (n < 2) return 0; if (n < 4) return 1;
    if (n%2==0 || n%3==0) return 0;
    for (u64 i=5; i*i<=n; i+=6) if (n%i==0 || n%(i+2)==0) return 0;
    return 1;
}
static u64 next_prime64(u64 n) { n|=1; while(!is_prime64(n)) n+=2; return n; }

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <n>\n", argv[0]); return 1; }
    int n = atoi(argv[1]);

    double el2 = n*(n-1)/2.0;
    for (int k=2; k<=n; k++) el2 -= log2(k);
    int ebits = (int)el2 + 20;
    int np = ebits / 62 + 2;
    if (np < 2) np = 2;
    np_global = np;

    /* Setup */
    num_odd = 0;
    for (int k=1; k<=n; k+=2) odd_parts[num_odd++] = k;
    for (int i=0; i<num_odd; i++)
        for (int j=0; j<num_odd; j++)
            gcd_tab[i][j] = gcd_int(odd_parts[i], odd_parts[j]);

    int me = n*(n-1)/2;
    max_e_g = me;

    primes_g = malloc(np * sizeof(u64));
    nfact_g = malloc(np * sizeof(u64));
    pow2_g = malloc(np * sizeof(u64*));
    buckets_g = malloc(np * sizeof(u64*));

    u64 p = (1ULL<<62)+1;
    for (int i=0; i<np; i++) {
        p = next_prime64(p); primes_g[i] = p; p += 2;
        nfact_g[i] = 1;
        for (int k=2; k<=n; k++) nfact_g[i] = mulmod(nfact_g[i], k, primes_g[i]);
        pow2_g[i] = calloc(me+1, sizeof(u64));
        pow2_g[i][0] = 1;
        for (int j=1; j<=me; j++) pow2_g[i][j] = mulmod(pow2_g[i][j-1], 2, primes_g[i]);
        buckets_g[i] = calloc(me+1, sizeof(u64));
    }

    printf("n=%d, %d primes, max_edges=%d\n", n, np, me);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* SINGLE partition enumeration, ALL primes accumulated */
    depth = 0;
    enum_vec(n, num_odd - 1);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double sec_enum = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;

    /* Merge buckets: residue[i] = Σ_e bucket[i][e] * 2^e mod p_i */
    u64 *resids = malloc(np * sizeof(u64));
    for (int i=0; i<np; i++) {
        u64 s = 0;
        for (int e=0; e<=me; e++)
            s = (s + mulmod(buckets_g[i][e], pow2_g[i][e], primes_g[i])) % primes_g[i];
        /* Multiply by (n!)^{-1} */
        resids[i] = mulmod(s, powmod(nfact_g[i], primes_g[i]-2, primes_g[i]), primes_g[i]);
    }

    /* CRT */
    mpz_t result, modulus, temp;
    mpz_init_set_ui(result, resids[0]);
    mpz_init_set_ui(modulus, primes_g[0]);
    mpz_init(temp);
    for (int i=1; i<np; i++) {
        u64 r = resids[i], pi = primes_g[i];
        u64 rm = mpz_fdiv_ui(result, pi);
        u64 d = (r+pi-rm)%pi;
        u64 mi = powmod(mpz_fdiv_ui(modulus,pi), pi-2, pi);
        u64 tv = mulmod(d, mi, pi);
        mpz_mul_ui(temp, modulus, tv);
        mpz_add(result, result, temp);
        mpz_mul_ui(modulus, modulus, pi);
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double sec_total = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;

    printf("Enum: %.3fs, Total: %.3fs\n", sec_enum, sec_total);
    printf("a(%d) = ", n);
    mpz_out_str(stdout, 10, result);
    printf("\n");

    /* Cleanup */
    mpz_clear(result); mpz_clear(modulus); mpz_clear(temp);
    for (int i=0; i<np; i++) { free(pow2_g[i]); free(buckets_g[i]); }
    free(primes_g); free(nfact_g); free(pow2_g); free(buckets_g); free(resids);
    return 0;
}
