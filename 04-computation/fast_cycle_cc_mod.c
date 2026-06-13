/*
 * fast_cycle_cc_mod.c — Mod-prime version of fast_cycle_cc for large n
 *
 * Like fast_cycle_cc.c but uses int32_t throughout (mod a prime P < 2^31),
 * cutting the paths array from 8 bytes/entry to 4 bytes/entry.
 *
 * Memory for n=25:
 *   Original (int64): paths = 6.4 GB, cc = 0.25 GB  → 6.65 GB total
 *   This version (int32 mod P): paths = 3.2 GB, cc = 0.125 GB → 3.33 GB total
 *
 * Usage:
 *   gcc -O3 -march=native -o fast_cycle_cc_mod fast_cycle_cc_mod.c
 *   ./fast_cycle_cc_mod 25 circulant 2147483647
 *   ./fast_cycle_cc_mod 25 circulant 2147483629
 *   ./fast_cycle_cc_mod 25 paley 2147483647
 *
 * Output: 05-knowledge/results/cc_n<n>_<type>_p<prime>.i32
 *         (int32_t array of 2^n entries, values in [0, P-1])
 *
 * To recover exact cc values: run twice with two primes, apply CRT in Python.
 * Or directly pipe into the SSC pipeline (alpha_from_cc_bin.py accepts --mod mode).
 *
 * Primes to use:
 *   P1 = 2147483647  (2^31 - 1, Mersenne prime)
 *   P2 = 2147483629  (2^31 - 19)
 *   CRT covers values up to P1 * P2 ≈ 4.6e18 (sufficient for alpha_k)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define MAXN 25

int n;
int32_t prime;
uint32_t adj[MAXN];  /* adj[v] = bitmask of out-neighbors */

uint32_t *paths_buf;  /* paths_buf[mask * n + v] mod prime, values in [0, prime-1] */
uint32_t *cc_buf;     /* cc_buf[mask] mod prime, values in [0, prime-1] */

/* Euler's criterion: is a a quadratic residue mod p? */
int is_QR(int a, int p) {
    a = ((a % p) + p) % p;
    if (a == 0) return 0;
    long long base = a, exp = (p-1)/2, mod = p, result = 1;
    while (exp > 0) {
        if (exp & 1) result = result * base % mod;
        base = base * base % mod;
        exp >>= 1;
    }
    return result == 1;
}

void build_circulant(void) {
    int S_len = (n - 1) / 2;
    memset(adj, 0, sizeof(adj));
    for (int v = 0; v < n; v++)
        for (int k = 1; k <= S_len; k++)
            adj[v] |= 1u << ((v + k) % n);
}

void build_paley(void) {
    memset(adj, 0, sizeof(adj));
    for (int v = 0; v < n; v++)
        for (int u = 0; u < n; u++)
            if (u != v && is_QR((u - v + n) % n, n))
                adj[v] |= 1u << u;
}

void compute_cycle_cc(void) {
    uint64_t N = 1ULL << n;
    uint32_t P = (uint32_t)prime;  /* prime < 2^31, safe as uint32_t */

    /* Seed: each vertex s seeds a path of length 1 at mask = {s}, v = s */
    for (int s = 0; s < n; s++) {
        uint64_t s_mask = 1ULL << s;
        paths_buf[s_mask * (uint64_t)n + s] = 1;
    }

    for (uint64_t mask = 1; mask < N; mask++) {
        int s = __builtin_ctzll(mask);
        int L = __builtin_popcountll(mask);
        uint32_t hi = ~((1u << s) - 1u);

        uint32_t vmask = (uint32_t)mask;
        while (vmask) {
            uint32_t vbit = vmask & (0u - vmask);  /* isolate LSB (two's complement) */
            int v = __builtin_ctz(vbit);
            vmask ^= vbit;

            uint32_t p = paths_buf[mask * (uint64_t)n + v];
            if (p == 0) continue;

            /* Close odd cycle back to s */
            if (L >= 3 && (L & 1) == 1 && ((adj[v] >> s) & 1)) {
                uint32_t sum = cc_buf[mask] + p;
                if (sum >= P) sum -= P;
                cc_buf[mask] = sum;
            }

            /* Extend path to u > s, u not in mask, reachable from v */
            uint32_t cands = adj[v] & ~((uint32_t)mask) & hi;
            while (cands) {
                uint32_t ubit = cands & (0u - cands);
                int u = __builtin_ctz(ubit);
                uint64_t nmask = mask | ubit;
                uint32_t sum = paths_buf[nmask * (uint64_t)n + u] + p;
                if (sum >= P) sum -= P;
                paths_buf[nmask * (uint64_t)n + u] = sum;
                cands ^= ubit;
            }
        }
    }
}

int main(int argc, char *argv[]) {
    n     = (argc > 1) ? atoi(argv[1]) : 23;
    const char *type = (argc > 2) ? argv[2] : "circulant";
    prime = (argc > 3) ? (int32_t)atol(argv[3]) : 2147483647;

    if (n > MAXN) { fprintf(stderr, "n=%d > MAXN=%d\n", n, MAXN); return 1; }
    if (n > 25)   { fprintf(stderr, "n > 25 not supported\n"); return 1; }

    uint64_t N = 1ULL << n;
    size_t paths_bytes = (size_t)N * n * sizeof(int32_t);
    size_t cc_bytes    = (size_t)N     * sizeof(int32_t);

    printf("n=%d, type=%s, prime=%d, N=%llu\n", n, type, prime, (unsigned long long)N);
    printf("Allocating paths (%.1f GB) + cc (%.0f MB)...\n",
           (double)paths_bytes / (1ULL << 30),
           (double)cc_bytes    / (1 << 20));
    fflush(stdout);

    paths_buf = (uint32_t *)calloc((size_t)N * n, sizeof(uint32_t));
    cc_buf    = (uint32_t *)calloc((size_t)N,     sizeof(uint32_t));

    if (!paths_buf || !cc_buf) {
        fprintf(stderr, "OOM: need %.2f GB\n",
                (double)(paths_bytes + cc_bytes) / (1ULL << 30));
        return 1;
    }

    if (strcmp(type, "paley") == 0) {
        if (n % 4 != 3) { fprintf(stderr, "Paley requires n ≡ 3 mod 4\n"); return 1; }
        build_paley();
        printf("Built Paley tournament (QR mod %d)\n", n);
    } else {
        build_circulant();
        printf("Built cyclic interval tournament S={1..%d}\n", (n-1)/2);
    }

    printf("Running cycle_cc mod %d...\n", prime);
    fflush(stdout);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    compute_cycle_cc();
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9;
    printf("cycle_cc done: %.2fs\n", elapsed);

    /* Save cc to binary file */
    char fname[256];
    snprintf(fname, sizeof(fname),
             "05-knowledge/results/cc_n%d_%s_p%d.i32", n, type, prime);
    FILE *f = fopen(fname, "wb");
    if (!f) { fprintf(stderr, "Cannot open %s\n", fname); return 1; }
    fwrite(cc_buf, sizeof(uint32_t), N, f);
    fclose(f);
    printf("Saved %llu int32 entries to %s (%.0f MB)\n",
           (unsigned long long)N, fname, (double)cc_bytes / (1 << 20));

    free(paths_buf);
    free(cc_buf);
    return 0;
}
