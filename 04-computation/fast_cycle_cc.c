/*
 * fast_cycle_cc.c — Fast C implementation of cycle_cc for tournament α-decomposition
 *
 * Computes cc[mask] = number of directed odd cycles with vertex-set exactly `mask`
 * in a tournament on n vertices (n ≤ 25).
 *
 * Algorithm:
 *   For each mask (in increasing integer order, which guarantees sub-before-super):
 *     Let s = lsb(mask) = minimum vertex (canonical start)
 *     For each v in mask:
 *       p = paths[mask * n + v]  // #paths from s to v using exactly vertices in mask
 *       if popcount(mask) >= 3 and odd and adj[v] has bit s:
 *         cc[mask] += p          // this closes an odd cycle back to s
 *       For each u in adj[v] & ~mask & hi_mask(s):
 *         paths[(mask|u_bit) * n + u] += p  // extend the path
 *
 * Memory:
 *   paths: 2^n * n * 8 bytes  (n=21: 352 MB, n=23: 1.47 GB, n=25: 6.4 GB)
 *   cc:    2^n * 8 bytes       (n=21: 16 MB, n=23: 64 MB, n=25: 256 MB)
 *
 * Output: binary file cc_n<n>.bin containing int64_t array of 2^n entries.
 *         Python can load this and run SSC separately.
 *
 * Usage:
 *   gcc -O3 -march=native -o fast_cycle_cc fast_cycle_cc.c
 *   ./fast_cycle_cc 23 circulant    # cyclic interval S={1,...,(n-1)/2}
 *   ./fast_cycle_cc 23 paley        # Paley QR tournament (n must be prime ≡ 3 mod 4)
 *
 * Output file: 05-knowledge/results/cc_n<n>_<type>.bin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define MAXN 25

int n;
uint32_t adj[MAXN];  /* adj[v] = bitmask of vertices u such that v -> u */

int64_t *paths_buf;  /* paths_buf[mask * n + v] */
int64_t *cc_buf;     /* cc_buf[mask] */

/* Build circulant tournament: adj[v] = {(v+k) mod n : k in S} */
void build_circulant(int *S, int S_len) {
    memset(adj, 0, sizeof(adj));
    for (int v = 0; v < n; v++) {
        for (int i = 0; i < S_len; i++) {
            int u = (v + S[i]) % n;
            adj[v] |= (1u << u);
        }
    }
}

/* Build Paley tournament: adj[v] has u iff (u-v) is a quadratic residue mod n */
int is_QR(int a, int p) {
    /* Euler's criterion: a^{(p-1)/2} ≡ 1 (mod p) */
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

void build_paley() {
    memset(adj, 0, sizeof(adj));
    for (int v = 0; v < n; v++) {
        for (int u = 0; u < n; u++) {
            if (u == v) continue;
            if (is_QR((u - v + n) % n, n)) {
                adj[v] |= (1u << u);
            }
        }
    }
}

void compute_cycle_cc() {
    uint64_t N = 1ULL << n;

    /* Clear arrays (calloc already zeroed) */
    /* Base case: paths from each s using just {s} */
    for (int s = 0; s < n; s++) {
        uint64_t s_mask = 1ULL << s;
        paths_buf[s_mask * (uint64_t)n + s] = 1;
    }

    /* Process all masks in increasing integer order
     * (guarantees sub-before-super: m ⊂ m' => m < m') */
    for (uint64_t mask = 1; mask < N; mask++) {
        int s = __builtin_ctzll(mask);           /* min vertex */
        int L = __builtin_popcountll(mask);       /* cycle length if closed */
        uint32_t hi = ~((1u << s) - 1u);         /* bits >= s */
        uint32_t cands_base = (uint32_t)(mask & hi); /* vertices in mask that are >= s */

        /* Iterate over all v in mask (only those >= s matter, since s = min) */
        uint32_t vmask = (uint32_t)mask;
        while (vmask) {
            uint32_t vbit = vmask & (-vmask);
            int v = __builtin_ctz(vbit);
            vmask ^= vbit;

            int64_t p = paths_buf[mask * (uint64_t)n + v];
            if (p == 0) continue;

            /* Try to close an odd cycle back to s */
            if (L >= 3 && (L & 1) == 1 && ((adj[v] >> s) & 1)) {
                cc_buf[mask] += p;
            }

            /* Extend path to u > s, not in mask, reachable from v */
            uint32_t cands = adj[v] & ~((uint32_t)mask) & hi;
            while (cands) {
                uint32_t ubit = cands & (-cands);
                int u = __builtin_ctz(ubit);
                paths_buf[(mask | ubit) * (uint64_t)n + u] += p;
                cands ^= ubit;
            }
        }
    }
}

int main(int argc, char *argv[]) {
    n = (argc > 1) ? atoi(argv[1]) : 23;
    const char *type = (argc > 2) ? argv[2] : "circulant";

    if (n > MAXN) {
        fprintf(stderr, "n=%d exceeds MAXN=%d\n", n, MAXN);
        return 1;
    }

    uint64_t N = 1ULL << n;
    size_t paths_bytes = (size_t)N * n * sizeof(int64_t);
    size_t cc_bytes    = (size_t)N * sizeof(int64_t);

    printf("n=%d, type=%s, N=%llu\n", n, type, (unsigned long long)N);
    printf("Allocating paths (%.0f MB) + cc (%.0f MB)...\n",
           (double)paths_bytes / (1 << 20),
           (double)cc_bytes    / (1 << 20));
    fflush(stdout);

    paths_buf = (int64_t *)calloc(N * (uint64_t)n, sizeof(int64_t));
    cc_buf    = (int64_t *)calloc(N, sizeof(int64_t));

    if (!paths_buf || !cc_buf) {
        fprintf(stderr, "Memory allocation failed (need %.1f GB)\n",
                (double)(paths_bytes + cc_bytes) / (1ULL << 30));
        return 1;
    }

    /* Build tournament */
    if (strcmp(type, "paley") == 0) {
        /* Check preconditions */
        if (n % 4 != 3) {
            fprintf(stderr, "Paley requires n ≡ 3 (mod 4)\n");
            return 1;
        }
        build_paley();
        printf("Built Paley tournament (QR mod %d)\n", n);
    } else {
        /* Default: cyclic interval S = {1, 2, ..., (n-1)/2} */
        int S_len = (n - 1) / 2;
        int S[MAXN];
        for (int i = 0; i < S_len; i++) S[i] = i + 1;
        build_circulant(S, S_len);
        printf("Built cyclic interval tournament S={1..%d}\n", S_len);
    }

    printf("Running cycle_cc...\n");
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
             "05-knowledge/results/cc_n%d_%s.bin", n, type);
    FILE *f = fopen(fname, "wb");
    if (!f) {
        fprintf(stderr, "Cannot open %s for writing\n", fname);
        return 1;
    }
    fwrite(cc_buf, sizeof(int64_t), N, f);
    fclose(f);
    printf("Saved cc[0..%llu] to %s (%.0f MB)\n",
           (unsigned long long)(N-1), fname, (double)cc_bytes / (1 << 20));

    free(paths_buf);
    free(cc_buf);
    return 0;
}
