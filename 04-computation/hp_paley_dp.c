/*
 * hp_paley_dp.c — Compute H(P_p) for Paley tournament via bitmask DP
 * opus-2026-03-14-S89c
 *
 * Usage: ./hp_paley_dp <p>
 * p must be prime ≡ 3 mod 4
 *
 * DP state: dp[mask][v] = number of Hamiltonian paths ending at v
 * using exactly the vertices in mask.
 * Memory: 2^p × p × 8 bytes. p=29: 2^29 × 29 × 8 = 125 GB — TOO BIG.
 * p=23: 2^23 × 23 × 8 = 1.5 GB — feasible.
 * p=29: need meet-in-the-middle or other tricks.
 *
 * For p ≤ 23, this works directly.
 * For p = 29, we'd need ~125 GB which is impractical.
 *
 * Compile: gcc -O3 -o hp_paley_dp hp_paley_dp.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

/* Legendre symbol (a/p) via Euler criterion */
int legendre(int a, int p) {
    a = ((a % p) + p) % p;
    if (a == 0) return 0;
    int result = 1;
    int base = a;
    int exp = (p - 1) / 2;
    base %= p;
    while (exp > 0) {
        if (exp & 1) result = (int)((long long)result * base % p);
        exp >>= 1;
        base = (int)((long long)base * base % p);
    }
    return result == 1 ? 1 : -1;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <p>\n", argv[0]);
        return 1;
    }
    int p = atoi(argv[1]);

    /* Verify p ≡ 3 mod 4 */
    if (p % 4 != 3) {
        fprintf(stderr, "p must be ≡ 3 mod 4 for Paley tournament\n");
        return 1;
    }

    printf("Computing H(P_%d) via bitmask DP\n", p);

    /* Build adjacency: adj[i] is bitmask of vertices i beats */
    uint32_t *adj = (uint32_t *)calloc(p, sizeof(uint32_t));
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            if (i != j && legendre(j - i, p) == 1) {
                adj[i] |= (1u << j);
            }
        }
    }

    /* Check adjacency */
    int total_edges = 0;
    for (int i = 0; i < p; i++)
        total_edges += __builtin_popcount(adj[i]);
    printf("Total directed edges: %d (expected %d)\n", total_edges, p * (p - 1) / 2);

    /* DP: dp[mask][v] = # of Ham paths through vertices in mask ending at v */
    long long full_mask = (1LL << p) - 1;
    long long n_states = 1LL << p;

    /* Memory estimate */
    double mem_gb = (double)n_states * p * sizeof(long long) / (1024.0 * 1024 * 1024);
    printf("Memory needed: %.2f GB\n", mem_gb);
    if (mem_gb > 16.0) {
        fprintf(stderr, "Too much memory (%.1f GB > 16 GB limit)\n", mem_gb);
        return 1;
    }

    /* Allocate dp arrays: current and next aren't practical for full table.
     * Instead, iterate over masks in order of popcount.
     * dp[v] for current mask. But we need ALL masks simultaneously...
     *
     * Standard approach: dp[mask][v] stored as dp[mask * p + v].
     * For p=23: 2^23 * 23 * 8 = 1.38 GB — tight but OK.
     * For p=29: 2^29 * 29 * 8 = 125 GB — impossible.
     */

    long long *dp = (long long *)calloc(n_states * p, sizeof(long long));
    if (!dp) {
        fprintf(stderr, "Failed to allocate %.1f GB\n", mem_gb);
        return 1;
    }

    /* Initialize: single vertex paths */
    for (int v = 0; v < p; v++) {
        dp[(1LL << v) * p + v] = 1;
    }

    printf("Running DP...\n");
    time_t t0 = time(NULL);
    int last_pc = 0;
    long long masks_done = 0;

    /* Iterate over all masks in order */
    for (long long mask = 1; mask < n_states; mask++) {
        int pc = __builtin_popcountll(mask);
        if (pc > last_pc) {
            time_t t1 = time(NULL);
            printf("  level %d/%d done, elapsed %lds, masks=%lld\n",
                   last_pc, p, t1 - t0, masks_done);
            fflush(stdout);
            last_pc = pc;
        }
        masks_done++;

        for (int v = 0; v < p; v++) {
            if (!(mask & (1LL << v))) continue;
            long long cnt = dp[mask * p + v];
            if (cnt == 0) continue;

            /* Extend to neighbor u not in mask */
            uint32_t neighbors = adj[v] & ~(uint32_t)mask;
            while (neighbors) {
                int u = __builtin_ctz(neighbors);
                neighbors &= neighbors - 1;
                dp[(mask | (1LL << u)) * p + u] += cnt;
            }
        }
    }

    time_t t2 = time(NULL);
    printf("DP done in %lds\n", t2 - t0);

    /* Sum over all endpoints */
    long long H = 0;
    for (int v = 0; v < p; v++) {
        long long hv = dp[full_mask * p + v];
        H += hv;
        if (p <= 23)
            printf("  h(*, %d) = %lld\n", v, hv);
    }

    printf("\nH(P_%d) = %lld\n", p, H);
    printf("H / p = %lld\n", H / p);
    printf("H mod p = %lld\n", H % p);
    printf("(p-1)/2 = %d\n", (p - 1) / 2);
    if (H % p == 0) {
        long long h0 = H / p;
        int half = (p - 1) / 2;
        if (h0 % half == 0)
            printf("H / (p × (p-1)/2) = %lld\n", h0 / half);
    }

    /* Factor H */
    printf("\nFactoring H...\n");
    long long rem = H;
    for (long long d = 2; d * d <= rem; d++) {
        while (rem % d == 0) {
            printf("  %lld", d);
            rem /= d;
        }
    }
    if (rem > 1) printf("  %lld", rem);
    printf("\n");

    free(dp);
    free(adj);
    return 0;
}
