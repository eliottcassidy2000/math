
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 19
#define FULL ((1 << N) - 1)

// dp[mask][last] = number of Hamiltonian paths using vertices in mask, ending at last
// We use 64-bit integers since the count can be large
// Memory: 2^19 * 19 * 8 bytes = ~80 MB -- feasible

static long long dp[1 << N][N];

// adj_out[v] = bitmask of vertices that v has an edge TO
static const int adj_out[N] = {
    0x30af2,  /* vertex 0 */
    0x615e4,  /* vertex 1 */
    0x42bc9,  /* vertex 2 */
    0x05793,  /* vertex 3 */
    0x0af26,  /* vertex 4 */
    0x15e4c,  /* vertex 5 */
    0x2bc98,  /* vertex 6 */
    0x57930,  /* vertex 7 */
    0x2f261,  /* vertex 8 */
    0x5e4c2,  /* vertex 9 */
    0x3c985,  /* vertex 10 */
    0x7930a,  /* vertex 11 */
    0x72615,  /* vertex 12 */
    0x64c2b,  /* vertex 13 */
    0x49857,  /* vertex 14 */
    0x130af,  /* vertex 15 */
    0x2615e,  /* vertex 16 */
    0x4c2bc,  /* vertex 17 */
    0x18579  /* vertex 18 */
};

// adj_in[v] = bitmask of vertices that have an edge TO v
static const int adj_in[N] = {
    0x4f50c,  /* vertex 0 */
    0x1ea19,  /* vertex 1 */
    0x3d432,  /* vertex 2 */
    0x7a864,  /* vertex 3 */
    0x750c9,  /* vertex 4 */
    0x6a193,  /* vertex 5 */
    0x54327,  /* vertex 6 */
    0x2864f,  /* vertex 7 */
    0x50c9e,  /* vertex 8 */
    0x2193d,  /* vertex 9 */
    0x4327a,  /* vertex 10 */
    0x064f5,  /* vertex 11 */
    0x0c9ea,  /* vertex 12 */
    0x193d4,  /* vertex 13 */
    0x327a8,  /* vertex 14 */
    0x64f50,  /* vertex 15 */
    0x49ea1,  /* vertex 16 */
    0x13d43,  /* vertex 17 */
    0x27a86  /* vertex 18 */
};


int main() {
    memset(dp, 0, sizeof(dp));

    // Base case: paths of length 1 (single vertex)
    for (int v = 0; v < N; v++) {
        dp[1 << v][v] = 1;
    }

    // Fill DP
    for (int mask = 1; mask <= FULL; mask++) {
        int popcount = __builtin_popcount(mask);
        if (popcount < 2) continue;  // need at least 2 vertices for a transition

        for (int last = 0; last < N; last++) {
            if (!(mask & (1 << last))) continue;  // last must be in mask

            // Predecessors of last that are in mask
            int prev_mask = mask ^ (1 << last);  // mask without last
            int candidates = adj_in[last] & prev_mask;  // vertices in prev_mask that point to last

            while (candidates) {
                int prev = __builtin_ctz(candidates);  // lowest set bit
                dp[mask][last] += dp[prev_mask][prev];
                candidates &= candidates - 1;  // clear lowest bit
            }
        }
    }

    // Sum over all ending vertices for the full mask
    long long total = 0;
    for (int v = 0; v < N; v++) {
        total += dp[FULL][v];
    }

    printf("H(T_%d) = %lld\n", N, total);

    // Also print per-endpoint counts for verification
    printf("\nPer-endpoint counts:\n");
    for (int v = 0; v < N; v++) {
        printf("  endpoint %d: %lld\n", v, dp[FULL][v]);
    }

    return 0;
}
