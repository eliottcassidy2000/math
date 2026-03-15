/*
 * nud_weight.c — Compute W(n) = Σ_{σ ∈ NUD(n)} 2^{adj1(σ)}
 * opus-2026-03-14-S89c
 *
 * NUD(n) = permutations of [0..n-1] with no unit descent
 * adj1(σ) = number of positions j with σ(j+1) = σ(j) + 1 (unit ascent)
 *
 * Uses transfer matrix DP:
 * State: (last_value, set_of_used_values)
 * Transition: extend by one position, track unit ascents
 *
 * dp[mask][v] = Σ 2^{adj1 so far} over all partial perms ending at v
 * using exactly the vertices in mask, with no unit descent.
 *
 * Memory: 2^n × n × 16 bytes (using __int128 for large counts)
 * n=20: 2^20 × 20 × 16 = 320 MB — feasible
 * n=25: 2^25 × 25 × 16 = 12.8 GB — tight
 *
 * Compile: gcc -O3 -o nud_weight nud_weight.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

/* Use unsigned __int128 for large counts */
typedef unsigned __int128 uint128;

void print_uint128(uint128 val) {
    if (val == 0) { printf("0"); return; }
    char buf[50];
    int pos = 0;
    while (val > 0) {
        buf[pos++] = '0' + (int)(val % 10);
        val /= 10;
    }
    for (int i = pos - 1; i >= 0; i--)
        putchar(buf[i]);
}

/* Compute n! as uint128 */
uint128 factorial128(int n) {
    uint128 f = 1;
    for (int i = 2; i <= n; i++) f *= i;
    return f;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    if (n < 2 || n > 25) {
        fprintf(stderr, "n must be 2..25\n");
        return 1;
    }

    long long n_states = 1LL << n;
    double mem_gb = (double)n_states * n * sizeof(uint128) / (1024.0*1024*1024);
    printf("n=%d, states=2^%d=%lld, memory=%.2f GB\n", n, n, n_states, mem_gb);
    if (mem_gb > 16.0) {
        fprintf(stderr, "Too much memory\n");
        return 1;
    }

    /* dp[mask * n + v] = sum of 2^{adj1} over partial NUD perms
     * using vertices in mask, ending at v */
    uint128 *dp = (uint128 *)calloc(n_states * n, sizeof(uint128));
    if (!dp) {
        fprintf(stderr, "Allocation failed\n");
        return 1;
    }

    /* Initialize: single-vertex paths */
    for (int v = 0; v < n; v++) {
        dp[((long long)1 << v) * n + v] = 1;
    }

    time_t t0 = time(NULL);
    int last_pc = 0;

    for (long long mask = 1; mask < n_states; mask++) {
        int pc = __builtin_popcountll(mask);
        if (pc > last_pc) {
            time_t t1 = time(NULL);
            printf("  level %d/%d, %lds\n", last_pc, n, t1 - t0);
            fflush(stdout);
            last_pc = pc;
        }

        for (int v = 0; v < n; v++) {
            if (!(mask & (1LL << v))) continue;
            uint128 cnt = dp[mask * n + v];
            if (cnt == 0) continue;

            /* Try extending with vertex u not in mask */
            for (int u = 0; u < n; u++) {
                if (mask & (1LL << u)) continue;

                /* Check: u must NOT be v-1 (that would be a unit descent) */
                if (u == v - 1) continue;

                /* Weight: multiply by 2 if u = v+1 (unit ascent) */
                uint128 weight = (u == v + 1) ? 2 * cnt : cnt;

                dp[(mask | (1LL << u)) * n + u] += weight;
            }
        }
    }

    time_t t2 = time(NULL);
    printf("  level %d/%d, %lds\n", n, n, t2 - t0);

    /* Sum over all endpoints for full mask */
    uint128 W = 0;
    uint128 NUD_count = 0;
    long long full = n_states - 1;

    /* We also need NUD count (without 2^adj1 weighting).
     * To get NUD count, we'd need a separate DP without the 2× factor.
     * But we can compute it from A000255 directly. */
    for (int v = 0; v < n; v++) {
        W += dp[full * n + v];
    }

    printf("\nResults for n=%d (elapsed %lds):\n", n, t2 - t0);
    printf("  W(n) = "); print_uint128(W); printf("\n");

    uint128 nf = factorial128(n);
    printf("  n! = "); print_uint128(nf); printf("\n");

    /* W/n! as decimal */
    /* CV² = W/n! - 1 */
    /* Since W and n! can be huge, compute ratio carefully */
    /* W = q * n! + r, then CV² = q + r/n! - 1 */
    uint128 q = W / nf;
    uint128 r = W % nf;
    printf("  W/n! = "); print_uint128(q); printf(" + ");
    print_uint128(r); printf("/"); print_uint128(nf); printf("\n");

    /* n * CV² = n * W/n! - n */
    /* For floating point approximation */
    double W_approx = 0;
    for (int v = 0; v < n; v++) {
        double val = (double)dp[full * n + v];
        W_approx += val;
    }
    double cv2 = W_approx / (double)nf - 1.0;
    printf("  CV² ≈ %.12f\n", cv2);
    printf("  n×CV² ≈ %.12f\n", n * cv2);
    printf("  2-n×CV² ≈ %.12f\n", 2.0 - n * cv2);
    printf("  n²(2/n-CV²) ≈ %.8f\n", (double)n*n * (2.0/n - cv2));

    free(dp);
    return 0;
}
