/*
 * h23_achievability.c — Is H=23 achievable? If so, at what n?
 *
 * Key insight from THM-079: cycle-rich n>=8 has H>=25.
 * So H=23 can only come from non-cycle-rich tournaments (which have
 * a removable vertex, reducing to smaller n).
 *
 * Actually: H=23 IS known to be achievable at n=6 and n=7.
 * From the n=7 spectrum: {1,3,5,9,11,...,23,...} includes 23.
 *
 * So this is NOT a gap. But what about H values near the boundaries
 * created by the cycle-rich min-H bounds?
 *
 * Cycle-rich min-H: n=8: 25, n=9: 45
 * Non-cycle-rich: has removable vertex, H(T) = H(T-v) for some v.
 * So H(T) can be ANY achievable H value at n-1.
 *
 * This means: H is achievable at n iff either:
 * (a) H is achievable by a cycle-rich n-tournament (H >= 25 for n=8), OR
 * (b) H is achievable at n-1 (via removable vertex induction)
 *
 * So the achievable H values at n are: achievable_at(n-1) UNION cycle_rich_H(n).
 * The achievable set is MONOTONICALLY GROWING in n.
 * A permanent gap must be absent from ALL n, i.e., never in cycle_rich_H(n)
 * for any n, AND never in achievable_at(n) for small n.
 *
 * Since achievable_at(7) includes {1,3,5,9,11,...,189} (77 values),
 * any value in this set is permanently achievable.
 *
 * A permanent gap must be: (a) not in achievable_at(7), AND
 * (b) not in cycle_rich_H(n) for any n >= 8.
 *
 * Known permanent gaps:
 * - H=7: not achievable at any n <= 7, not in any cycle_rich_H(n)
 * - H=21: not achievable at any n <= 8, not in any cycle_rich_H(n >= 8)
 *
 * Candidates for next gap: must be NOT in achievable_at(7).
 * From n=7 spectrum, gaps in [1,189] are:
 * 7, 21, 63, 107, 119, 149, 161-169, 173, 177-187
 *
 * 63 fills at n=8 (known). What about 107, 119, 149, etc.?
 *
 * Let's sample cycle-rich n=8 tournaments to see which of these fill.
 *
 * Author: opus-2026-03-07-S43
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 8

static unsigned int rng_state = 12345;
static inline unsigned int xorshift() {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

static inline void random_tournament(unsigned short out[N]) {
    memset(out, 0, N * sizeof(unsigned short));
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++) {
            if (xorshift() & 1)
                out[i] |= (1 << j);
            else
                out[j] |= (1 << i);
        }
}

static long long hamiltonian_paths(const unsigned short out[N]) {
    static long long dp[1 << N][N];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++)
        dp[1 << v][v] = 1;
    int full = (1 << N) - 1;
    for (int mask = 1; mask <= full; mask++)
        for (int v = 0; v < N; v++) {
            if (!(mask & (1 << v)) || dp[mask][v] == 0) continue;
            unsigned short nexts = out[v] & ~mask;
            while (nexts) {
                int u = __builtin_ctz(nexts);
                nexts &= nexts - 1;
                dp[mask | (1 << u)][u] += dp[mask][v];
            }
        }
    long long total = 0;
    for (int v = 0; v < N; v++) total += dp[full][v];
    return total;
}

#define MAX_H 700

int main() {
    setvbuf(stdout, NULL, _IONBF, 0);

    /* n=7 achievable H values (known) */
    int n7_achievable[] = {1,3,5,9,11,13,15,17,19,23,25,27,29,31,33,35,37,39,41,
        43,45,47,49,51,53,55,57,59,61,65,67,69,71,73,75,77,79,81,83,85,87,89,91,
        93,95,97,99,101,103,105,109,111,113,115,117,121,123,125,127,129,131,133,
        135,137,139,141,143,145,147,151,153,155,157,159,171,175,189};

    int found_at_n8[MAX_H + 1];
    memset(found_at_n8, 0, sizeof(found_at_n8));

    /* Mark n=7 achievable */
    for (int i = 0; i < (int)(sizeof(n7_achievable)/sizeof(int)); i++)
        if (n7_achievable[i] <= MAX_H)
            found_at_n8[n7_achievable[i]] = 1; /* inherited */

    long long trials = 100000000LL;
    unsigned short out[N];

    for (long long t = 0; t < trials; t++) {
        random_tournament(out);
        long long H = hamiltonian_paths(out);
        if (H <= MAX_H && !found_at_n8[(int)H]) {
            found_at_n8[(int)H] = 2; /* NEW at n=8 */
            printf("NEW at n=8: H=%lld (trial %lld)\n", H, t);
        }

        if ((t + 1) % 20000000 == 0) {
            int missing = 0;
            for (int h = 1; h <= MAX_H; h += 2)
                if (!found_at_n8[h]) missing++;
            fprintf(stderr, "t=%lldM, missing odd in [1,%d]: %d\n",
                    (t+1)/1000000, MAX_H, missing);
        }
    }

    printf("\n=== n=8 SPECTRUM (sampling + n=7 inheritance) ===\n");
    printf("Trials: %lld\n", trials);

    printf("\nStill missing odd values in [1, %d]:\n", MAX_H);
    for (int h = 1; h <= MAX_H; h += 2) {
        if (!found_at_n8[h])
            printf("  H=%d\n", h);
    }

    printf("\nNew values first appearing at n=8:\n");
    for (int h = 1; h <= MAX_H; h += 2) {
        if (found_at_n8[h] == 2)
            printf("  H=%d\n", h);
    }

    return 0;
}
