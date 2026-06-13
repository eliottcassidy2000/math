/*
 * simplicial_redei_n8.c — Exhaustive simplicial Rédei at n=8
 * opus-2026-03-15-S90b
 *
 * Highly optimized: bitmask encoding, fast transitive core via lookup,
 * early DAG/total-order termination.
 *
 * Tournament encoding: 28 bits for C(8,2) = 28 arcs.
 * Bit k corresponds to arc (i,j) with i<j, in lexicographic order.
 * Bit=1 means i→j, bit=0 means j→i.
 *
 * For each of 2^28 = 268,435,456 tournaments:
 *   1. Compute transitive core edges (check all C(8,3)=56 triples)
 *   2. Check if core is a DAG (topological sort)
 *   3. If DAG, check if total order (all pairs comparable)
 *   4. Count H for sim_H=1 tournaments
 *
 * Expected: 2·8! = 80640 with sim_H=1
 *   40320 transitive (H=1) + 40320 near-transitive (H=65 = 2^6+1)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

#define N 8
#define NPAIRS 28
#define TOTAL_TOURNAMENTS (1U << NPAIRS) /* 268435456 */

/* Arc index: for pair (i,j) with i<j, the bit position */
static int arc_idx[N][N];
static int arc_i[NPAIRS], arc_j[NPAIRS];

static void init_arc_indices(void) {
    int idx = 0;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++) {
            arc_idx[i][j] = idx;
            arc_idx[j][i] = idx;
            arc_i[idx] = i;
            arc_j[idx] = j;
            idx++;
        }
}

/* Check if i→j in tournament encoded as bits */
static inline int has_arc(uint32_t bits, int i, int j) {
    if (i < j) return (bits >> arc_idx[i][j]) & 1;
    else return !((bits >> arc_idx[j][i]) & 1);
}

/* All C(8,3) = 56 triples, precomputed */
static int triples[56][3];
static int ntriples;

static void init_triples(void) {
    ntriples = 0;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++)
            for (int k = j+1; k < N; k++) {
                triples[ntriples][0] = i;
                triples[ntriples][1] = j;
                triples[ntriples][2] = k;
                ntriples++;
            }
}

/* Process one tournament. Returns:
 * 0 = cyclic core (sim_H = 0)
 * 1 = total order (sim_H = 1)
 * 2 = DAG but not total order (sim_H > 1 — should not occur)
 */
static int process_tournament(uint32_t bits) {
    /* Core edges as adjacency bitmask per vertex */
    uint8_t core[N] = {0}; /* core[i] bit j set means i→j is in core */

    /* Check all 56 triples */
    for (int t = 0; t < ntriples; t++) {
        int i = triples[t][0], j = triples[t][1], k = triples[t][2];

        int ij = has_arc(bits, i, j);
        int jk = has_arc(bits, j, k);
        int ik = has_arc(bits, i, k);

        /* Check all 6 orderings for transitive triple */
        if (ij && jk && ik) {
            /* i→j→k, i→k: transitive (i > j > k) */
            core[i] |= (1 << j) | (1 << k);
            core[j] |= (1 << k);
        } else if (!ij && !jk && !ik) {
            /* j→i, k→j, k→i: transitive (k > j > i) */
            core[k] |= (1 << j) | (1 << i);
            core[j] |= (1 << i);
        } else if (ij && !jk && ik) {
            /* i→j, k→j, i→k: transitive (i > k > j) */
            core[i] |= (1 << k) | (1 << j);
            core[k] |= (1 << j);
        } else if (!ij && jk && !ik) {
            /* j→i, j→k, k→i: transitive (j > k > i) */
            core[j] |= (1 << k) | (1 << i);
            core[k] |= (1 << i);
        } else if (!ij && jk && ik) {
            /* j→i, j→k, i→k: transitive (j > i > k) */
            core[j] |= (1 << i) | (1 << k);
            core[i] |= (1 << k);
        } else if (ij && !jk && !ik) {
            /* i→j, k→j, k→i: transitive (k > i > j) */
            core[k] |= (1 << i) | (1 << j);
            core[i] |= (1 << j);
        }
        /* else: it's a 3-cycle (3 of 6 arcs match = cycle) */
    }

    /* DAG check via topological sort (Kahn's algorithm) */
    uint8_t in_deg[N] = {0};
    for (int i = 0; i < N; i++) {
        uint8_t c = core[i];
        while (c) {
            int j = __builtin_ctz(c);
            in_deg[j]++;
            c &= c - 1;
        }
    }

    uint8_t queue[N];
    int qfront = 0, qback = 0;
    for (int v = 0; v < N; v++)
        if (in_deg[v] == 0)
            queue[qback++] = v;

    int topo_count = 0;
    uint8_t topo[N];
    uint8_t in_d[N];
    memcpy(in_d, in_deg, sizeof(in_deg));

    while (qfront < qback) {
        int v = queue[qfront++];
        topo[topo_count++] = v;
        uint8_t c = core[v];
        while (c) {
            int u = __builtin_ctz(c);
            in_d[u]--;
            if (in_d[u] == 0)
                queue[qback++] = u;
            c &= c - 1;
        }
    }

    if (topo_count < N)
        return 0; /* cyclic */

    /* Total order check: every pair must be comparable in transitive closure */
    /* Use Floyd-Warshall on the core */
    uint8_t reach[N];
    memcpy(reach, core, sizeof(core));

    for (int k = 0; k < N; k++)
        for (int i = 0; i < N; i++)
            if (reach[i] & (1 << k))
                reach[i] |= reach[k];

    /* Check all pairs comparable */
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++)
            if (!(reach[i] & (1 << j)) && !(reach[j] & (1 << i)))
                return 2; /* partial order, not total */

    return 1; /* total order */
}

/* Count Hamiltonian paths (only called for sim_H=1 tournaments) */
static int count_ham_paths(uint32_t bits) {
    /* DP with bitmask: dp[mask][last] = # paths visiting vertices in mask, ending at last */
    int dp[1 << N][N];
    memset(dp, 0, sizeof(dp));

    for (int v = 0; v < N; v++)
        dp[1 << v][v] = 1;

    for (int mask = 1; mask < (1 << N); mask++) {
        for (int last = 0; last < N; last++) {
            if (!dp[mask][last]) continue;
            if (!(mask & (1 << last))) continue;
            for (int next = 0; next < N; next++) {
                if (mask & (1 << next)) continue;
                if (has_arc(bits, last, next))
                    dp[mask | (1 << next)][next] += dp[mask][last];
            }
        }
    }

    int full = (1 << N) - 1;
    int total = 0;
    for (int v = 0; v < N; v++)
        total += dp[full][v];
    return total;
}

int main(void) {
    init_arc_indices();
    init_triples();

    time_t start = time(NULL);

    long long total_cyclic = 0;
    long long total_total_order = 0;
    long long total_partial = 0;
    long long total_transitive = 0; /* H=1 */
    long long total_nt65 = 0;       /* H=65 */
    long long total_other = 0;

    int h_counts[256] = {0};

    printf("Simplicial Rédei exhaustive at n=%d (%u tournaments)\n", N, TOTAL_TOURNAMENTS);
    fflush(stdout);

    for (uint32_t bits = 0; bits < TOTAL_TOURNAMENTS; bits++) {
        if ((bits & 0xFFFFFFF) == 0 && bits > 0) {
            double pct = 100.0 * bits / TOTAL_TOURNAMENTS;
            time_t elapsed = time(NULL) - start;
            double rate = (double)bits / (elapsed > 0 ? elapsed : 1);
            double eta = (TOTAL_TOURNAMENTS - bits) / rate;
            printf("  %u/%u (%.1f%%) elapsed=%lds rate=%.0f/s ETA=%.0fs\n",
                   bits, TOTAL_TOURNAMENTS, pct, elapsed, rate, eta);
            fflush(stdout);
        }

        int result = process_tournament(bits);
        if (result == 0) {
            total_cyclic++;
        } else if (result == 1) {
            total_total_order++;
            /* Count H for statistics */
            int H = count_ham_paths(bits);
            if (H < 256) h_counts[H]++;
            if (H == 1) total_transitive++;
            else if (H == 65) total_nt65++;
            else total_other++;
        } else {
            total_partial++;
            if (total_partial <= 5) {
                printf("  *** UNEXPECTED: DAG but not total order at bits=%u ***\n", bits);
                int H = count_ham_paths(bits);
                printf("      H=%d\n", H);
            }
        }
    }

    time_t elapsed = time(NULL) - start;

    printf("\n======================================================================\n");
    printf("RESULTS: Simplicial Rédei at n=%d (EXHAUSTIVE)\n", N);
    printf("======================================================================\n\n");

    printf("sim_H = 0 (cyclic core):     %lld\n", total_cyclic);
    printf("sim_H = 1 (total order):     %lld\n", total_total_order);
    printf("sim_H > 1 (partial DAG):     %lld\n", total_partial);

    if (total_partial > 0) {
        printf("\n  *** SIMPLICIAL REDEI FAILS at n=%d! ***\n", N);
    } else {
        printf("\n  *** SIMPLICIAL REDEI HOLDS: sim_H in {0,1} for ALL %u tournaments! ***\n",
               TOTAL_TOURNAMENTS);
    }

    printf("\nBreakdown of sim_H = 1 (%lld total):\n", total_total_order);
    printf("  Transitive (H=1):    %lld  (expected: %d = %d!)\n",
           total_transitive, 40320, N);
    printf("  H=65 (=2^%d+1):      %lld  (expected: %d = %d!)\n",
           N-2, total_nt65, 40320, N);
    printf("  Other H values:      %lld\n", total_other);
    printf("  Total sim_H=1:       %lld  (expected: %d = 2*%d!)\n",
           total_total_order, 80640, N);

    printf("\nH value distribution for sim_H=1:\n");
    for (int h = 0; h < 256; h++)
        if (h_counts[h] > 0)
            printf("  H=%3d: %d tournaments\n", h, h_counts[h]);

    printf("\nFraction with sim_H=1: %lld/%u = %.8f\n",
           total_total_order, TOTAL_TOURNAMENTS,
           (double)total_total_order / TOTAL_TOURNAMENTS);

    double predicted = 2.0 * 40320.0 / TOTAL_TOURNAMENTS;
    printf("Predicted fraction (2*n!/2^C(n,2)): %.8f\n", predicted);
    printf("\nElapsed: %ld seconds\n", elapsed);

    return 0;
}
