/*
 * simplicial_redei_n7.c — Exhaustive simplicial Rédei verification at n=7
 * opus-2026-03-15-S90
 *
 * For each of 2^21 = 2,097,152 tournaments on 7 vertices:
 *   1. Compute transitive core (edges in at least one transitive triple)
 *   2. Check if transitive core is a DAG
 *   3. If DAG, check if it's a total order (sim_H = 1) or partial order (sim_H > 1)
 *   4. If cyclic, sim_H = 0
 *
 * Also counts Hamiltonian paths H(T) for sim_H=1 tournaments to check
 * the pattern H = 2^{n-2} + 1 = 33 for non-transitive sim_H=1 class.
 *
 * Predictions:
 *   - sim_H ∈ {0,1} for all tournaments (simplicial Rédei)
 *   - 2·7! = 10080 tournaments have sim_H = 1
 *   - 5040 transitive (H=1) + 5040 with H=33
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N 7
#define NPAIRS (N*(N-1)/2)  /* 21 */
#define NTOURNAMENTS (1 << NPAIRS)  /* 2097152 */

/* Tournament adjacency matrix stored as bitmask per row */
static int T[N][N];

/* Transitive core: core[i][j] = 1 if edge i→j is in some transitive triple */
static int core[N][N];

/* Build tournament from bit pattern */
static void build_tournament(int bits) {
    int idx = 0;
    memset(T, 0, sizeof(T));
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            if ((bits >> idx) & 1) {
                T[i][j] = 1;
            } else {
                T[j][i] = 1;
            }
            idx++;
        }
    }
}

/* Compute transitive core */
static void compute_core(void) {
    memset(core, 0, sizeof(core));
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            for (int k = j+1; k < N; k++) {
                /* Check all 6 permutations of (i,j,k) for transitive triple */
                int verts[3] = {i, j, k};
                /* Try all orderings a→b→c with a→c */
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        if (b == a) continue;
                        int c = 3 - a - b;
                        int va = verts[a], vb = verts[b], vc = verts[c];
                        if (T[va][vb] && T[vb][vc] && T[va][vc]) {
                            /* Transitive triple found: va > vb > vc */
                            core[va][vb] = 1;
                            core[vb][vc] = 1;
                            core[va][vc] = 1;
                            goto next_triple;
                        }
                    }
                }
                next_triple:;
            }
        }
    }
}

/* Check if core is a DAG using DFS cycle detection */
/* Returns: 0 = has cycle, 1 = DAG */
static int color[N]; /* 0=white, 1=gray, 2=black */

static int dfs_acyclic(int v) {
    color[v] = 1; /* gray */
    for (int u = 0; u < N; u++) {
        if (core[v][u]) {
            if (color[u] == 1) return 0; /* back edge = cycle */
            if (color[u] == 0 && !dfs_acyclic(u)) return 0;
        }
    }
    color[v] = 2; /* black */
    return 1;
}

static int is_dag(void) {
    memset(color, 0, sizeof(color));
    for (int v = 0; v < N; v++) {
        if (color[v] == 0) {
            if (!dfs_acyclic(v)) return 0;
        }
    }
    return 1;
}

/* Check if DAG is a total order (every pair comparable) */
/* Uses transitive closure via Floyd-Warshall */
static int is_total_order(void) {
    int reach[N][N];
    memcpy(reach, core, sizeof(core));

    /* Transitive closure */
    for (int k = 0; k < N; k++)
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                if (reach[i][k] && reach[k][j])
                    reach[i][j] = 1;

    /* Check every pair is comparable */
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++)
            if (!reach[i][j] && !reach[j][i])
                return 0; /* incomparable pair */

    return 1;
}

/* Count linear extensions of the core (brute force, only called if DAG but not total order) */
static int count_linear_extensions(void) {
    /* Generate all permutations and check against core */
    int perm[N], used[N], count = 0;
    memset(used, 0, sizeof(used));

    /* Iterative permutation generator with pruning */
    int pos = 0;
    perm[0] = -1;

    while (pos >= 0) {
        perm[pos]++;
        if (perm[pos] >= N) {
            if (pos > 0) used[perm[pos-1]] = 0;
            pos--;
            continue;
        }
        if (used[perm[pos]]) continue;

        /* Check: perm[pos] must not have a core predecessor placed after it */
        int ok = 1;
        for (int earlier = 0; earlier < pos; earlier++) {
            if (core[perm[pos]][perm[earlier]]) {
                /* perm[pos] should come before perm[earlier], but it's after */
                ok = 0;
                break;
            }
        }
        if (!ok) continue;

        if (pos == N-1) {
            count++;
        } else {
            used[perm[pos]] = 1;
            pos++;
            perm[pos] = -1;
        }
    }
    return count;
}

/* Count Hamiltonian paths */
static int count_ham_paths(void) {
    int perm[N], used[N], count = 0;
    memset(used, 0, sizeof(used));
    int pos = 0;
    perm[0] = -1;

    while (pos >= 0) {
        perm[pos]++;
        if (perm[pos] >= N) {
            if (pos > 0) used[perm[pos-1]] = 0;
            pos--;
            continue;
        }
        if (used[perm[pos]]) continue;

        /* Check arc from previous vertex */
        if (pos > 0 && !T[perm[pos-1]][perm[pos]]) continue;

        if (pos == N-1) {
            count++;
        } else {
            used[perm[pos]] = 1;
            pos++;
            perm[pos] = -1;
        }
    }
    return count;
}

/* Count 3-cycles */
static int count_3cycles(void) {
    int c = 0;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++)
            for (int k = j+1; k < N; k++) {
                int s = T[i][j] + T[j][k] + T[k][i];
                if (s == 0 || s == 3) c++;
            }
    return c;
}

int main(void) {
    time_t start = time(NULL);

    long long total_cyclic = 0;       /* sim_H = 0 */
    long long total_total_order = 0;  /* sim_H = 1 (DAG = total order) */
    long long total_partial = 0;      /* sim_H > 1 (DAG, not total order) — unexpected! */
    long long total_transitive = 0;   /* H = 1 (transitive tournaments) */
    long long total_sim1_H33 = 0;     /* sim_H=1, H=33 */
    long long total_sim1_other = 0;   /* sim_H=1, H != 1 and H != 33 */

    /* H value distribution for sim_H = 1 */
    int h_counts[256];
    memset(h_counts, 0, sizeof(h_counts));

    /* Score sequence tracking for sim_H = 1 non-transitive */
    int max_le = 0;

    printf("Simplicial Rédei exhaustive verification at n=%d\n", N);
    printf("Total tournaments: %d\n\n", NTOURNAMENTS);

    for (int bits = 0; bits < NTOURNAMENTS; bits++) {
        if (bits % 500000 == 0 && bits > 0) {
            printf("  %d/%d processed (%.1f%%) ...\n", bits, NTOURNAMENTS,
                   100.0 * bits / NTOURNAMENTS);
            fflush(stdout);
        }

        build_tournament(bits);
        compute_core();

        if (!is_dag()) {
            total_cyclic++;
            continue;
        }

        if (is_total_order()) {
            total_total_order++;
            int H = count_ham_paths();
            if (H < 256) h_counts[H]++;
            if (H == 1) {
                total_transitive++;
            } else if (H == 33) {
                total_sim1_H33++;
            } else {
                total_sim1_other++;
            }
        } else {
            /* DAG but not total order — count linear extensions */
            int le = count_linear_extensions();
            total_partial++;
            if (le > max_le) max_le = le;
            printf("  UNEXPECTED: DAG but not total order! bits=%d, LE=%d\n", bits, le);
            if (total_partial <= 5) {
                int H = count_ham_paths();
                int c3 = count_3cycles();
                printf("    H=%d, c3=%d\n", H, c3);
                /* Print score sequence */
                printf("    scores=(");
                for (int v = 0; v < N; v++) {
                    int s = 0;
                    for (int u = 0; u < N; u++) s += T[v][u];
                    printf("%d%s", s, v < N-1 ? "," : "");
                }
                printf(")\n");
            }
        }
    }

    time_t elapsed = time(NULL) - start;

    printf("\n======================================================================\n");
    printf("RESULTS: Simplicial Rédei at n=%d\n", N);
    printf("======================================================================\n\n");

    printf("sim_H = 0 (cyclic core):     %lld\n", total_cyclic);
    printf("sim_H = 1 (total order):     %lld\n", total_total_order);
    printf("sim_H > 1 (partial DAG):     %lld\n", total_partial);
    if (total_partial > 0) {
        printf("  *** SIMPLICIAL REDEI FAILS! max LE = %d ***\n", max_le);
    } else {
        printf("\n  *** SIMPLICIAL REDEI HOLDS: sim_H in {0,1} for ALL %d tournaments! ***\n",
               NTOURNAMENTS);
    }

    printf("\nBreakdown of sim_H = 1 (%lld total):\n", total_total_order);
    printf("  Transitive (H=1):    %lld  (expected: %d = %d!)\n",
           total_transitive, 5040, N);
    printf("  H=33 (=2^%d+1):      %lld  (expected: %d = %d!)\n",
           N-2, total_sim1_H33, 5040, N);
    printf("  Other H values:      %lld\n", total_sim1_other);
    printf("  Total sim_H=1:       %lld  (expected: %d = 2*%d!)\n",
           total_total_order, 10080, N);

    printf("\nH value distribution for sim_H=1:\n");
    for (int h = 0; h < 256; h++) {
        if (h_counts[h] > 0) {
            printf("  H=%3d: %d tournaments\n", h, h_counts[h]);
        }
    }

    printf("\nFraction with sim_H=1: %lld/%d = %.6f\n",
           total_total_order, NTOURNAMENTS,
           (double)total_total_order / NTOURNAMENTS);

    /* Pattern check: 2*n!/2^C(n,2) */
    double predicted = 2.0 * 5040.0 / NTOURNAMENTS;
    printf("Predicted fraction (2*n!/2^C(n,2)): %.6f\n", predicted);

    printf("\nElapsed: %ld seconds\n", elapsed);

    return 0;
}
