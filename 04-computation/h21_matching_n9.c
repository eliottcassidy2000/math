/*
 * h21_matching_n9.c — At n=9, for cycle-rich tournaments with t3 <= 10,
 * check if 3 pairwise-disjoint 3-cycles always exist (alpha_3 >= 1).
 *
 * If yes: Part C immediately blocks H=21 for ALL n >= 9.
 *
 * Strategy: Sample random n=9 tournaments, filter to cycle-rich + t3<=10,
 * then find maximum matching of disjoint 3-cycle vertex sets.
 *
 * Author: opus-2026-03-07-S42
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N 9

static int adj[N][N];

void random_tournament() {
    for (int i = 0; i < N; i++) {
        adj[i][i] = 0;
        for (int j = i+1; j < N; j++) {
            if (rand() & 1) {
                adj[i][j] = 1; adj[j][i] = 0;
            } else {
                adj[i][j] = 0; adj[j][i] = 1;
            }
        }
    }
}

int has_source_sink() {
    for (int i = 0; i < N; i++) {
        int s = 0;
        for (int j = 0; j < N; j++) s += adj[i][j];
        if (s == 0 || s == N-1) return 1;
    }
    return 0;
}

int count_3cycles() {
    int scores[N];
    for (int i = 0; i < N; i++) {
        scores[i] = 0;
        for (int j = 0; j < N; j++) scores[i] += adj[i][j];
    }
    int sum_c2 = 0;
    for (int i = 0; i < N; i++) sum_c2 += scores[i] * (scores[i]-1) / 2;
    return 84 - sum_c2;  /* C(9,3) = 84 */
}

int vertex_in_3cycle(int v) {
    for (int a = 0; a < N; a++) {
        if (a == v) continue;
        for (int b = a+1; b < N; b++) {
            if (b == v) continue;
            /* Check v->a->b->v or v->b->a->v */
            if (adj[v][a] && adj[a][b] && adj[b][v]) return 1;
            if (adj[v][b] && adj[b][a] && adj[a][v]) return 1;
        }
    }
    return 0;
}

int all_in_3cycle() {
    for (int v = 0; v < N; v++)
        if (!vertex_in_3cycle(v)) return 0;
    return 1;
}

/* Find all 3-cycle vertex sets */
typedef struct { int a, b, c; } Triple;
Triple cycles3[200];
int n_cycles3;

void find_3cycle_sets() {
    n_cycles3 = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            for (int c = b+1; c < N; c++) {
                if ((adj[a][b] && adj[b][c] && adj[c][a]) ||
                    (adj[a][c] && adj[c][b] && adj[b][a])) {
                    if (n_cycles3 < 200) {
                        cycles3[n_cycles3].a = a;
                        cycles3[n_cycles3].b = b;
                        cycles3[n_cycles3].c = c;
                        n_cycles3++;
                    }
                }
            }
}

int disjoint(Triple *x, Triple *y) {
    if (x->a == y->a || x->a == y->b || x->a == y->c) return 0;
    if (x->b == y->a || x->b == y->b || x->b == y->c) return 0;
    if (x->c == y->a || x->c == y->b || x->c == y->c) return 0;
    return 1;
}

/* Find maximum matching of disjoint 3-cycles */
int max_matching() {
    int best = 0;
    for (int i = 0; i < n_cycles3 && best < 3; i++) {
        for (int j = i+1; j < n_cycles3 && best < 3; j++) {
            if (!disjoint(&cycles3[i], &cycles3[j])) continue;
            if (best < 2) best = 2;
            for (int k = j+1; k < n_cycles3; k++) {
                if (disjoint(&cycles3[i], &cycles3[k]) &&
                    disjoint(&cycles3[j], &cycles3[k])) {
                    return 3;  /* Found 3 disjoint! */
                }
            }
        }
    }
    if (best == 0 && n_cycles3 > 0) best = 1;
    return best;
}

/* Also compute alpha_1 = number of odd cycles total */
int count_alpha1() {
    /* Just count 3-cycle vertex sets + 5-cycle vertex sets + 7-cycle + 9-cycle */
    int alpha1 = n_cycles3;

    /* 5-cycles via Held-Karp on each 5-subset */
    for (int a = 0; a < N; a++)
      for (int b = a+1; b < N; b++)
        for (int c = b+1; c < N; c++)
          for (int d = c+1; d < N; d++)
            for (int e = d+1; e < N; e++) {
                int v[5] = {a,b,c,d,e};
                /* Count directed Hamiltonian cycles in this sub-tournament */
                /* dp[mask][last]: number of directed paths from v[0] visiting mask ending at v[last] */
                int dp[32][5];
                memset(dp, 0, sizeof(dp));
                dp[1][0] = 1;
                for (int mask = 1; mask < 32; mask++) {
                    for (int i = 0; i < 5; i++) {
                        if (!(mask & (1<<i))) continue;
                        if (dp[mask][i] == 0) continue;
                        for (int j = 0; j < 5; j++) {
                            if (mask & (1<<j)) continue;
                            if (adj[v[i]][v[j]])
                                dp[mask|(1<<j)][j] += dp[mask][i];
                        }
                    }
                }
                int cyc5 = 0;
                for (int j = 1; j < 5; j++)
                    if (dp[31][j] && adj[v[j]][v[0]])
                        cyc5 += dp[31][j];
                /* Each directed 5-cycle counted 5 times (starting from each vertex in sorted order?
                   No — we fix v[0] as start, so each directed cycle is counted once.
                   But there are 2 directed cycles per undirected: one CW, one CCW.
                   A "directed cycle" in a tournament on 5 vertices: either the 5 arcs form a cycle
                   or they don't. We count directed Hamiltonian cycles = paths returning to start.
                   Each undirected 5-cycle gives exactly 1 directed cycle (not 2 — in a tournament,
                   the arcs are fixed). Wait, for a 5-vertex sub-tournament, the number of
                   Hamiltonian circuits is 0, 1, or 2 (one for each direction around the cycle),
                   but in a tournament only one direction can be a cycle.
                   Actually no: a tournament on 5 vertices can have multiple Hamiltonian cycles
                   (up to 24 = 4!). We fix start = v[0], so we count directed H-cycles / 1.

                   The number of VERTEX SETS that support a 5-cycle is what we want for alpha_1.
                   A vertex set supports a 5-cycle iff the sub-tournament has at least one
                   Hamiltonian cycle. So: alpha_1 += (cyc5 > 0 ? 1 : 0). */
                if (cyc5 > 0) alpha1++;
            }

    /* 7-cycles and 9-cycle: skip for now, they only increase alpha_1 */
    return alpha1;
}

int main() {
    srand(time(NULL));
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    long long total = 0, cycle_rich = 0;
    long long match_counts[4] = {0};  /* match_counts[k] = how many have max matching = k */
    long long match2_low_alpha = 0;
    int match2_examples = 0;

    long long TRIALS = 10000000;

    for (long long trial = 0; trial < TRIALS; trial++) {
        random_tournament();

        if (has_source_sink()) continue;
        int t3 = count_3cycles();
        if (t3 > 10) continue;
        if (!all_in_3cycle()) continue;

        total++;
        cycle_rich++;

        find_3cycle_sets();
        int mm = max_matching();
        if (mm < 4) match_counts[mm]++;

        if (mm <= 2) {
            /* This is the interesting case — check alpha_1 */
            int a1 = count_alpha1();
            if (a1 <= 10) {
                match2_low_alpha++;
                if (match2_examples < 5) {
                    match2_examples++;
                    printf("INTERESTING: mm=%d, t3=%d, alpha1>=%d, n_cycles3=%d\n",
                           mm, t3, a1, n_cycles3);
                    printf("Adjacency:\n");
                    for (int i = 0; i < N; i++) {
                        for (int j = 0; j < N; j++) printf("%d", adj[i][j]);
                        printf("\n");
                    }
                    printf("\n");
                }
            }
        }

        if ((trial & 0xFFFFF) == 0xFFFFF) {
            fprintf(stderr, "Trial %lld: cycle_rich=%lld, mm[0]=%lld mm[1]=%lld mm[2]=%lld mm[3]=%lld low_alpha=%lld\n",
                    trial, cycle_rich, match_counts[0], match_counts[1], match_counts[2], match_counts[3], match2_low_alpha);
        }
    }

    printf("\n=== RESULTS (n=9, %lld trials) ===\n", TRIALS);
    printf("Cycle-rich (no src/sink, t3<=10, all in 3cyc): %lld\n", cycle_rich);
    printf("Max matching of disjoint 3-cycles:\n");
    for (int k = 0; k <= 3; k++)
        printf("  mm=%d: %lld (%.4f%%)\n", k, match_counts[k],
               cycle_rich > 0 ? 100.0*match_counts[k]/cycle_rich : 0);
    printf("mm<=2 with alpha1<=10: %lld\n", match2_low_alpha);

    return 0;
}
