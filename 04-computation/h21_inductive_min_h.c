/*
 * h21_inductive_min_h.c — For cycle-rich n=9 tournaments, check:
 * Does every cycle-rich n=9 tournament contain a vertex v such that
 * T-v is cycle-rich at n=8?
 *
 * If yes: H(T) >= H(T-v) (by I(G,2) monotonicity under subgraph inclusion?
 *         Actually, Omega(T-v) is a subgraph of Omega(T), so I grows.)
 *         And H(T-v) >= 25 (exhaustive at n=8).
 *         But 25 > 21, so H(T) >= 25 > 21.
 *
 * Wait, Omega(T-v) is NOT necessarily a subgraph of Omega(T).
 * Removing v removes all cycles through v. But cycles NOT through v
 * are the SAME in T and T-v. So Omega(T-v) = Omega(T) restricted
 * to cycles not through v. That IS an induced subgraph.
 *
 * And I(G,2) >= I(G-S,2) for any S (subset deletion increases
 * independent sets... wait, deleting vertices from G DECREASES
 * independent sets? No — removing vertex u from G means removing
 * all independent sets that include u. So I(G,2) >= I(G-u,2).
 *
 * Therefore: H(T) = I(Omega(T),2) >= I(Omega(T-v),2) = H(T-v).
 *
 * So if T-v is cycle-rich at n=8, then H(T) >= H(T-v) >= 25 > 21.
 *
 * KEY QUESTION: Does every cycle-rich n=9 tournament have a vertex
 * whose deletion gives a cycle-rich tournament?
 *
 * "Cycle-rich" means: no source/sink, every vertex in 3-cycle.
 * After deleting v from T:
 * - Source/sink: vertex w could become source/sink in T-v if all its
 *   in/out neighbors (other than v) go one way.
 * - Every vertex in 3-cycle: some vertex might lose all its 3-cycles
 *   (if all its 3-cycles went through v).
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

int has_source_sink_n(int n, int skip) {
    for (int i = 0; i < N; i++) {
        if (i == skip) continue;
        int s = 0;
        for (int j = 0; j < N; j++) {
            if (j == skip) continue;
            s += adj[i][j];
        }
        if (s == 0 || s == n-1) return 1;
    }
    return 0;
}

int vertex_in_3cycle_skip(int v, int skip) {
    for (int a = 0; a < N; a++) {
        if (a == v || a == skip) continue;
        for (int b = a+1; b < N; b++) {
            if (b == v || b == skip) continue;
            if (adj[v][a] && adj[a][b] && adj[b][v]) return 1;
            if (adj[v][b] && adj[b][a] && adj[a][v]) return 1;
        }
    }
    return 0;
}

int all_in_3cycle_skip(int skip) {
    for (int v = 0; v < N; v++) {
        if (v == skip) continue;
        if (!vertex_in_3cycle_skip(v, skip)) return 0;
    }
    return 1;
}

int count_3cycles_n9() {
    int scores[N];
    for (int i = 0; i < N; i++) {
        scores[i] = 0;
        for (int j = 0; j < N; j++) scores[i] += adj[i][j];
    }
    int sum_c2 = 0;
    for (int i = 0; i < N; i++) sum_c2 += scores[i] * (scores[i]-1) / 2;
    return 84 - sum_c2;
}

int main() {
    srand(time(NULL));
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    long long TRIALS = 10000000;
    long long cycle_rich = 0;
    long long has_good_deletion = 0;
    long long no_good_deletion = 0;

    for (long long trial = 0; trial < TRIALS; trial++) {
        random_tournament();
        if (has_source_sink_n(N, -1)) continue;
        int t3 = count_3cycles_n9();
        if (t3 > 10) continue;

        /* Check all vertices in 3-cycle */
        int all_cyc = 1;
        for (int v = 0; v < N; v++)
            if (!vertex_in_3cycle_skip(v, -1)) { all_cyc = 0; break; }
        if (!all_cyc) continue;

        cycle_rich++;

        /* Try each vertex deletion */
        int found_good = 0;
        for (int v = 0; v < N; v++) {
            /* Check if T-v is cycle-rich at n=8 */
            if (has_source_sink_n(N-1, v)) continue;
            if (!all_in_3cycle_skip(v)) continue;
            /* Also check t3 <= 10 in T-v */
            /* T-v has C(8,3) = 56 triples */
            int scores_tv[N];
            for (int i = 0; i < N; i++) {
                if (i == v) continue;
                scores_tv[i] = 0;
                for (int j = 0; j < N; j++) {
                    if (j == v) continue;
                    scores_tv[i] += adj[i][j];
                }
            }
            int sum_c2 = 0;
            for (int i = 0; i < N; i++) {
                if (i == v) continue;
                sum_c2 += scores_tv[i] * (scores_tv[i]-1) / 2;
            }
            int t3_tv = 56 - sum_c2;
            if (t3_tv <= 10) {
                found_good = 1;
                break;
            }
        }

        if (found_good) {
            has_good_deletion++;
        } else {
            no_good_deletion++;
            if (no_good_deletion <= 5) {
                printf("NO GOOD DELETION! t3=%d\n", t3);
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) printf("%d", adj[i][j]);
                    printf("\n");
                }
                printf("\n");
            }
        }

        if ((trial & 0xFFFFF) == 0xFFFFF) {
            fprintf(stderr, "Trial %lld: cycle_rich=%lld good=%lld bad=%lld\n",
                    trial, cycle_rich, has_good_deletion, no_good_deletion);
        }
    }

    printf("\n=== RESULTS (n=%d) ===\n", N);
    printf("Cycle-rich: %lld\n", cycle_rich);
    printf("Has good deletion: %lld (%.4f%%)\n", has_good_deletion,
           cycle_rich > 0 ? 100.0*has_good_deletion/cycle_rich : 0);
    printf("No good deletion: %lld (%.4f%%)\n", no_good_deletion,
           cycle_rich > 0 ? 100.0*no_good_deletion/cycle_rich : 0);

    return 0;
}
