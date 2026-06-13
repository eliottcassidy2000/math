/*
 * h21_no_deletion_no_matching.c — At n=9, find cycle-rich tournaments
 * that have BOTH:
 * 1. No vertex deletion giving cycle-rich n=8
 * 2. No 3 pairwise-disjoint 3-cycles
 *
 * If no such tournament exists, the proof is complete!
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

int has_source_sink_skip(int skip) {
    for (int i = 0; i < N; i++) {
        if (i == skip) continue;
        int s = 0;
        for (int j = 0; j < N; j++) {
            if (j == skip) continue;
            s += adj[i][j];
        }
        int n_eff = (skip >= 0) ? N-1 : N;
        if (s == 0 || s == n_eff - 1) return 1;
    }
    return 0;
}

int vertex_in_3cycle(int v, int skip) {
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

int all_in_3cycle(int skip) {
    for (int v = 0; v < N; v++) {
        if (v == skip) continue;
        if (!vertex_in_3cycle(v, skip)) return 0;
    }
    return 1;
}

typedef struct { int a, b, c; } Triple;
Triple cycles3[100];
int n_cycles3;

void find_3cycle_sets() {
    n_cycles3 = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            for (int c = b+1; c < N; c++)
                if ((adj[a][b] && adj[b][c] && adj[c][a]) ||
                    (adj[a][c] && adj[c][b] && adj[b][a])) {
                    cycles3[n_cycles3].a = a;
                    cycles3[n_cycles3].b = b;
                    cycles3[n_cycles3].c = c;
                    n_cycles3++;
                }
}

int disjoint(Triple *x, Triple *y) {
    unsigned m1 = (1<<x->a)|(1<<x->b)|(1<<x->c);
    unsigned m2 = (1<<y->a)|(1<<y->b)|(1<<y->c);
    return (m1 & m2) == 0;
}

int has_3_disjoint() {
    for (int i = 0; i < n_cycles3; i++)
        for (int j = i+1; j < n_cycles3; j++) {
            if (!disjoint(&cycles3[i], &cycles3[j])) continue;
            for (int k = j+1; k < n_cycles3; k++)
                if (disjoint(&cycles3[i], &cycles3[k]) &&
                    disjoint(&cycles3[j], &cycles3[k]))
                    return 1;
        }
    return 0;
}

int has_good_deletion() {
    for (int v = 0; v < N; v++) {
        if (has_source_sink_skip(v)) continue;
        if (!all_in_3cycle(v)) continue;
        /* Check t3 <= 10 */
        int scores[N];
        for (int i = 0; i < N; i++) {
            if (i == v) continue;
            scores[i] = 0;
            for (int j = 0; j < N; j++) {
                if (j == v) continue;
                scores[i] += adj[i][j];
            }
        }
        int sum_c2 = 0;
        for (int i = 0; i < N; i++) {
            if (i == v) continue;
            sum_c2 += scores[i] * (scores[i]-1) / 2;
        }
        int t3_tv = 56 - sum_c2;
        if (t3_tv <= 10)
            return 1;  /* This is a valid cycle-rich deletion */
    }
    return 0;
}

long long hamiltonian_paths() {
    long long dp[512][9];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++) dp[1<<v][v] = 1;
    for (int mask = 1; mask < 512; mask++)
        for (int v = 0; v < N; v++) {
            if (!(mask & (1<<v)) || dp[mask][v] == 0) continue;
            for (int u = 0; u < N; u++) {
                if (mask & (1<<u)) continue;
                if (adj[v][u])
                    dp[mask|(1<<u)][u] += dp[mask][v];
            }
        }
    long long total = 0;
    for (int v = 0; v < N; v++) total += dp[511][v];
    return total;
}

int main() {
    srand(time(NULL));
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    long long TRIALS = 50000000;
    long long cycle_rich = 0;
    long long bad_both = 0;

    for (long long trial = 0; trial < TRIALS; trial++) {
        random_tournament();
        if (has_source_sink_skip(-1)) continue;

        /* Check cycle-rich */
        int t3_quick = 0;
        {
            int scores[N];
            for (int i = 0; i < N; i++) {
                scores[i] = 0;
                for (int j = 0; j < N; j++) scores[i] += adj[i][j];
            }
            int sum_c2 = 0;
            for (int i = 0; i < N; i++) sum_c2 += scores[i]*(scores[i]-1)/2;
            t3_quick = 84 - sum_c2;
        }
        if (t3_quick > 10) continue;
        if (!all_in_3cycle(-1)) continue;

        cycle_rich++;

        find_3cycle_sets();
        int has_3d = has_3_disjoint();
        int has_gd = has_good_deletion();

        if (!has_3d && !has_gd) {
            bad_both++;
            long long H = hamiltonian_paths();
            printf("BAD BOTH! t3=%d n_cyc=%d H=%lld\n", n_cycles3, n_cycles3, H);
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) printf("%d", adj[i][j]);
                printf("\n");
            }
            printf("3-cycles: ");
            for (int i = 0; i < n_cycles3; i++)
                printf("{%d,%d,%d} ", cycles3[i].a, cycles3[i].b, cycles3[i].c);
            printf("\n\n");
        }

        if ((trial & 0xFFFFF) == 0xFFFFF) {
            fprintf(stderr, "Trial %lld: cr=%lld bad=%lld\n",
                    trial, cycle_rich, bad_both);
        }
    }

    printf("\n=== RESULTS ===\n");
    printf("Cycle-rich: %lld\n", cycle_rich);
    printf("No 3-disjoint AND no good deletion: %lld\n", bad_both);

    return 0;
}
