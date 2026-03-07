/*
 * h21_match2_h_values.c — At n=9, for cycle-rich tournaments with
 * max 3-cycle matching <= 2, what is the minimum H?
 *
 * If min H > 21, Sub-case 2b is resolved.
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

int vertex_in_3cycle(int v) {
    for (int a = 0; a < N; a++) {
        if (a == v) continue;
        for (int b = a+1; b < N; b++) {
            if (b == v) continue;
            if (adj[v][a] && adj[a][b] && adj[b][v]) return 1;
            if (adj[v][b] && adj[b][a] && adj[a][v]) return 1;
        }
    }
    return 0;
}

typedef struct { int a, b, c; } Triple;
Triple cycles3[200];
int n_cycles3;

void find_3cycle_sets() {
    n_cycles3 = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            for (int c = b+1; c < N; c++)
                if ((adj[a][b] && adj[b][c] && adj[c][a]) ||
                    (adj[a][c] && adj[c][b] && adj[b][a]))
                    cycles3[n_cycles3++] = (Triple){a, b, c};
}

int disjoint(Triple *x, Triple *y) {
    unsigned m1 = (1u<<x->a)|(1u<<x->b)|(1u<<x->c);
    unsigned m2 = (1u<<y->a)|(1u<<y->b)|(1u<<y->c);
    return (m1 & m2) == 0;
}

int max_matching() {
    int best = 0;
    for (int i = 0; i < n_cycles3 && best < 3; i++)
        for (int j = i+1; j < n_cycles3 && best < 3; j++) {
            if (!disjoint(&cycles3[i], &cycles3[j])) continue;
            if (best < 2) best = 2;
            for (int k = j+1; k < n_cycles3; k++)
                if (disjoint(&cycles3[i], &cycles3[k]) &&
                    disjoint(&cycles3[j], &cycles3[k]))
                    return 3;
        }
    if (best == 0 && n_cycles3 > 0) best = 1;
    return best;
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

    long long TRIALS = 20000000;
    long long cycle_rich = 0;
    long long match_le2 = 0;
    long long min_h_m2 = 999999;
    long long h_hist[200] = {0};

    for (long long trial = 0; trial < TRIALS; trial++) {
        random_tournament();
        if (has_source_sink()) continue;
        {
            int scores[N], sum_c2 = 0;
            for (int i = 0; i < N; i++) {
                scores[i] = 0;
                for (int j = 0; j < N; j++) scores[i] += adj[i][j];
            }
            for (int i = 0; i < N; i++) sum_c2 += scores[i]*(scores[i]-1)/2;
            if (84 - sum_c2 > 10) continue;
        }
        int all_cyc = 1;
        for (int v = 0; v < N; v++)
            if (!vertex_in_3cycle(v)) { all_cyc = 0; break; }
        if (!all_cyc) continue;

        cycle_rich++;
        find_3cycle_sets();
        int mm = max_matching();
        if (mm >= 3) continue;

        match_le2++;
        long long H = hamiltonian_paths();
        if (H < min_h_m2) {
            min_h_m2 = H;
            printf("NEW MIN H=%lld (mm=%d, t3=%d)\n", H, mm, n_cycles3);
        }
        if (H < 200) h_hist[H]++;

        if ((trial & 0xFFFFF) == 0xFFFFF) {
            fprintf(stderr, "Trial %lld: cr=%lld m<=2=%lld min_h=%lld\n",
                    trial, cycle_rich, match_le2, min_h_m2);
        }
    }

    printf("\n=== RESULTS ===\n");
    printf("Cycle-rich: %lld\n", cycle_rich);
    printf("Max matching <= 2: %lld\n", match_le2);
    printf("Min H (mm<=2): %lld\n", min_h_m2);
    printf("\nH distribution for mm<=2 (first 100):\n");
    for (int h = 0; h < 100; h++)
        if (h_hist[h]) printf("  H=%d: %lld\n", h, h_hist[h]);

    return 0;
}
