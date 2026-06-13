/*
 * h21_alpha10_h_values_n9.c — At n=9, for no-src-sink tournaments with
 * alpha_1 <= 12, compute full H and i_2 to understand why H != 21.
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

typedef unsigned short CycleMask;
CycleMask all_cycles[500];
int n_all_cycles;

void add_cycle(CycleMask m) {
    for (int i = 0; i < n_all_cycles; i++)
        if (all_cycles[i] == m) return;
    if (n_all_cycles < 500)
        all_cycles[n_all_cycles++] = m;
}

int count_3cycle_sets() {
    int count = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            for (int c = b+1; c < N; c++) {
                if ((adj[a][b] && adj[b][c] && adj[c][a]) ||
                    (adj[a][c] && adj[c][b] && adj[b][a])) {
                    add_cycle((1<<a)|(1<<b)|(1<<c));
                    count++;
                }
            }
    return count;
}

int count_5cycle_sets() {
    int count = 0;
    for (int a = 0; a < N; a++)
      for (int b = a+1; b < N; b++)
        for (int c = b+1; c < N; c++)
          for (int d = c+1; d < N; d++)
            for (int e = d+1; e < N; e++) {
                int v[5] = {a,b,c,d,e};
                int dp[32][5];
                memset(dp, 0, sizeof(dp));
                dp[1][0] = 1;
                for (int mask = 1; mask < 32; mask++)
                    for (int i = 0; i < 5; i++) {
                        if (!(mask & (1<<i)) || dp[mask][i] == 0) continue;
                        for (int j = 0; j < 5; j++) {
                            if (mask & (1<<j)) continue;
                            if (adj[v[i]][v[j]])
                                dp[mask|(1<<j)][j] += dp[mask][i];
                        }
                    }
                int cyc = 0;
                for (int j = 1; j < 5; j++)
                    if (dp[31][j] && adj[v[j]][v[0]]) cyc += dp[31][j];
                if (cyc > 0) {
                    add_cycle((1<<a)|(1<<b)|(1<<c)|(1<<d)|(1<<e));
                    count++;
                }
            }
    return count;
}

int count_7cycle_sets() {
    int count = 0;
    for (int a = 0; a < N; a++)
      for (int b = a+1; b < N; b++)
        for (int c = b+1; c < N; c++)
          for (int d = c+1; d < N; d++)
            for (int e = d+1; e < N; e++)
              for (int f = e+1; f < N; f++)
                for (int g = f+1; g < N; g++) {
                    int v[7] = {a,b,c,d,e,f,g};
                    int dp[128][7];
                    memset(dp, 0, sizeof(dp));
                    dp[1][0] = 1;
                    for (int mask = 1; mask < 128; mask++)
                        for (int i = 0; i < 7; i++) {
                            if (!(mask & (1<<i)) || dp[mask][i] == 0) continue;
                            for (int j = 0; j < 7; j++) {
                                if (mask & (1<<j)) continue;
                                if (adj[v[i]][v[j]])
                                    dp[mask|(1<<j)][j] += dp[mask][i];
                            }
                        }
                    int cyc = 0;
                    for (int j = 1; j < 7; j++)
                        if (dp[127][j] && adj[v[j]][v[0]]) cyc += dp[127][j];
                    if (cyc > 0) {
                        CycleMask m = 0;
                        for (int x = 0; x < 7; x++) m |= (1 << v[x]);
                        add_cycle(m);
                        count++;
                    }
                }
    return count;
}

int count_9cycle() {
    int dp[512][9];
    memset(dp, 0, sizeof(dp));
    dp[1][0] = 1;
    for (int mask = 1; mask < 512; mask++)
        for (int i = 0; i < 9; i++) {
            if (!(mask & (1<<i)) || dp[mask][i] == 0) continue;
            for (int j = 0; j < 9; j++) {
                if (mask & (1<<j)) continue;
                if (adj[i][j])
                    dp[mask|(1<<j)][j] += dp[mask][i];
            }
        }
    int cyc = 0;
    for (int j = 1; j < 9; j++)
        if (dp[511][j] && adj[j][0]) cyc += dp[511][j];
    if (cyc > 0) { add_cycle(511); return 1; }
    return 0;
}

int count_disjoint_pairs() {
    int i2 = 0;
    for (int i = 0; i < n_all_cycles; i++)
        for (int j = i+1; j < n_all_cycles; j++)
            if ((all_cycles[i] & all_cycles[j]) == 0) i2++;
    return i2;
}

int count_disjoint_triples() {
    int i3 = 0;
    for (int i = 0; i < n_all_cycles; i++)
        for (int j = i+1; j < n_all_cycles; j++) {
            if ((all_cycles[i] & all_cycles[j]) != 0) continue;
            for (int k = j+1; k < n_all_cycles; k++)
                if ((all_cycles[i] & all_cycles[k]) == 0 &&
                    (all_cycles[j] & all_cycles[k]) == 0)
                    i3++;
        }
    return i3;
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

    long long TRIALS = 10000000;
    long long no_srcsink = 0;
    int found = 0;

    for (long long trial = 0; trial < TRIALS; trial++) {
        random_tournament();
        if (has_source_sink()) continue;
        no_srcsink++;

        n_all_cycles = 0;
        int t3 = count_3cycle_sets();
        if (t3 > 12) continue;

        int t5 = count_5cycle_sets();
        int alpha1 = t3 + t5;
        if (alpha1 > 12) continue;

        int t7 = count_7cycle_sets();
        alpha1 += t7;
        if (alpha1 > 12) continue;

        int t9 = count_9cycle();
        alpha1 += t9;
        if (alpha1 > 12) continue;

        /* Found a low-alpha1 tournament! Compute everything */
        int i2 = count_disjoint_pairs();
        int i3 = count_disjoint_triples();
        long long H = hamiltonian_paths();

        int scores[N];
        for (int i = 0; i < N; i++) {
            scores[i] = 0;
            for (int j = 0; j < N; j++) scores[i] += adj[i][j];
        }

        found++;
        printf("a1=%d (t3=%d t5=%d t7=%d t9=%d) i2=%d i3=%d H=%lld scores=(",
               alpha1, t3, t5, t7, t9, i2, i3, H);
        for (int i = 0; i < N; i++) printf("%d%s", scores[i], i<N-1?",":")\n");

        if ((trial & 0x3FFFF) == 0x3FFFF) {
            fprintf(stderr, "Trial %lld: no_srcsink=%lld found=%d\n",
                    trial, no_srcsink, found);
        }
    }

    printf("\n=== RESULTS (n=%d, %lld trials) ===\n", N, TRIALS);
    printf("No source/sink: %lld\n", no_srcsink);
    printf("alpha1<=12 found: %d\n", found);

    return 0;
}
