/*
 * h21_10_0_n9_check.c — At n=9, check if (10,0) decomposition is achievable:
 * alpha_1=10, all pairwise-sharing (i_2=0), no source/sink.
 *
 * Strategy: Sample random n=9 tournaments. For each:
 * 1. Skip source/sink
 * 2. Count all 3-cycle vertex sets
 * 3. If t3 <= 10, count 5-cycle vertex sets
 * 4. If t3+t5 <= 10, count 7-cycle vertex sets + 9-cycle
 * 5. If total alpha_1 = 10, check if all pairwise-sharing
 *
 * Also track: min H for tournaments with alpha_1 around 10
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

/* Store cycles as vertex bitmasks for fast disjointness check */
typedef unsigned short CycleMask;
CycleMask all_cycles[500];
int n_all_cycles;

/* Add cycle if not already present (by vertex set) */
void add_cycle(CycleMask m) {
    for (int i = 0; i < n_all_cycles; i++)
        if (all_cycles[i] == m) return;  /* already have this vertex set */
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
    if (cyc > 0) {
        add_cycle(511);
        return 1;
    }
    return 0;
}

/* Count disjoint pairs */
int count_disjoint_pairs() {
    int i2 = 0;
    for (int i = 0; i < n_all_cycles; i++)
        for (int j = i+1; j < n_all_cycles; j++)
            if ((all_cycles[i] & all_cycles[j]) == 0)
                i2++;
    return i2;
}

/* Held-Karp for full H(T) */
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

    long long TRIALS = 2000000;
    long long no_srcsink = 0;
    long long alpha_10_count = 0;
    long long alpha_10_no_disjoint = 0;  /* (10,0) candidates */
    long long h21_count = 0;
    int min_alpha1_no_srcsink = 999;

    /* Track alpha_1 histogram for no-src-sink */
    long long a1_hist[100] = {0};

    for (long long trial = 0; trial < TRIALS; trial++) {
        random_tournament();
        if (has_source_sink()) continue;
        no_srcsink++;

        n_all_cycles = 0;
        int t3 = count_3cycle_sets();

        if (t3 > 12) {
            if (t3 < 100) a1_hist[t3]++;  /* lower bound */
            continue;  /* alpha_1 >= t3 > 12, can't be 10 exactly (well, t3 is vertex-set count) */
        }

        int t5 = count_5cycle_sets();
        int alpha1 = t3 + t5;

        if (alpha1 > 12) {
            if (alpha1 < 100) a1_hist[alpha1]++;
            continue;
        }

        int t7 = count_7cycle_sets();
        alpha1 += t7;

        if (alpha1 > 12) {
            if (alpha1 < 100) a1_hist[alpha1]++;
            continue;
        }

        int t9 = count_9cycle();
        alpha1 += t9;

        if (alpha1 < 100) a1_hist[alpha1]++;
        if (alpha1 < min_alpha1_no_srcsink) {
            min_alpha1_no_srcsink = alpha1;
            printf("NEW MIN alpha1=%d (t3=%d t5=%d t7=%d t9=%d) no src/sink\n",
                   alpha1, t3, t5, t7, t9);
        }

        if (alpha1 == 10) {
            alpha_10_count++;
            int i2 = count_disjoint_pairs();
            if (i2 == 0) {
                alpha_10_no_disjoint++;
                long long H = hamiltonian_paths();
                printf("(10,0) CANDIDATE: alpha1=10 i2=0 H=%lld\n", H);
                if (H == 21) {
                    h21_count++;
                    printf("!!! H=21 FOUND !!!\n");
                }
            }
        }

        /* Also check H=21 for any alpha1 <= 10 */
        if (alpha1 <= 10) {
            long long H = hamiltonian_paths();
            if (H == 21) {
                h21_count++;
                printf("H=21! alpha1=%d\n", alpha1);
            }
        }

        if ((trial & 0x1FFFF) == 0x1FFFF) {
            fprintf(stderr, "Trial %lld: no_srcsink=%lld min_a1=%d a10=%lld (10,0)=%lld h21=%lld\n",
                    trial, no_srcsink, min_alpha1_no_srcsink, alpha_10_count,
                    alpha_10_no_disjoint, h21_count);
        }
    }

    printf("\n=== RESULTS (n=%d, %lld trials) ===\n", N, TRIALS);
    printf("No source/sink: %lld\n", no_srcsink);
    printf("Min alpha1 (no src/sink): %d\n", min_alpha1_no_srcsink);
    printf("alpha1=10: %lld\n", alpha_10_count);
    printf("(10,0) candidates: %lld\n", alpha_10_no_disjoint);
    printf("H=21 found: %lld\n", h21_count);
    printf("\nalpha1 distribution (first 30):\n");
    for (int k = 0; k < 30; k++)
        if (a1_hist[k]) printf("  a1=%d: %lld\n", k, a1_hist[k]);

    return 0;
}
