/*
 * h21_alpha10_structure_n9.c — At n=9, check if alpha_1=10 always implies
 * a source or sink exists.
 *
 * For H=21 impossibility: the (10,0) decomposition requires alpha_1=10, i_2=0.
 * At n=7,8: every tournament with alpha_1=10 has a source or sink.
 * Does this hold at n=9?
 *
 * Strategy: Sample random n=9 tournaments. For those with alpha_1 near 10
 * and no source/sink, output details.
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
    return 84 - sum_c2;
}

/* Count 3-cycle vertex SETS (not directed) */
int count_3cycle_sets() {
    int count = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            for (int c = b+1; c < N; c++) {
                if ((adj[a][b] && adj[b][c] && adj[c][a]) ||
                    (adj[a][c] && adj[c][b] && adj[b][a]))
                    count++;
            }
    return count;
}

/* Count 5-cycle vertex sets */
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
                int cyc = 0;
                for (int j = 1; j < 5; j++)
                    if (dp[31][j] && adj[v[j]][v[0]]) cyc += dp[31][j];
                if (cyc > 0) count++;
            }
    return count;
}

/* Count 7-cycle vertex sets */
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
                    for (int mask = 1; mask < 128; mask++) {
                        for (int i = 0; i < 7; i++) {
                            if (!(mask & (1<<i))) continue;
                            if (dp[mask][i] == 0) continue;
                            for (int j = 0; j < 7; j++) {
                                if (mask & (1<<j)) continue;
                                if (adj[v[i]][v[j]])
                                    dp[mask|(1<<j)][j] += dp[mask][i];
                            }
                        }
                    }
                    int cyc = 0;
                    for (int j = 1; j < 7; j++)
                        if (dp[127][j] && adj[v[j]][v[0]]) cyc += dp[127][j];
                    if (cyc > 0) count++;
                }
    return count;
}

/* Count 9-cycle (full Hamiltonian cycle) */
int count_9cycle() {
    int dp[512][9];
    memset(dp, 0, sizeof(dp));
    dp[1][0] = 1;
    for (int mask = 1; mask < 512; mask++) {
        for (int i = 0; i < 9; i++) {
            if (!(mask & (1<<i))) continue;
            if (dp[mask][i] == 0) continue;
            for (int j = 0; j < 9; j++) {
                if (mask & (1<<j)) continue;
                if (adj[i][j])
                    dp[mask|(1<<j)][j] += dp[mask][i];
            }
        }
    }
    int cyc = 0;
    for (int j = 1; j < 9; j++)
        if (dp[511][j] && adj[j][0]) cyc += dp[511][j];
    return (cyc > 0) ? 1 : 0;
}

int main() {
    srand(time(NULL));
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    long long trials_done = 0;
    int found_alpha10_no_srcsink = 0;
    /* Distribution of alpha_1 for no-src-sink tournaments */
    long long alpha1_hist[300] = {0};
    long long no_srcsink_count = 0;
    long long TRIALS = 5000000;

    for (long long trial = 0; trial < TRIALS; trial++) {
        random_tournament();
        trials_done++;

        if (has_source_sink()) continue;
        no_srcsink_count++;

        int t3 = count_3cycle_sets();
        int t5 = count_5cycle_sets();
        int alpha1 = t3 + t5;

        /* Quick check: if alpha1 already > 20, skip 7/9-cycle counting */
        if (alpha1 <= 15) {
            int t7 = count_7cycle_sets();
            alpha1 += t7;
            if (alpha1 <= 15) {
                int t9 = count_9cycle();
                alpha1 += t9;
            }
        }

        if (alpha1 < 300) alpha1_hist[alpha1]++;

        if (alpha1 <= 10) {
            found_alpha10_no_srcsink++;
            printf("FOUND: alpha1=%d (t3=%d, t5=%d) no src/sink\n", alpha1, t3, t5);
            int scores[N];
            for (int i = 0; i < N; i++) {
                scores[i] = 0;
                for (int j = 0; j < N; j++) scores[i] += adj[i][j];
            }
            printf("Scores: ");
            for (int i = 0; i < N; i++) printf("%d ", scores[i]);
            printf("\n");
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) printf("%d", adj[i][j]);
                printf("\n");
            }
            printf("\n");
        }

        if ((trial & 0x3FFFF) == 0x3FFFF) {
            fprintf(stderr, "Trial %lld: no_srcsink=%lld, alpha1<=10=%d\n",
                    trial, no_srcsink_count, found_alpha10_no_srcsink);
        }
    }

    printf("\n=== RESULTS (n=%d, %lld trials) ===\n", N, trials_done);
    printf("No source/sink: %lld\n", no_srcsink_count);
    printf("alpha1<=10 with no src/sink: %d\n", found_alpha10_no_srcsink);
    printf("\nalpha1 distribution (no src/sink, first 50):\n");
    for (int k = 0; k < 50; k++) {
        if (alpha1_hist[k] > 0)
            printf("  alpha1=%d: %lld (%.4f%%)\n", k, alpha1_hist[k],
                   100.0*alpha1_hist[k]/no_srcsink_count);
    }

    return 0;
}
