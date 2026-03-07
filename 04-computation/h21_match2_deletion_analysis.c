/*
 * h21_match2_deletion_analysis.c — For mm<=2 cycle-rich n=9 tournaments,
 * analyze which vertices can be deleted to maintain cycle-rich at n=8.
 *
 * Key question: is there always at least one good deletion, OR
 * can we show H > 21 directly via 5-cycle counting?
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

int has_source_sink_skip(int skip) {
    for (int i = 0; i < N; i++) {
        if (i == skip) continue;
        int s = 0;
        for (int j = 0; j < N; j++) {
            if (j == skip) continue;
            s += adj[i][j];
        }
        if (s == 0 || s == N-2) return 1;  /* n-1 = 8-1 = 7 */
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
                if (cyc > 0) count++;
            }
    return count;
}

int main() {
    srand(time(NULL));
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    long long TRIALS = 20000000;
    long long cycle_rich = 0;
    long long match_le2 = 0;
    int examples = 0;

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
        if (!all_in_3cycle_skip(-1)) continue;

        cycle_rich++;
        find_3cycle_sets();
        int mm = max_matching();
        if (mm >= 3) continue;

        match_le2++;
        if (examples >= 20) continue;

        int t5 = count_5cycle_sets();
        int scores[N];
        for (int i = 0; i < N; i++) {
            scores[i] = 0;
            for (int j = 0; j < N; j++) scores[i] += adj[i][j];
        }

        /* Check each deletion */
        examples++;
        printf("=== Example %d: mm=%d t3=%d t5=%d scores=(", examples, mm, n_cycles3, t5);
        for (int i = 0; i < N; i++) printf("%d%s", scores[i], i<N-1?",":")\n");
        printf("3-cycles: ");
        for (int i = 0; i < n_cycles3; i++)
            printf("{%d,%d,%d} ", cycles3[i].a, cycles3[i].b, cycles3[i].c);
        printf("\n");

        for (int v = 0; v < N; v++) {
            int srcsink = has_source_sink_skip(v);
            int all_cyc = all_in_3cycle_skip(v);
            /* Count 3-cycles in T-v */
            int tv_scores[N], sum_c2 = 0;
            for (int i = 0; i < N; i++) {
                if (i == v) continue;
                tv_scores[i] = 0;
                for (int j = 0; j < N; j++) {
                    if (j == v) continue;
                    tv_scores[i] += adj[i][j];
                }
                sum_c2 += tv_scores[i] * (tv_scores[i]-1)/2;
            }
            int t3_tv = 56 - sum_c2;
            printf("  Del v=%d (s=%d): srcsink=%d all_cyc=%d t3=%d %s\n",
                   v, scores[v], srcsink, all_cyc, t3_tv,
                   (!srcsink && all_cyc && t3_tv <= 10) ? "GOOD" : "");
        }
        printf("\n");

        if ((trial & 0xFFFFF) == 0xFFFFF) {
            fprintf(stderr, "Trial %lld: cr=%lld m<=2=%lld ex=%d\n",
                    trial, cycle_rich, match_le2, examples);
        }
    }

    printf("\n=== RESULTS ===\n");
    printf("Cycle-rich: %lld, mm<=2: %lld\n", cycle_rich, match_le2);

    return 0;
}
