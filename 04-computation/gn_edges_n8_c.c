/*
 * gn_edges_n8_c.c — Compute |E(G_8)| using C for speed
 * opus-2026-03-22-S213
 *
 * n=8: 6880 iso classes, 2^28 = 268M tournaments, 28 arc positions
 * Strategy:
 *   Phase 1: Build tournament fingerprints (score + c3 + deletion-H)
 *            using only CLASS REPRESENTATIVES (not all 268M tournaments)
 *   Phase 2: Use nauty-like partition refinement for fast canonical form
 *
 * ACTUALLY: at n=8, enumerating all 2^28 tournaments is infeasible even in C
 * for the full fingerprint approach. Instead, use a SAMPLING approach:
 *   1. Enumerate by score sequence (much fewer than 2^28)
 *   2. Within each score class, sample to find iso classes
 *   3. For each class, compute neighbor classes
 *
 * OR: use the theoretical bounds from Burnside to guide exploration.
 *
 * For now: compute T_8 = 188288 (already done), predict E_8, and
 * verify with partial computation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define N 8
#define M 28  /* C(8,2) */
#define TOTAL (1 << M)  /* 268435456 */

/* Tournament is a 28-bit integer */
typedef uint32_t tournament_t;

/* Adjacency: bit k corresponds to pair (i,j) where i < j */
/* Pair ordering: (0,1),(0,2),(0,3),...,(0,7),(1,2),(1,3),...,(6,7) */

static int pair_i[M], pair_j[M];
static int pair_idx[N][N]; /* pair_idx[i][j] = index of pair (min(i,j), max(i,j)) */

void init_pairs() {
    int k = 0;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++) {
            pair_i[k] = i;
            pair_j[k] = j;
            pair_idx[i][j] = k;
            pair_idx[j][i] = k;
            k++;
        }
}

/* Get adjacency: returns 1 if i->j */
static inline int adj(tournament_t T, int i, int j) {
    if (i < j) return (T >> pair_idx[i][j]) & 1;
    else return 1 - ((T >> pair_idx[j][i]) & 1);
}

/* Score sequence (sorted) */
void score_seq(tournament_t T, int *scores) {
    for (int i = 0; i < N; i++) {
        int s = 0;
        for (int j = 0; j < N; j++)
            if (i != j) s += adj(T, i, j);
        scores[i] = s;
    }
    /* Sort */
    for (int i = 0; i < N-1; i++)
        for (int j = i+1; j < N; j++)
            if (scores[i] > scores[j]) {
                int tmp = scores[i];
                scores[i] = scores[j];
                scores[j] = tmp;
            }
}

/* Count 3-cycles */
int count_c3(tournament_t T) {
    int c3 = 0;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++)
            for (int k = j+1; k < N; k++) {
                if (adj(T,i,j) && adj(T,j,k) && adj(T,k,i)) c3++;
                if (adj(T,i,k) && adj(T,k,j) && adj(T,j,i)) c3++;
            }
    return c3;
}

/* Hamiltonian path count via DP */
int64_t count_hp(tournament_t T) {
    int64_t dp[1 << N][N];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++)
        dp[1 << v][v] = 1;
    for (uint32_t S = 1; S < (1u << N); S++) {
        for (int v = 0; v < N; v++) {
            if (!(S & (1 << v))) continue;
            if (dp[S][v] == 0) continue;
            for (int u = 0; u < N; u++) {
                if (S & (1 << u)) continue;
                if (adj(T, v, u))
                    dp[S | (1 << u)][u] += dp[S][v];
            }
        }
    }
    int64_t total = 0;
    uint32_t full = (1u << N) - 1;
    for (int v = 0; v < N; v++)
        total += dp[full][v];
    return total;
}

/* H of tournament with vertex v removed */
int64_t hp_minus_v(tournament_t T, int v) {
    /* Build sub-tournament on N-1 vertices */
    int verts[N-1];
    int k = 0;
    for (int i = 0; i < N; i++)
        if (i != v) verts[k++] = i;

    int nn = N-1;
    int64_t dp[1 << (N-1)][N-1];
    memset(dp, 0, sizeof(dp));
    for (int i = 0; i < nn; i++)
        dp[1 << i][i] = 1;
    for (uint32_t S = 1; S < (1u << nn); S++) {
        for (int i = 0; i < nn; i++) {
            if (!(S & (1 << i))) continue;
            if (dp[S][i] == 0) continue;
            for (int j = 0; j < nn; j++) {
                if (S & (1 << j)) continue;
                if (adj(T, verts[i], verts[j]))
                    dp[S | (1 << j)][j] += dp[S][i];
            }
        }
    }
    int64_t total = 0;
    uint32_t full = (1u << nn) - 1;
    for (int i = 0; i < nn; i++)
        total += dp[full][i];
    return total;
}

int main() {
    init_pairs();

    printf("G_8 Edge Count Computation\n");
    printf("n=%d, m=%d, total=%u\n\n", N, M, TOTAL);

    /* Phase 1: Sample random tournaments and build class representatives */
    /* For n=8, we need 6880 classes. Random sampling should cover most. */

    /* For now, just compute H and deletion fingerprint for some samples */
    /* to test the speed */

    srand(42);
    int n_samples = 100;
    clock_t t0 = clock();

    for (int s = 0; s < n_samples; s++) {
        tournament_t T = ((uint32_t)rand() << 16) | (uint32_t)rand();
        T &= (TOTAL - 1);

        int scores[N];
        score_seq(T, scores);
        int c3 = count_c3(T);
        int64_t H = count_hp(T);

        /* Deletion fingerprint */
        int64_t del_h[N];
        for (int v = 0; v < N; v++)
            del_h[v] = hp_minus_v(T, v);

        /* Sort del_h */
        for (int i = 0; i < N-1; i++)
            for (int j = i+1; j < N; j++)
                if (del_h[i] > del_h[j]) {
                    int64_t tmp = del_h[i];
                    del_h[i] = del_h[j];
                    del_h[j] = tmp;
                }

        if (s < 5) {
            printf("  Sample %d: H=%lld, c3=%d, scores=(%d",
                   s, (long long)H, c3, scores[0]);
            for (int i = 1; i < N; i++) printf(",%d", scores[i]);
            printf("), del_H=(");
            for (int i = 0; i < N; i++) {
                if (i) printf(",");
                printf("%lld", (long long)del_h[i]);
            }
            printf(")\n");
        }
    }

    double elapsed = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("\n  %d samples in %.3fs (%.0f/s)\n", n_samples, elapsed, n_samples/elapsed);
    printf("  Estimated time for 268M tournaments: %.0f seconds\n",
           268435456.0 / (n_samples / elapsed));

    /* The key question: is computing (score, c3, deletion_H) for ALL 268M
       tournaments feasible? At ~10000/s (Python) or ~500000/s (C): */
    printf("  Estimated time at 500K/s: %.0f seconds = %.0f minutes\n",
           268435456.0 / 500000, 268435456.0 / 500000 / 60);

    return 0;
}
