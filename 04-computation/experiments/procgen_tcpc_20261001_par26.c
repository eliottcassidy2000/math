/*
 * procgen_tcpc_20261001_par26.c  (TCPC lane, collatz-procgen-20260922, 2026-10-01)
 *
 * Memory-lean parity engine for the arc-HP counts of a digraph on n <= 27 vertices,
 * using only a FORWARD subset DP (one 2^n array of n-bit masks):
 *     H(D) mod 2   and   H(D - e) mod 2   for a list of arcs e,
 * so that c(e) = H(D) - H(D - e) (mod 2).  Independent of procgen_tcpc_20261001_hp.c
 * (which uses forward x backward convolution).
 *
 * stdin:  n, the n x n 0/1 matrix, then m, then m lines "u v" (arcs to test).
 * stdout: "H <H mod 2>" and one line "C u v <c(u->v) mod 2>" per tested arc.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static int n;
static int A[32][32];

static int hp_parity(uint32_t *F, int skip_u, int skip_v) {
    size_t N = (size_t)1 << n;
    uint32_t inA[32];
    for (int w = 0; w < n; w++) {
        inA[w] = 0;
        for (int u = 0; u < n; u++)
            if (u != w && (A[u][w] & 1) && !(u == skip_u && w == skip_v)) inA[w] |= 1u << u;
    }
    for (size_t S = 0; S < N; S++) F[S] = 0;
    for (int v = 0; v < n; v++) F[(size_t)1 << v] = 1u << v;
    for (size_t S = 1; S < N; S++) {
        uint32_t f = F[S];
        if (!f) continue;
        for (int w = 0; w < n; w++) {
            if (S >> w & 1) continue;
            if (__builtin_popcount(f & inA[w]) & 1) F[S | ((size_t)1 << w)] ^= 1u << w;
        }
    }
    return __builtin_popcount(F[N - 1]) & 1;
}

int main(void) {
    if (scanf("%d", &n) != 1 || n < 1 || n > 27) { fprintf(stderr, "bad n\n"); return 1; }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        if (scanf("%d", &A[i][j]) != 1) { fprintf(stderr, "bad matrix\n"); return 1; }
    int m;
    if (scanf("%d", &m) != 1) m = 0;
    uint32_t *F = malloc(((size_t)1 << n) * sizeof(uint32_t));
    if (!F) { fprintf(stderr, "alloc failed\n"); return 3; }
    int H = hp_parity(F, -1, -1);
    printf("H %d\n", H);
    fflush(stdout);
    for (int t = 0; t < m; t++) {
        int u, v;
        if (scanf("%d %d", &u, &v) != 2) break;
        int He = hp_parity(F, u, v);
        printf("C %d %d %d\n", u, v, (H ^ He) & 1);
        fflush(stdout);
    }
    free(F);
    return 0;
}
