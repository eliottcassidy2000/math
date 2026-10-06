/* Enumerate ALL undirected closed knight's tours of an R x C board
 * (chessboard-weave, 2026-10-06).
 *
 * Corner a1 = (0,0) has degree 2 (neighbours (1,2),(2,1)), so every closed tour
 * contains (1,2)-(0,0)-(2,1).  Undirected closed tours <-> Hamiltonian paths
 * from (1,2) to (2,1) in G - (0,0) covering all other squares (one each).
 * DFS with exact pruning: an unvisited square v != target needs >= 2 available
 * neighbours (unvisited or the current head), the target needs >= 1; an
 * unvisited non-target neighbour of the head with exactly 2 available
 * neighbours forces the next move.
 *
 * usage: ./tours R C [outfile]   (writes one tour per line: square indices i*C+j,
 *        starting 0, (1,2), ..., (2,1)) ; prints the count.
 *        ./tours R C -e        prints the count, then "u v c(e)" for every knight
 *        edge u<v (c(e) = number of undirected closed tours through e).
 * build: cc -O2 -o tours chessboard_weave_20261006_knight_tours.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static int R, C, NV, A, B;
static uint64_t nbr[64];
static int path[64], plen;
static long long count = 0;
static FILE *out = NULL;
static int emode = 0;
static long long ec[64][64];

static inline int pc(uint64_t x) { return __builtin_popcountll(x); }

static void dfs(int head, uint64_t U) {
    /* U = unvisited squares (excluding head, which is visited) */
    if (U == 0) {
        /* head must be B (B is the last vertex) */
        if (head == B) {
            count++;
            if (emode) {
                int prev = 0;
                for (int k = 0; k < plen; k++) { int a = prev, b = path[k]; if (a > b) { int t = a; a = b; b = t; } ec[a][b]++; prev = path[k]; }
                { int a = 0, b = path[plen - 1]; if (a > b) { int t = a; a = b; b = t; } ec[a][b]++; }
            }
            if (out) {
                fprintf(out, "0");
                for (int k = 0; k < plen; k++) fprintf(out, " %d", path[k]);
                fprintf(out, "\n");
            }
        }
        return;
    }
    if (head == B) return; /* reached target too early */
    /* pruning scan */
    int forced = -1;
    uint64_t hb = 1ULL << head;
    for (uint64_t W = U; W; W &= W - 1) {
        int v = __builtin_ctzll(W);
        int eff = pc(nbr[v] & U) + ((nbr[v] & hb) ? 1 : 0);
        if (v == B) {
            if (eff < 1) return;
            if (eff == 1 && (nbr[v] & hb) && (U & ~(1ULL << B))) return;
        } else {
            if (eff < 2) return;
            if (eff == 2 && (nbr[v] & hb)) {
                if (forced >= 0 && forced != v) return;
                forced = v;
            }
        }
    }
    uint64_t cand = nbr[head] & U;
    if (forced >= 0) cand = 1ULL << forced;
    for (uint64_t W = cand; W; W &= W - 1) {
        int v = __builtin_ctzll(W);
        path[plen++] = v;
        dfs(v, U & ~(1ULL << v));
        plen--;
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s R C [outfile]\n", argv[0]); return 1; }
    R = atoi(argv[1]); C = atoi(argv[2]); NV = R * C;
    if (NV > 64 || R < 3 || C < 3) { fprintf(stderr, "bad size\n"); return 1; }
    if (argc > 3) { if (argv[3][0] == '-' && argv[3][1] == 'e') emode = 1; else out = fopen(argv[3], "w"); }
    int d[8][2] = {{1,2},{2,1},{-1,2},{-2,1},{1,-2},{2,-1},{-1,-2},{-2,-1}};
    for (int i = 0; i < R; i++) for (int j = 0; j < C; j++) {
        uint64_t m = 0;
        for (int k = 0; k < 8; k++) {
            int a = i + d[k][0], b = j + d[k][1];
            if (a >= 0 && a < R && b >= 0 && b < C) m |= 1ULL << (a * C + b);
        }
        nbr[i * C + j] = m;
    }
    A = 1 * C + 2; B = 2 * C + 1;
    uint64_t all = (NV == 64) ? ~0ULL : ((1ULL << NV) - 1);
    uint64_t U = all & ~1ULL & ~(1ULL << A);
    plen = 0; path[plen++] = A;
    dfs(A, U);
    printf("%dx%d undirected closed knight's tours: %lld\n", R, C, count);
    if (emode)
        for (int u = 0; u < NV; u++) for (int v = u + 1; v < NV; v++)
            if (nbr[u] >> v & 1ULL) printf("%d %d %lld\n", u, v, ec[u][v]);
    if (out) fclose(out);
    return 0;
}
