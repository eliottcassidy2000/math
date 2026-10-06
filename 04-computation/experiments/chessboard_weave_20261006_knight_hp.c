/* All directed Hamiltonian paths (open knight's tours) of an R x C board,
 * tabulated by (start x, second square y, end z)  (chessboard-weave, 2026-10-06).
 *
 * DFS from every start with exact pruning: an unvisited square needs >= 1
 * available neighbour (unvisited or head); at most one unvisited square may
 * have exactly 1 (it must be the end); a square whose only available
 * neighbour is the head must be the last square.
 *
 * Also tabulates, for every knight edge {u,v}, the number of directed Hamiltonian
 * paths using it (subtree leaf counts accumulated on DFS tree edges).
 *
 * usage: ./hp R C outfile [edgefile]  (outfile lines: "x y z count", squares i*C+j;
 *                                      edgefile lines: "u v count" for u<v, directed paths)
 * build: cc -O2 -o hp chessboard_weave_20261006_knight_hp.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int R, C, NV;
static uint64_t nbr[64];
static long long cnt[64][64];   /* for current start: [second][end] */
static int second;
static long long total = 0;
static long long ecnt[64][64];  /* directed-path usage of undirected edge (min,max) */

static inline int pc(uint64_t x) { return __builtin_popcountll(x); }

static long long dfs(int head, uint64_t U) {
    /* returns the number of completed paths below this node */
    if (U == 0) { cnt[second][head]++; total++; return 1; }
    uint64_t hb = 1ULL << head;
    int ones = 0;
    int nU = pc(U);
    for (uint64_t W = U; W; W &= W - 1) {
        int v = __builtin_ctzll(W);
        uint64_t au = nbr[v] & U;
        int eff = pc(au) + ((nbr[v] & hb) ? 1 : 0);
        if (eff == 0) return 0;
        if (eff == 1) {
            if (++ones > 1) return 0;
            if (au == 0 && nU > 1) return 0;   /* only the head: must be last now */
        }
    }
    long long sub = 0;
    for (uint64_t W = nbr[head] & U; W; W &= W - 1) {
        int v = __builtin_ctzll(W);
        long long s = dfs(v, U & ~(1ULL << v));
        if (s) { if (head < v) ecnt[head][v] += s; else ecnt[v][head] += s; }
        sub += s;
    }
    return sub;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s R C outfile\n", argv[0]); return 1; }
    R = atoi(argv[1]); C = atoi(argv[2]); NV = R * C;
    if (NV > 64) return 1;
    FILE *out = fopen(argv[3], "w");
    int d[8][2] = {{1,2},{2,1},{-1,2},{-2,1},{1,-2},{2,-1},{-1,-2},{-2,-1}};
    for (int i = 0; i < R; i++) for (int j = 0; j < C; j++) {
        uint64_t m = 0;
        for (int k = 0; k < 8; k++) {
            int a = i + d[k][0], b = j + d[k][1];
            if (a >= 0 && a < R && b >= 0 && b < C) m |= 1ULL << (a * C + b);
        }
        nbr[i * C + j] = m;
    }
    uint64_t all = (NV == 64) ? ~0ULL : ((1ULL << NV) - 1);
    for (int x = 0; x < NV; x++) {
        memset(cnt, 0, sizeof cnt);
        uint64_t U0 = all & ~(1ULL << x);
        for (uint64_t W = nbr[x]; W; W &= W - 1) {
            int y = __builtin_ctzll(W);
            second = y;
            long long s = dfs(y, U0 & ~(1ULL << y));
            if (s) { if (x < y) ecnt[x][y] += s; else ecnt[y][x] += s; }
        }
        for (int y = 0; y < NV; y++) for (int z = 0; z < NV; z++)
            if (cnt[y][z]) fprintf(out, "%d %d %d %lld\n", x, y, z, cnt[y][z]);
    }
    fclose(out);
    if (argc > 4) {
        FILE *eo = fopen(argv[4], "w");
        for (int u = 0; u < NV; u++) for (int v = u + 1; v < NV; v++)
            if (nbr[u] >> v & 1ULL) fprintf(eo, "%d %d %lld\n", u, v, ecnt[u][v]);
        fclose(eo);
    }
    printf("%dx%d directed Hamiltonian paths (open tours): %lld\n", R, C, total);
    return 0;
}
