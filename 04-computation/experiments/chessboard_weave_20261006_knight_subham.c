/* Count Hamiltonian cycles / paths of the knight graph induced on a subset of
 * the 8x8 board (chessboard-weave, 2026-10-06).
 *
 * usage: ./subham MASKHEX c     -> number of undirected Hamiltonian cycles
 *        ./subham MASKHEX p     -> number of directed Hamiltonian paths
 * Square (i,j) is bit 8*i+j of the mask.
 * Cycles: anchor v0 = a square of minimum degree; for each unordered pair {a<b}
 * of its neighbours count Hamiltonian a->b paths of G - v0 (exact DFS with
 * the same pruning as ..._knight_tours.c).  Each undirected cycle is counted once.
 * build: cc -O2 -o subham chessboard_weave_20261006_knight_subham.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static uint64_t nbr[64];
static int TGT;
static long long cnt;

static inline int pc(uint64_t x) { return __builtin_popcountll(x); }

static void dfs_cyc(int head, uint64_t U) {
    if (U == 0) { if (head == TGT) cnt++; return; }
    if (head == TGT) return;
    int forced = -1;
    uint64_t hb = 1ULL << head;
    for (uint64_t W = U; W; W &= W - 1) {
        int v = __builtin_ctzll(W);
        int eff = pc(nbr[v] & U) + ((nbr[v] & hb) ? 1 : 0);
        if (v == TGT) {
            if (eff < 1) return;
            if (eff == 1 && (nbr[v] & hb) && (U & ~(1ULL << TGT))) return;
        } else {
            if (eff < 2) return;
            if (eff == 2 && (nbr[v] & hb)) { if (forced >= 0 && forced != v) return; forced = v; }
        }
    }
    uint64_t cand = nbr[head] & U;
    if (forced >= 0) cand = 1ULL << forced;
    for (uint64_t W = cand; W; W &= W - 1) {
        int v = __builtin_ctzll(W);
        dfs_cyc(v, U & ~(1ULL << v));
    }
}

static void dfs_path(int head, uint64_t U) {
    if (U == 0) { cnt++; return; }
    uint64_t hb = 1ULL << head;
    int ones = 0, nU = pc(U);
    for (uint64_t W = U; W; W &= W - 1) {
        int v = __builtin_ctzll(W);
        uint64_t au = nbr[v] & U;
        int eff = pc(au) + ((nbr[v] & hb) ? 1 : 0);
        if (eff == 0) return;
        if (eff == 1) { if (++ones > 1) return; if (au == 0 && nU > 1) return; }
    }
    for (uint64_t W = nbr[head] & U; W; W &= W - 1) {
        int v = __builtin_ctzll(W);
        dfs_path(v, U & ~(1ULL << v));
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s MASKHEX c|p\n", argv[0]); return 1; }
    uint64_t M = strtoull(argv[1], NULL, 16);
    int d[8][2] = {{1,2},{2,1},{-1,2},{-2,1},{1,-2},{2,-1},{-1,-2},{-2,-1}};
    for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) {
        uint64_t m = 0;
        if (M >> (8 * i + j) & 1ULL)
            for (int k = 0; k < 8; k++) {
                int a = i + d[k][0], b = j + d[k][1];
                if (a >= 0 && a < 8 && b >= 0 && b < 8 && (M >> (8 * a + b) & 1ULL)) m |= 1ULL << (8 * a + b);
            }
        nbr[8 * i + j] = m;
    }
    cnt = 0;
    if (argv[2][0] == 'c') {
        int v0 = -1;
        for (int v = 0; v < 64; v++) if (M >> v & 1ULL) if (v0 < 0 || pc(nbr[v]) < pc(nbr[v0])) v0 = v;
        uint64_t rest = M & ~(1ULL << v0);
        for (uint64_t W = nbr[v0]; W; W &= W - 1) {
            int a = __builtin_ctzll(W);
            for (uint64_t W2 = nbr[v0] & ~((2ULL << a) - 1); W2; W2 &= W2 - 1) {
                int b = __builtin_ctzll(W2);
                TGT = b;
                dfs_cyc(a, rest & ~(1ULL << a));
            }
        }
        printf("mask %016llx: %d squares, anchor %d (deg %d): undirected Hamiltonian cycles = %lld\n",
               (unsigned long long)M, pc(M), v0, pc(nbr[v0]), cnt);
    } else {
        for (int x = 0; x < 64; x++) if (M >> x & 1ULL) dfs_path(x, M & ~(1ULL << x));
        printf("mask %016llx: %d squares: directed Hamiltonian paths = %lld\n", (unsigned long long)M, pc(M), cnt);
    }
    return 0;
}
