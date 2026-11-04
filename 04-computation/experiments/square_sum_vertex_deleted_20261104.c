/* square_sum_vertex_deleted_20261104.c (kind-pasteur-2026-11-04-S10): square-sum graphs Q_n, Hamiltonian cycles
 * (capped count) and, for every v, whether Q_n - v is Hamiltonian / traceable; run structure along the unique cycle
 * of Q_32, Q_33. Usage: square_sum_vertex_deleted_20261104 NMIN NMAX [CAP] */
/* sq_hypo.c -- square-sum graphs Q_n on [1,n] (i ~ j iff i + j is a square).
 * For each n in [NMIN, NMAX]: count Hamiltonian cycles of Q_n (capped), and for every v test whether
 * Q_n - v has a Hamiltonian cycle ("v non-critical") and a Hamiltonian path; report min degree.
 * For uniquely Hamiltonian Q_n, test whether the non-critical vertices form one run of the unique cycle.
 * Usage: sq_hypo NMIN NMAX [CAP] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64;
static int N; static u64 adj[70];
static int issq(int x) { int r = 0; while ((r + 1) * (r + 1) <= x) r++; return r * r == x; }
static long cnt, cap; static int path[70], best[70];
/* count Hamiltonian cycles through vertex s0 (each undirected cycle counted twice) in induced graph on mask */
static u64 M;
static void dfs_cyc(int v, u64 vis, int len) {
    if (cnt >= cap) return;
    if (len == __builtin_popcountll(M)) { if (adj[v] >> path[0] & 1) { if (cnt == 0) memcpy(best, path, sizeof(int) * len); cnt++; } return; }
    /* prune: every unvisited vertex needs >= 2 neighbours among (unvisited + v + start) */
    u64 un = M & ~vis;
    u64 ends = (1ULL << v) | (1ULL << path[0]);
    for (u64 t = un; t; t &= t - 1) { int w = __builtin_ctzll(t); if (__builtin_popcountll(adj[w] & (un | ends)) < 2) return; }
    u64 nx = adj[v] & un;
    while (nx) { int w = __builtin_ctzll(nx); nx &= nx - 1; path[len] = w; dfs_cyc(w, vis | (1ULL << w), len + 1); if (cnt >= cap) return; }
}
static long hcycles(u64 mask, long c) {
    M = mask; cnt = 0; cap = c; if (__builtin_popcountll(mask) < 3) return 0;
    int s = __builtin_ctzll(mask), md = 99;  /* start at a min-degree vertex */
    for (u64 t = mask; t; t &= t - 1) { int w = __builtin_ctzll(t); int d = __builtin_popcountll(adj[w] & mask); if (d < md) { md = d; s = w; } }
    if (md < 2) return 0;
    path[0] = s; dfs_cyc(s, 1ULL << s, 1); return cnt;
}
static int found;
static void dfs_path(int v, u64 vis, int len) {
    if (found) return;
    if (len == __builtin_popcountll(M)) { found = 1; return; }
    u64 nx = adj[v] & M & ~vis;
    while (nx) { int w = __builtin_ctzll(nx); nx &= nx - 1; dfs_path(w, vis | (1ULL << w), len + 1); if (found) return; }
}
static int traceable(u64 mask) {
    M = mask; found = 0; int ones = 0;
    for (u64 t = mask; t; t &= t - 1) { int w = __builtin_ctzll(t); int d = __builtin_popcountll(adj[w] & mask); if (d == 0 && __builtin_popcountll(mask) > 1) return 0; if (d == 1) ones++; }
    if (ones > 2) return 0;
    for (u64 t = mask; t && !found; t &= t - 1) { int w = __builtin_ctzll(t);
        if (ones > 0 && __builtin_popcountll(adj[w] & mask) != 1) continue;
        dfs_path(w, 1ULL << w, 1); }
    return found;
}
int main(int argc, char **argv) {
    int lo = atoi(argv[1]), hi = atoi(argv[2]); long CAP = argc > 3 ? atol(argv[3]) : 1000000;
    for (int n = lo; n <= hi; n++) {
        N = n; memset(adj, 0, sizeof adj);
        for (int i = 1; i <= n; i++) for (int j = 1; j <= n; j++) if (i != j && issq(i + j)) adj[i - 1] |= 1ULL << (j - 1);
        u64 full = (n == 64) ? ~0ULL : ((1ULL << n) - 1);
        int mindeg = 99; for (int i = 0; i < n; i++) { int d = __builtin_popcountll(adj[i]); if (d < mindeg) mindeg = d; }
        long hc = hcycles(full, 2 * CAP) / 2;   /* undirected cycles */
        int cyc[70]; int haveCyc = hc > 0; if (haveCyc) memcpy(cyc, best, sizeof(int) * n);
        int tr = traceable(full);
        u64 NCc = 0, NCp = 0;
        for (int v = 0; v < n; v++) { u64 m = full & ~(1ULL << v); if (hcycles(m, 1) > 0) NCc |= 1ULL << v; if (traceable(m)) NCp |= 1ULL << v; }
        printf("n=%2d mindeg=%d traceable=%d Hcycles=%ld%s  #v:Q-v Hamiltonian=%d  #v:Q-v traceable=%d",
               n, mindeg, tr, hc, hc >= CAP ? "+" : "", __builtin_popcountll(NCc), __builtin_popcountll(NCp));
        if (hc == 1) {   /* run structure along the unique cycle */
            int runs = 0; for (int i = 0; i < n; i++) { int a = NCc >> cyc[i] & 1, b = NCc >> cyc[(i + n - 1) % n] & 1; if (a && !b) runs++; }
            printf("  [UH] runs of non-critical along cycle=%d; cycle:", runs);
            for (int i = 0; i < n; i++) printf(" %d%s", cyc[i] + 1, (NCc >> cyc[i] & 1) ? "*" : "");
        }
        if (!haveCyc && __builtin_popcountll(NCc) == n) printf("  HYPOHAMILTONIAN");
        if (!tr && __builtin_popcountll(NCp) == n) printf("  HYPOTRACEABLE");
        printf("\n"); fflush(stdout);
    }
    return 0;
}
