/* procgen_brackets_20260924_squaresum.c -- square-sum graphs and their target-set variants.
 * G_n: vertices 1..n, edge {a,b} (a != b) iff a+b lies in the target set:
 *   target 0: squares   1: odd squares (degenerate: a+b = 1 mod 8 splits the graph by residues mod 8)
 *   2: triangular numbers k(k+1)/2 (k >= 2)   3: squares and doubled squares   10+s: the squares k^2 with
 *   hash(k, s) even (a pseudo-random half of the squares, no congruence structure)
 * Modes:
 *   count LO HI TARGET LIMIT      number of undirected Hamiltonian paths and cycles for n in [LO,HI]
 *                                 (exhaustive DFS with dead-end pruning; stops counting paths at LIMIT)
 *   exist LO HI TARGET NODES      decide existence of a Hamiltonian path / cycle (exhaustive up to NODES DFS nodes,
 *                                 Warnsdorff order); prints YES / NO / UNKNOWN for path and cycle
 * n <= 127 (two 64-bit words per vertex set).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef struct { uint64_t w[2]; } Set;
static int n, target;
static Set adj[130];
static int deg0[130];
static inline int in(const Set *s, int v) { return (s->w[v >> 6] >> (v & 63)) & 1; }
static inline void add(Set *s, int v) { s->w[v >> 6] |= 1ULL << (v & 63); }
static inline void del(Set *s, int v) { s->w[v >> 6] &= ~(1ULL << (v & 63)); }
static inline int popc(Set a, Set b) { return __builtin_popcountll(a.w[0] & b.w[0]) + __builtin_popcountll(a.w[1] & b.w[1]); }

static uint64_t mix(uint64_t x) { x += 0x9E3779B97F4A7C15ULL; x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL; x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL; return x ^ (x >> 31); }
static int is_target(int s) {
    int r = 0; while ((r + 1) * (r + 1) <= s) r++;
    int sq = r * r == s;
    if (target == 0) return sq;
    if (target == 1) return sq && (r & 1);
    if (target == 2) { for (int k = 2; k * (k + 1) / 2 <= s; k++) if (k * (k + 1) / 2 == s) return 1; return 0; }
    if (target == 3) { if (sq) return 1; if (s % 2) return 0; int h = s / 2, t = 0; while ((t + 1) * (t + 1) <= h) t++; return t * t == h; }
    if (target >= 10) return sq && !(mix((uint64_t)r * 1000003ULL + (uint64_t)(target - 10)) & 1);
    return 0;
}
static void build(void) {
    for (int v = 0; v <= n + 1; v++) { adj[v].w[0] = adj[v].w[1] = 0; deg0[v] = 0; }
    for (int a = 1; a <= n; a++) for (int b = a + 1; b <= n; b++) if (is_target(a + b)) { add(&adj[a], b); add(&adj[b], a); deg0[a]++; deg0[b]++; }
}

static Set unv; static int path[130]; static int64_t npaths, ncycles, nodes, node_limit; static int want_cycle_only;
static int stop_first, found_path, found_cycle, aborted;
/* dead-end test: every unvisited vertex needs an available neighbour (unvisited or the current end);
   vertices with exactly one available neighbour must be the final vertex: at most one of them. */
static int dead(int cur, int remaining) {
    int ones = 0;
    Set av = unv; add(&av, cur);
    for (int wi = 0; wi < 2; wi++) {
        uint64_t x = unv.w[wi];
        while (x) {
            int v = (wi << 6) + __builtin_ctzll(x); x &= x - 1;
            int a = popc(adj[v], av);
            if (a == 0) return 1;
            if (a == 1) { if (in(&adj[v], cur) && remaining > 1) { /* v must come next and be last */ return 1; } if (++ones > 1) return 1; }
        }
    }
    return 0;
}
static void dfs(int cur, int depth) {
    if (aborted) return;
    if (++nodes > node_limit) { aborted = 1; return; }
    if (depth == n) {
        npaths++; found_path = 1;
        if (in(&adj[cur], path[0])) { ncycles++; found_cycle = 1; }
        return;
    }
    if (dead(cur, n - depth)) return;
    /* Warnsdorff order: fewest onward options first */
    int cand[130], sc[130], nc = 0;
    for (int wi = 0; wi < 2; wi++) { uint64_t x = adj[cur].w[wi] & unv.w[wi]; while (x) { int v = (wi << 6) + __builtin_ctzll(x); x &= x - 1; cand[nc] = v; sc[nc] = popc(adj[v], unv); nc++; } }
    for (int a = 1; a < nc; a++) { int v = cand[a], s = sc[a], b = a - 1; while (b >= 0 && sc[b] > s) { cand[b + 1] = cand[b]; sc[b + 1] = sc[b]; b--; } cand[b + 1] = v; sc[b + 1] = s; }
    for (int t = 0; t < nc; t++) {
        int v = cand[t]; del(&unv, v); path[depth] = v;
        dfs(v, depth + 1);
        add(&unv, v);
        if (stop_first && found_path && (!want_cycle_only || found_cycle)) return;
        if (aborted) return;
    }
}
static void reset_unv(void) { unv.w[0] = unv.w[1] = 0; for (int v = 1; v <= n; v++) add(&unv, v); }

int main(int argc, char **argv) {
    if (argc < 2) return 1;
    if (!strcmp(argv[1], "count")) {
        int lo = atoi(argv[2]), hi = atoi(argv[3]); target = atoi(argv[4]); int64_t limit = atoll(argv[5]);
        for (n = lo; n <= hi; n++) {
            build();
            npaths = ncycles = 0; nodes = 0; node_limit = limit; aborted = 0; stop_first = 0;
            for (int s = 1; s <= n && !aborted; s++) { reset_unv(); del(&unv, s); path[0] = s; dfs(s, 1); }
            /* ordered paths: each undirected path counted twice; each undirected cycle counted 2n times */
            printf("COUNT target %d n %d ham_paths %lld ham_cycles %lld %s nodes %lld\n", target, n, (long long)(npaths / 2), (long long)(ncycles / (2 * n)), aborted ? "ABORTED" : "exact", (long long)nodes);
            fflush(stdout);
        }
        return 0;
    }
    if (!strcmp(argv[1], "exist")) {
        int lo = atoi(argv[2]), hi = atoi(argv[3]); target = atoi(argv[4]); int64_t lim = atoll(argv[5]);
        for (n = lo; n <= hi; n++) {
            build();
            int mindeg = 1 << 30, leaves = 0, iso = 0; for (int v = 1; v <= n; v++) { if (deg0[v] < mindeg) mindeg = deg0[v]; leaves += deg0[v] == 1; iso += deg0[v] == 0; }
            const char *res[2];
            for (int mode = 0; mode < 2; mode++) {        /* 0 path, 1 cycle */
                nodes = 0; node_limit = lim; aborted = 0; stop_first = 1; want_cycle_only = mode; found_path = found_cycle = 0;
                npaths = ncycles = 0;
                if (mode == 1 && (mindeg < 2)) { res[1] = "NO(deg<2)"; continue; }
                if (mode == 0 && (iso > 0 && n > 1)) { res[0] = "NO(isolated)"; continue; }
                if (mode == 0 && leaves > 2) { res[0] = "NO(leaves>2)"; continue; }
                /* start vertices: for a path, a leaf if any, else every vertex; for a cycle, vertex 1 */
                int starts[130], ns = 0;
                if (mode == 1) starts[ns++] = 1;
                else { for (int v = 1; v <= n; v++) if (deg0[v] == 1) starts[ns++] = v; if (ns == 0) for (int v = 1; v <= n; v++) starts[ns++] = v; else ns = 1; }
                for (int t = 0; t < ns && !aborted; t++) {
                    reset_unv(); del(&unv, starts[t]); path[0] = starts[t]; dfs(starts[t], 1);
                    if (mode == 0 && found_path) break;
                    if (mode == 1 && found_cycle) break;
                }
                int ok = mode == 0 ? found_path : found_cycle;
                res[mode] = ok ? "YES" : (aborted ? "UNKNOWN" : "NO(search)");
            }
            printf("EXIST target %d n %d mindeg %d leaves %d isolated %d path %s cycle %s\n", target, n, mindeg, leaves, iso, res[0], res[1]);
            fflush(stdout);
        }
        return 0;
    }
    return 1;
}
