/*
 * procgen_smallgraph_20260926_ham.c
 * Session collatz-procgen-20260922, lane "smallgraph" (2026-09-26).
 *
 * Hamiltonian path / cycle search on an explicit undirected graph read from stdin.
 * Used for the target-sum graphs G_S(n) (x ~ y iff x != y and x + y in S).
 * Every witness printed here is re-verified independently in Python by the runner;
 * only the EXHAUSTIVE modes ("NONE", counts marked COMPLETE) are used as proofs of
 * non-existence / exact counts, and they rely only on the sound pruning rules below.
 *
 * Input (stdin):  n m   then m lines "u v" (1-indexed, u != v, no duplicates).
 * Usage:
 *   ham exist_path  [-L nodelimit]              exhaustive existence of a Hamiltonian path
 *   ham count_paths [-L nodelimit]              exact number of Hamiltonian paths (up to reversal)
 *                                               + histogram of endpoint pairs
 *   ham exist_cycle [-L nodelimit]              exhaustive existence of a Hamiltonian cycle
 *   ham count_cycles [-L nodelimit]             exact number of Hamiltonian cycles (up to rotation, reversal)
 *   ham heur_path  [-S seed] [-R restarts] [-B budget] [-s start] [-f v1,v2,...]
 *                                               randomized Warnsdorff DFS; -f = vertices that must be INTERIOR
 *   ham heur_cycle [-S seed] [-R restarts] [-B budget]
 *
 * Sound pruning rules for a partial path ending at c with unvisited set U (|U| = r):
 *   du[u] = number of neighbours of u inside U (u in U).
 *   (P1) a vertex u in U with du[u] <= 1 is either the next vertex (adjacent to c) or the final endpoint;
 *        hence at most two such vertices, and if two, one of them is adjacent to c and is taken next.
 *   (P2) u in U with du[u] = 0 must be the next AND last vertex: allowed only if r = 1 and u ~ c.
 *   (P3) U must be connected (induced subgraph) and c must have a neighbour in U (r >= 1).
 *   (P4) a vertex flagged interior-only (or, in cycle mode, the closing condition) can never be the final
 *        endpoint: if du[u] <= 1 it must be the next vertex and needs du[u] = 1.
 *   (P5) must_end t: the final vertex must be t; any other u with du[u] <= 1 must be next.
 * All five rules only discard partial paths that cannot extend to a Hamiltonian path with the required
 * endpoint properties, so exhaustive modes are exact when they finish within the node limit.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n, m;
static int *deg0, *adjs, *adj;        /* CSR adjacency */
static char *isadj_c;                 /* scratch */
static char *vis;
static int *du;
static int *path;
static int plen;
static char *intonly;                 /* interior-only flags */
static int must_end = -1;
static long long nodes = 0, node_limit = 0;
static int aborted = 0;
static int mode_count = 0;            /* 1: count all, 0: stop at first */
static int found = 0;
static long long count = 0;
static long long count_cap = 0;         /* -K: stop once the raw count reaches the cap */
static int capped = 0;
static long long *ep;                 /* endpoint pair histogram (n+1)*(n+1) when n <= 400 */
static int cycle_mode = 0, cycle_v0 = -1;
static int *bestpath;
static int *stk, *stamp; static int stampv = 0;

static int adjacent(int a, int b) {
    for (int i = adjs[a]; i < adjs[a + 1]; i++) if (adj[i] == b) return 1;
    return 0;
}

static void visit(int w) { vis[w] = 1; for (int i = adjs[w]; i < adjs[w + 1]; i++) du[adj[i]]--; path[plen++] = w; }
static void unvisit(int w) { plen--; for (int i = adjs[w]; i < adjs[w + 1]; i++) du[adj[i]]++; vis[w] = 0; }

/* connectivity of U plus requirement that c has a neighbour in U */
static int connected_ok(int c, int r) {
    if (r == 0) return 1;
    int has = 0;
    for (int i = adjs[c]; i < adjs[c + 1]; i++) if (!vis[adj[i]]) { has = 1; break; }
    if (!has) return 0;
    int s = -1;
    for (int i = adjs[c]; i < adjs[c + 1]; i++) if (!vis[adj[i]]) { s = adj[i]; break; }
    stampv++;
    int top = 0, cnt = 0; stk[top++] = s; stamp[s] = stampv;
    while (top) {
        int u = stk[--top]; cnt++;
        for (int i = adjs[u]; i < adjs[u + 1]; i++) { int v = adj[i]; if (!vis[v] && stamp[v] != stampv) { stamp[v] = stampv; stk[top++] = v; } }
    }
    return cnt == r;
}

static void record_complete(void) {
    int a = path[0], b = path[n - 1];
    if (cycle_mode) {
        if (!adjacent(b, cycle_v0)) return;
        if (mode_count) { if (path[1] < b) { count++; } return; }
        found = 1; memcpy(bestpath, path, sizeof(int) * n); return;
    }
    if (must_end >= 0 && b != must_end) return;
    if (intonly[b]) return;
    if (mode_count) {
        count++;
        if (ep) { int x = a < b ? a : b, y = a < b ? b : a; ep[(long long)x * (n + 1) + y]++; }
        if (count_cap && count >= count_cap) { capped = 1; aborted = 1; }
        return;
    }
    found = 1; memcpy(bestpath, path, sizeof(int) * n);
}

/* candidate filter: returns number of candidates written to cand[], or -1 if dead */
static int candidates(int c, int *cand) {
    int r = n - plen;
    int low[3], nlow = 0;
    int forcednext = -1;
    for (int u = 1; u <= n; u++) {
        if (vis[u]) continue;
        int d = du[u];
        int finalok = 1;
        if (intonly[u]) finalok = 0;
        if (must_end >= 0 && u != must_end) finalok = 0;
        if (cycle_mode && !adjacent(u, cycle_v0)) finalok = 0;
        if (d == 0) {
            if (!(r == 1 && adjacent(u, c) && finalok)) return -1;
        }
        if (d <= 1) {
            if (!finalok) {
                /* must be next */
                if (d == 0 && r > 1) return -1;
                if (!adjacent(u, c)) return -1;
                if (forcednext >= 0 && forcednext != u) return -1;
                forcednext = u;
            }
            if (nlow < 3) low[nlow] = u;
            nlow++;
            if (nlow >= 3) return -1;
        }
    }
    if (cycle_mode && r >= 1) {
        int ok = 0;
        for (int i = adjs[cycle_v0]; i < adjs[cycle_v0 + 1]; i++) if (!vis[adj[i]]) { ok = 1; break; }
        if (!ok) return -1;
    }
    int nc = 0;
    if (forcednext >= 0) { cand[nc++] = forcednext; return nc; }
    if (nlow == 2) {
        /* one of them must be next; the other is the end */
        for (int k = 0; k < 2; k++) if (adjacent(low[k], c)) cand[nc++] = low[k];
        return nc;
    }
    for (int i = adjs[c]; i < adjs[c + 1]; i++) { int v = adj[i]; if (!vis[v]) cand[nc++] = v; }
    return nc;
}

static void dfs(int c) {
    if (aborted) return;
    nodes++;
    if (node_limit && nodes > node_limit) { aborted = 1; return; }
    if (plen == n) { record_complete(); return; }
    int r = n - plen;
    if (!connected_ok(c, r)) return;
    int *cand = (int *)malloc(sizeof(int) * (deg0[c] + 2));
    int nc = candidates(c, cand);
    for (int k = 0; k < nc; k++) {
        int w = cand[k];
        visit(w);
        dfs(w);
        unvisit(w);
        if (aborted || (found && !mode_count)) break;
    }
    free(cand);
}

/* ---------------- heuristic (randomized Warnsdorff with backtracking budget) ---------------- */
static uint64_t rng = 88172645463325252ULL;
static uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static long long budget, bnodes;

static int hdfs(int c) {
    if (plen == n) {
        int b = path[n - 1];
        if (cycle_mode) return adjacent(b, cycle_v0);
        if (intonly[b]) return 0;
        return 1;
    }
    bnodes++;
    if (bnodes > budget) return 0;
    int r = n - plen;
    if ((bnodes & 15) == 0 && !connected_ok(c, r)) return 0;
    int *cand = (int *)malloc(sizeof(int) * (deg0[c] + 2));
    int nc = candidates(c, cand);
    if (nc <= 0) { free(cand); return 0; }
    /* order by du ascending with random tie-break */
    int *key = (int *)malloc(sizeof(int) * nc);
    for (int k = 0; k < nc; k++) key[k] = du[cand[k]] * 1024 + (int)(xr() % 1024);
    for (int i = 1; i < nc; i++) { int kk = key[i], cc = cand[i], j = i - 1; while (j >= 0 && key[j] > kk) { key[j + 1] = key[j]; cand[j + 1] = cand[j]; j--; } key[j + 1] = kk; cand[j + 1] = cc; }
    int ok = 0;
    for (int k = 0; k < nc && !ok; k++) {
        int w = cand[k];
        visit(w);
        if (hdfs(w)) ok = 1; else unvisit(w);
        if (bnodes > budget) break;
    }
    free(cand); free(key);
    return ok;
}

static void usage(void) { fprintf(stderr, "usage: ham MODE [opts] < graph\n"); exit(2); }

int main(int argc, char **argv) {
    if (argc < 2) usage();
    const char *mode = argv[1];
    long long seed = 1, restarts = 100, bud = 200000; int start = -1; char *forbid = NULL;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-L") && i + 1 < argc) node_limit = atoll(argv[++i]);
        else if (!strcmp(argv[i], "-S") && i + 1 < argc) seed = atoll(argv[++i]);
        else if (!strcmp(argv[i], "-R") && i + 1 < argc) restarts = atoll(argv[++i]);
        else if (!strcmp(argv[i], "-B") && i + 1 < argc) bud = atoll(argv[++i]);
        else if (!strcmp(argv[i], "-s") && i + 1 < argc) start = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-f") && i + 1 < argc) forbid = argv[++i];
        else if (!strcmp(argv[i], "-e") && i + 1 < argc) must_end = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-K") && i + 1 < argc) count_cap = atoll(argv[++i]);
        else usage();
    }
    if (scanf("%d %d", &n, &m) != 2) { fprintf(stderr, "bad input\n"); return 2; }
    int *eu = (int *)malloc(sizeof(int) * (m + 1)), *ev = (int *)malloc(sizeof(int) * (m + 1));
    deg0 = (int *)calloc(n + 2, sizeof(int));
    for (int i = 0; i < m; i++) { if (scanf("%d %d", &eu[i], &ev[i]) != 2) { fprintf(stderr, "bad edge\n"); return 2; } deg0[eu[i]]++; deg0[ev[i]]++; }
    adjs = (int *)calloc(n + 2, sizeof(int));
    for (int u = 1; u <= n; u++) adjs[u + 1] = adjs[u] + deg0[u];
    adj = (int *)malloc(sizeof(int) * (2 * m + 1));
    int *fill = (int *)calloc(n + 2, sizeof(int));
    for (int i = 0; i < m; i++) { adj[adjs[eu[i]] + fill[eu[i]]++] = ev[i]; adj[adjs[ev[i]] + fill[ev[i]]++] = eu[i]; }
    vis = (char *)calloc(n + 2, 1); du = (int *)calloc(n + 2, sizeof(int)); path = (int *)calloc(n + 2, sizeof(int));
    bestpath = (int *)calloc(n + 2, sizeof(int)); intonly = (char *)calloc(n + 2, 1);
    stk = (int *)calloc(n + 2, sizeof(int)); stamp = (int *)calloc(n + 2, sizeof(int));
    isadj_c = NULL;
    for (int u = 1; u <= n; u++) du[u] = deg0[u];
    if (forbid) { char *p = forbid; while (*p) { int v = (int)strtol(p, &p, 10); if (v >= 1 && v <= n) intonly[v] = 1; if (*p == ',') p++; else if (*p) p++; } }

    if (n == 1) { printf("RESULT PATH 1 1\n"); return 0; }

    if (!strcmp(mode, "exist_path") || !strcmp(mode, "count_paths")) {
        mode_count = !strcmp(mode, "count_paths");
        if (mode_count && n <= 400) ep = (long long *)calloc((size_t)(n + 1) * (n + 1), sizeof(long long));
        int leaves[3], nl = 0, iso = 0;
        for (int u = 1; u <= n; u++) { if (deg0[u] == 0) iso = 1; if (deg0[u] == 1) { if (nl < 3) leaves[nl] = u; nl++; } }
        if (iso || nl >= 3) {
            if (mode_count) printf("COUNT 0 COMPLETE nodes=0 reason=%s\n", iso ? "isolated" : "3leaves");
            else printf("RESULT NONE nodes=0 reason=%s\n", iso ? "isolated" : "3leaves");
            return 0;
        }
        if (nl == 2) {
            if (!intonly[leaves[0]] && !intonly[leaves[1]] && (must_end < 0 || must_end == leaves[0] || must_end == leaves[1])) {
                int a = leaves[0], b = leaves[1];
                if (must_end == a) { a = leaves[1]; b = leaves[0]; }
                must_end = b;
                visit(a); dfs(a); unvisit(a);
            }
        } else if (nl == 1) {
            if (!intonly[leaves[0]]) {
                if (must_end == leaves[0]) { /* search from the other end: reverse roles */
                    for (int s = 1; s <= n && !aborted && !(found && !mode_count); s++) {
                        if (s == leaves[0] || intonly[s]) continue;
                        visit(s); dfs(s); unvisit(s);
                    }
                } else { visit(leaves[0]); dfs(leaves[0]); unvisit(leaves[0]); }
            }
        } else {
            /* all starts; every undirected path is met once from each end */
            int fixed_end = must_end;
            for (int s = 1; s <= n && !aborted && !(found && !mode_count); s++) {
                if (intonly[s]) continue;
                if (fixed_end >= 0 && s == fixed_end) continue;
                visit(s); dfs(s); unvisit(s);
            }
            if (fixed_end >= 0 && mode_count) count *= 2; /* compensate the halving below */
            if (fixed_end >= 0 && mode_count && ep) for (long long i = 0; i < (long long)(n + 1) * (n + 1); i++) ep[i] *= 2;
            if (mode_count) {
                /* every undirected path was found twice (once from each end) */
                count /= 2;
                if (ep) for (long long i = 0; i < (long long)(n + 1) * (n + 1); i++) ep[i] /= 2;
            }
        }
        if (mode_count) {
            printf("COUNT %lld %s nodes=%lld\n", count, capped ? "CAPPED" : (aborted ? "ABORTED" : "COMPLETE"), nodes);
            if (ep && !aborted) for (int x = 1; x <= n; x++) for (int y = x + 1; y <= n; y++) { long long c = ep[(long long)x * (n + 1) + y]; if (c) printf("EP %d %d %lld\n", x, y, c); }
        } else {
            if (found) { printf("RESULT PATH %d", n); for (int i = 0; i < n; i++) printf(" %d", bestpath[i]); printf("\n"); }
            else printf("RESULT %s nodes=%lld\n", aborted ? "UNKNOWN" : "NONE", nodes);
        }
        return 0;
    }
    if (!strcmp(mode, "exist_cycle") || !strcmp(mode, "count_cycles")) {
        mode_count = !strcmp(mode, "count_cycles");
        cycle_mode = 1;
        int v0 = 1;
        for (int u = 1; u <= n; u++) { if (deg0[u] < 2) { if (mode_count) printf("COUNT 0 COMPLETE nodes=0 reason=mindeg\n"); else printf("RESULT NONE nodes=0 reason=mindeg\n"); return 0; } if (deg0[u] < deg0[v0]) v0 = u; }
        if (n < 3) { printf(mode_count ? "COUNT 0 COMPLETE nodes=0 reason=small\n" : "RESULT NONE nodes=0 reason=small\n"); return 0; }
        cycle_v0 = v0;
        visit(v0); dfs(v0); unvisit(v0);
        if (mode_count) printf("COUNT %lld %s nodes=%lld\n", count, aborted ? "ABORTED" : "COMPLETE", nodes);
        else {
            if (found) { printf("RESULT CYCLE %d", n); for (int i = 0; i < n; i++) printf(" %d", bestpath[i]); printf("\n"); }
            else printf("RESULT %s nodes=%lld\n", aborted ? "UNKNOWN" : "NONE", nodes);
        }
        return 0;
    }
    if (!strcmp(mode, "heur_path") || !strcmp(mode, "heur_cycle")) {
        rng ^= (uint64_t)seed * 0x9E3779B97F4A7C15ULL; if (!rng) rng = 1;
        for (int k = 0; k < 20; k++) xr();
        budget = bud;
        if (!strcmp(mode, "heur_cycle")) {
            cycle_mode = 1; int v0 = 1;
            for (int u = 1; u <= n; u++) { if (deg0[u] < 2) { printf("RESULT FAIL reason=mindeg\n"); return 0; } if (deg0[u] < deg0[v0]) v0 = u; }
            cycle_v0 = v0;
        }
        for (long long t = 0; t < restarts; t++) {
            int s;
            if (cycle_mode) s = cycle_v0;
            else if (start >= 1) s = start;
            else {
                /* prefer a leaf; otherwise random non-interior-only vertex of small degree */
                int lf = -1; for (int u = 1; u <= n; u++) if (deg0[u] == 1 && !intonly[u]) { lf = u; break; }
                if (lf > 0 && (xr() % 4)) s = lf;
                else { do { s = 1 + (int)(xr() % n); } while (intonly[s]); }
            }
            if (intonly[s]) continue;
            bnodes = 0;
            visit(s);
            int ok = hdfs(s);
            if (ok) {
                printf("RESULT %s %d", cycle_mode ? "CYCLE" : "PATH", n); for (int i = 0; i < n; i++) printf(" %d", path[i]); printf(" restarts=%lld\n", t + 1);
                return 0;
            }
            while (plen > 0) unvisit(path[plen - 1]);
        }
        printf("RESULT FAIL restarts=%lld\n", restarts);
        return 0;
    }
    usage();
    return 0;
}
