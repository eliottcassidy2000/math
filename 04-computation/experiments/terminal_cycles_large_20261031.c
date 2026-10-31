/* opus-2026-10-31-S69: terminal-form Hamiltonian cycles of the square-sum graph for large N (HYP-9438).

   Terminal form (session124): a Hamiltonian cycle of Q_N containing 8-1-3 (N even) or 15-1-3 (N odd), i.e. a
   Hamiltonian path of Q_N starting 1,3 and ending at 8 (N even) or 15 (N odd).
   Search: depth-first with Warnsdorff ordering (fewest onward options first), dead-end pruning, and a node budget.
   Output per N: FOUND (with the path verified independently), or UNDECIDED (budget exhausted). Never "NO".

   usage: terminal_cycles_large LO HI STEP BUDGET [END_OVERRIDE]
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXV 2048
static int N, deg[MAXV], *adj[MAXV], path[MAXV], vis[MAXV], avail[MAXV];
static long long nodes, budget;
static int target;

static int is_square(int x) { int r = 0; while ((long long)(r + 1) * (r + 1) <= x) r++; return r * r == x; }

static void build(int n) {
    N = n;
    for (int v = 1; v <= n; v++) {
        int c = 0, tmp[64];
        for (int k = 1; k * k <= 2 * n; k++) { int w = k * k - v; if (w >= 1 && w <= n && w != v) tmp[c++] = w; }
        free(adj[v]); adj[v] = malloc(sizeof(int) * (c ? c : 1)); memcpy(adj[v], tmp, sizeof(int) * c); deg[v] = c;
    }
}

static int dfs(int len) {
    if (++nodes > budget) return -1;
    int cur = path[len - 1];
    if (len == N) return cur == target;
    if (cur == target) return 0;
    /* dead-end pruning: unvisited vertices with 0 available neighbours, or more than one (besides target) with 1 */
    int dead = 0;
    for (int k = 0; k < deg[cur]; k++) { int w = adj[cur][k]; if (!vis[w]) avail[w]--; }
    for (int v = 1; v <= N; v++) if (!vis[v]) {
        int a = avail[v] + 0;   /* available = unvisited neighbours (cur now visited) ; cur counts as one entry option */
        int ent = a + 0;
        /* cur is the current end: a neighbour of cur that is unvisited can still be entered from cur */
        int adj_cur = 0;
        for (int k = 0; k < deg[v]; k++) if (adj[v][k] == cur) { adj_cur = 1; break; }
        ent = a + adj_cur;
        if (ent == 0) { dead = 2; break; }
        if (ent == 1 && v != target) dead++;
    }
    int res = 0;
    if (dead == 0) {   /* any non-target vertex with <= 1 connection, or a stranded target, is fatal */
        int cand[64], sc[64], c = 0;
        for (int k = 0; k < deg[cur]; k++) { int w = adj[cur][k]; if (!vis[w]) { cand[c] = w; sc[c] = avail[w]; c++; } }
        /* Warnsdorff: fewest onward options first; the target last unless it completes the path */
        for (int i = 1; i < c; i++) for (int j = i; j > 0 && (sc[j] < sc[j - 1]); j--) {
            int t = sc[j]; sc[j] = sc[j - 1]; sc[j - 1] = t; t = cand[j]; cand[j] = cand[j - 1]; cand[j - 1] = t; }
        for (int i = 0; i < c && res == 0; i++) {
            int w = cand[i];
            if (w == target && len + 1 != N) continue;
            vis[w] = 1; path[len] = w;
            res = dfs(len + 1);
            vis[w] = 0;
        }
    }
    for (int k = 0; k < deg[cur]; k++) { int w = adj[cur][k]; if (!vis[w]) avail[w]++; }
    return res;
}

static int verify(void) {
    int seen[MAXV]; memset(seen, 0, sizeof(seen));
    if (path[0] != 1 || path[1] != 3 || path[N - 1] != target) return 0;
    for (int i = 0; i < N; i++) { if (path[i] < 1 || path[i] > N || seen[path[i]]) return 0; seen[path[i]] = 1; }
    for (int i = 0; i + 1 < N; i++) if (!is_square(path[i] + path[i + 1])) return 0;
    return is_square(path[N - 1] + 1);   /* closing edge to 1 */
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s LO HI STEP BUDGET [END]\n", argv[0]); return 1; }
    int lo = atoi(argv[1]), hi = atoi(argv[2]), step = atoi(argv[3]); budget = atoll(argv[4]);
    int override = argc > 5 ? atoi(argv[5]) : 0;
    for (int n = lo; n <= hi; n += step) {
        build(n);
        target = override ? override : ((n % 2 == 0) ? 8 : 15);
        memset(vis, 0, sizeof(int) * (n + 2));
        for (int v = 1; v <= n; v++) avail[v] = deg[v];
        path[0] = 1; path[1] = 3; vis[1] = vis[3] = 1;
        for (int k = 0; k < deg[1]; k++) avail[adj[1][k]]--;
        nodes = 0;
        int r = dfs(2);
        if (r == 1) printf("N=%d target=%d FOUND verified=%s nodes=%lld\n", n, target, verify() ? "yes" : "NO", nodes);
        else if (r == -1) printf("N=%d target=%d UNDECIDED (budget %lld)\n", n, target, budget);
        else printf("N=%d target=%d EXHAUSTED (no path; complete search) nodes=%lld\n", n, target, nodes);
        fflush(stdout);
    }
    return 0;
}
