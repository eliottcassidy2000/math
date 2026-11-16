/* hampath_count.c: Hamiltonian paths of the square-sum graph on 1..n (i ~ j iff i + j is a perfect square), n = 15..25.
   Every such graph in this range has a vertex of degree 1, which must be an endpoint; a depth-first search from it counts
   each Hamiltonian path exactly once.  Prints the count and the multiset of far endpoints. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int n, adj[32][32], deg[32];
static long long cnt, endc[32];
static int vis[32];
static void dfs(int v, int depth) {
    if (depth == n) { cnt++; endc[v]++; return; }
    for (int w = 1; w <= n; w++) if (adj[v][w] && !vis[w]) { vis[w] = 1; dfs(w, depth + 1); vis[w] = 0; }
}
int main(void) {
    for (n = 15; n <= 25; n++) {
        memset(adj, 0, sizeof adj); memset(deg, 0, sizeof deg); cnt = 0; memset(endc, 0, sizeof endc); memset(vis, 0, sizeof vis);
        for (int i = 1; i <= n; i++) for (int j = i + 1; j <= n; j++) { int s = i + j, r = 0; while (r * r < s) r++; if (r * r == s) { adj[i][j] = adj[j][i] = 1; deg[i]++; deg[j]++; } }
        int start = -1;
        for (int i = 1; i <= n; i++) if (deg[i] == 1) { start = i; break; }
        if (start < 0) { printf("n=%d: no degree-1 vertex (not handled)\n", n); continue; }
        vis[start] = 1; dfs(start, 1);
        printf("n=%2d: start %2d (degree 1), Hamiltonian paths %lld, far endpoints:", n, start, cnt);
        for (int i = 1; i <= n; i++) if (endc[i]) printf(" %d x%lld", i, endc[i]);
        printf("\n");
    }
    return 0;
}
