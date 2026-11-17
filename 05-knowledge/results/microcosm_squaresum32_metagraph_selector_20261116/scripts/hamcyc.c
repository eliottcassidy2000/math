/* hamcyc.c (M-4): Hamiltonian paths (n = 15..34) and Hamiltonian cycles (n = 28..34) of the square-sum graph on 1..n.
   Paths: DFS from a degree-1 vertex if one exists (each path once), else from every start and halved.  Cycles: DFS from
   vertex 1 closing back to 1, halved.  Pruning: a vertex of degree 1 among the unvisited must be the last vertex (paths). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int n, adj[40][40], nb[40][40], deg[40], vis[40]; static long long cnt;
static void dfsP(int v, int depth) { if (depth == n) { cnt++; return; } for (int i = 0; i < deg[v]; i++) { int w = nb[v][i]; if (!vis[w]) { vis[w] = 1; dfsP(w, depth + 1); vis[w] = 0; } } }
static void dfsC(int v, int depth) { if (depth == n) { if (adj[v][1]) cnt++; return; } for (int i = 0; i < deg[v]; i++) { int w = nb[v][i]; if (!vis[w]) { vis[w] = 1; dfsC(w, depth + 1); vis[w] = 0; } } }
int main(void) {
    for (n = 15; n <= 34; n++) {
        memset(adj, 0, sizeof adj); memset(deg, 0, sizeof deg); memset(vis, 0, sizeof vis);
        for (int i = 1; i <= n; i++) for (int j = i + 1; j <= n; j++) { int s = i + j, r = 0; while (r * r < s) r++; if (r * r == s) { adj[i][j] = adj[j][i] = 1; nb[i][deg[i]++] = j; nb[j][deg[j]++] = i; } }
        int start = -1; for (int i = 1; i <= n; i++) if (deg[i] == 1) { start = i; break; }
        long long paths;
        if (start > 0) { cnt = 0; vis[start] = 1; dfsP(start, 1); vis[start] = 0; paths = cnt; }
        else { cnt = 0; for (int s = 1; s <= n; s++) { vis[s] = 1; dfsP(s, 1); vis[s] = 0; } paths = cnt / 2; }
        long long cycles = -1;
        if (n >= 28) { cnt = 0; vis[1] = 1; dfsC(1, 1); vis[1] = 0; cycles = cnt / 2; }
        printf("n = %2d: Hamiltonian paths %lld%s\n", n, paths, n >= 28 ? ({ static char b[64]; snprintf(b, 64, ", Hamiltonian cycles %lld", cycles); b; }) : "");
        fflush(stdout);
    }
    return 0;
}
