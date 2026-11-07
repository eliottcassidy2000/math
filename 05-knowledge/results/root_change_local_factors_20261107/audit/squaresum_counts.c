/* squaresum_counts.c - auditor's check of the square-sum counts quoted in ATLAS.md: undirected Hamiltonian paths of
 * Q_n (vertices 1..n, edges x+y a square) for n = 15..32 and undirected Hamiltonian cycles for n = 30..37.
 * Plain DFS with a dead-end pruning (an unvisited vertex with no unvisited neighbour and not adjacent to the current end
 * kills the branch when more than one vertex is left). usage: squaresum_counts */
#include <stdio.h>
#include <string.h>
typedef unsigned long long ull;
static int n; static ull adj[64];
static ull paths, cycles; static int start;
static int issq(int x){ int r = 0; while (r*r < x) r++; return r*r == x; }
static int dead(ull vis, int cur, int left){
  if (left <= 1) return 0;
  ull un = ~vis & ((n == 64) ? ~0ULL : ((1ULL << n) - 1));
  for (ull m = un; m; m &= m - 1){ int v = __builtin_ctzll(m);
    if (!(adj[v] & un) && !(adj[cur] >> v & 1)) return 1; }
  return 0;
}
static void dfsp(int cur, ull vis, int left){
  if (left == 0){ paths++; return; }
  if (dead(vis, cur, left)) return;
  for (ull m = adj[cur] & ~vis; m; m &= m - 1){ int v = __builtin_ctzll(m); dfsp(v, vis | (1ULL << v), left - 1); }
}
static void dfsc(int cur, ull vis, int left){
  if (left == 0){ if (adj[cur] >> start & 1) cycles++; return; }
  for (ull m = adj[cur] & ~vis; m; m &= m - 1){ int v = __builtin_ctzll(m); dfsc(v, vis | (1ULL << v), left - 1); }
}
int main(void){
  for (n = 15; n <= 37; n++){
    memset(adj, 0, sizeof adj);
    for (int x = 1; x <= n; x++) for (int y = 1; y <= n; y++) if (x != y && issq(x + y)) adj[x-1] |= 1ULL << (y-1);
    long long pr = -1, cy = -1;
    if (n <= 32){ paths = 0; for (int s = 0; s < n; s++) dfsp(s, 1ULL << s, n - 1); pr = (long long)(paths / 2); }
    if (n >= 30){ cycles = 0; start = 0; dfsc(0, 1ULL, n - 1); cy = (long long)(cycles / 2); }
    printf("n=%d undirected Hamiltonian paths %lld cycles %lld\n", n, pr, cy); fflush(stdout);
  }
  return 0;
}
