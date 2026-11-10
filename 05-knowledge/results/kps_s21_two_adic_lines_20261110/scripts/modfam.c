/* modfam.c - two-module exponential families of Hamiltonian path counts.
 * For every tournament class K of order m (canonical upper-triangle code, first field of each input line) and every
 * ordered pair (x,y) of distinct vertices, substitute transitive modules TT_a at x and TT_b at y and compute
 *   h(a,b) = H(K[TT_a at x, TT_b at y]),  1 <= a,b <= AMAX.
 * Fit h = al + be*2^a + ga*2^b + de*2^(a+b) on a,b in {1,2} and test exactness on the rest of the grid.
 * usage: modfam m AMAX < cls_m.txt
 * output line: code x y al be ga de exact(0/1)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
typedef long long i64;
#define MAXV 16
static int adj[MAXV][MAXV];
static u64 dp[1<<MAXV][MAXV];           /* too big for MAXV=16 on stack; static: 65536*16*8 = 8MB */
static u64 hp(int n, int a[MAXV][MAXV]){
  int full = (1<<n) - 1; unsigned out[MAXV];
  for (int v=0; v<n; v++){ out[v]=0; for (int u=0; u<n; u++) if (a[v][u]) out[v] |= 1u<<u; }
  for (int m=1; m<=full; m++) for (int v=0; v<n; v++) dp[m][v]=0;
  for (int v=0; v<n; v++) dp[1<<v][v]=1;
  for (int m=1; m<=full; m++) for (int v=0; v<n; v++){ u64 c = dp[m][v]; if (!c) continue;
    unsigned nx = out[v] & ~m; while (nx){ int u = __builtin_ctz(nx); nx &= nx-1; dp[m|1u<<u][u] += c; } }
  u64 s=0; for (int v=0; v<n; v++) s += dp[full][v]; return s;
}
static int K[MAXV][MAXV];
static void decode(u64 code, int m){ int b=0; memset(K,0,sizeof K);
  for (int i=0;i<m;i++) for (int j=i+1;j<m;j++){ if (code>>b & 1) K[i][j]=1; else K[j][i]=1; b++; } }
/* build K with x -> TT_a, y -> TT_b (y<0: only x) */
static int build(int m, int x, int a, int y, int b){
  int orig[MAXV], idx[MAXV], n=0;   /* orig[v] = core vertex, idx = position inside its module */
  for (int v=0; v<m; v++){ int rep = (v==x)? a : (v==y)? b : 1; for (int i=0;i<rep;i++){ orig[n]=v; idx[n]=i; n++; } }
  for (int p=0;p<n;p++) for (int q=0;q<n;q++){ if (p==q){ adj[p][q]=0; continue; }
    if (orig[p]==orig[q]) adj[p][q] = idx[p] < idx[q]; else adj[p][q] = K[orig[p]][orig[q]]; }
  return n;
}
int main(int argc, char **argv){
  int m = atoi(argv[1]), AM = atoi(argv[2]); char line[4096];
  while (fgets(line, sizeof line, stdin)){
    u64 code = strtoull(line, NULL, 10); decode(code, m);
    for (int x=0;x<m;x++) for (int y=0;y<m;y++){ if (x==y) continue;
      i64 h[5][5];
      for (int a=1;a<=AM;a++) for (int b=1;b<=AM;b++){ int n = build(m,x,a,y,b); h[a][b] = (i64)hp(n, adj); }
      i64 D = h[2][2]-h[2][1]-h[1][2]+h[1][1];
      int exact = (D % 4 == 0);
      i64 de = D/4, be = (h[2][1]-h[1][1])/2 - 2*de, ga = (h[1][2]-h[1][1])/2 - 2*de, al = h[1][1]-2*be-2*ga-4*de;
      for (int a=1;a<=AM && exact;a++) for (int b=1;b<=AM;b++){
        i64 pred = al + be*(1LL<<a) + ga*(1LL<<b) + de*(1LL<<(a+b)); if (pred != h[a][b]){ exact=0; break; } }
      printf("%llu %d %d %lld %lld %lld %lld %d\n", code, x, y, al, be, ga, de, exact);
    }
  }
  return 0;
}
