/* nu_lines.c - opus-2026-11-09-S87.
 * Every line {t, t + 1} of the tiling line metagraph at order n, classified by the packing numbers
 * (nu = max number of vertex-disjoint cyclic triangles = of disjoint cycles) of its two endpoints and by colour.
 * Also records which endpoint is strong (closed tilings contain the Hamiltonian cycle; THM-4729).
 * usage: nu_lines n SLICE NSLICES      (prints counts "nuA nuB blue count" with nuA <= nuB)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
#define MAXN 12
static int n, m, ti[66], tj[66], sg[66];
static int nu_of(const uint16_t *adj){
  uint16_t tri[300]; int t=0;
  for (int a=0;a<n;a++) for (int b=a+1;b<n;b++){ int ab = adj[a]>>b&1; for (int c=b+1;c<n;c++){
    int bc=adj[b]>>c&1, ca=adj[c]>>a&1; if (ab==bc && bc==ca) tri[t++]=(1u<<a)|(1u<<b)|(1u<<c); } }
  if (!t) return 0; int best=1;
  for (int i=0;i<t;i++) for (int j=i+1;j<t;j++) if (!(tri[i]&tri[j])){ if (best<2) best=2; if (n < 9) continue;
    for (int k=j+1;k<t;k++) if (!(tri[k]&(tri[i]|tri[j]))){ if (best<3) best=3; if (n < 12) continue;
      for (int l=k+1;l<t;l++) if (!(tri[l]&(tri[i]|tri[j]|tri[k]))) return 4; } }
  return best; }
int main(int argc, char **argv){
  n = atoi(argv[1]); int slice = atoi(argv[2]), ns = atoi(argv[3]);
  int idx[MAXN][MAXN]; m=0; for (int i=0;i<n;i++) for (int j=i+2;j<n;j++){ ti[m]=i; tj[m]=j; idx[i][j]=m; m++; }
  for (int b=0;b<m;b++) sg[b] = idx[n-1-tj[b]][n-1-ti[b]];
  u64 all = (1ULL<<m)-1, half = 1ULL<<(m-1);
  static long cnt[5][5][2];
  for (u64 t = slice; t < half; t += ns){
    u64 e[2] = { t, t ^ all }; int nu[2];
    for (int k=0;k<2;k++){ uint16_t adj[MAXN]; memset(adj,0,sizeof adj);
      for (int i=0;i+1<n;i++) adj[i] |= 1u<<(i+1);
      for (int b=0;b<m;b++){ if (e[k]>>b & 1) adj[tj[b]] |= 1u<<ti[b]; else adj[ti[b]] |= 1u<<tj[b]; }
      nu[k] = nu_of(adj); }
    u64 st=0; for (int b=0;b<m;b++) if (t>>b&1) st |= 1ULL<<sg[b];
    int blue = (st==t); int a=nu[0], c=nu[1]; if (a>c){int z=a;a=c;c=z;}
    cnt[a][c][blue]++;
  }
  for (int a=0;a<5;a++) for (int c=a;c<5;c++) for (int bl=0;bl<2;bl++) if (cnt[a][c][bl]) printf("%d %d %d %ld\n", a, c, bl, cnt[a][c][bl]);
  return 0;
}
