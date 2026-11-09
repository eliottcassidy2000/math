/* blue_lines.c - opus-2026-11-09-S87.
 * Every BLUE line {t, t + 1} (t in Fix(sigma)) at order n: exact H(t), H(t + 1), the packing numbers nu, and the
 * number of cyclic triangles; aggregated as counts of (H(t) mod 8, H(t+1) mod 8) and (nu(t), nu(t+1)).
 * Fix(sigma) has a basis of sigma-orbits of tiles (fixed tiles and 2-cycles), so blue tilings are enumerated directly.
 * mode 0: all blue lines; mode 1: SAMPLES random blue lines (for large n).
 * usage: blue_lines n mode [samples seed]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
#define MAXN 13
static int n, m, ti[80], tj[80], sg[80];
static u64 orb[80]; static int norb;
static u64 hpaths(const uint16_t *adj){ static u64 f[1<<MAXN][MAXN]; int F=1<<n;
  for (int S=0;S<F;S++) for (int v=0;v<n;v++) f[S][v]=0; for (int v=0;v<n;v++) f[1<<v][v]=1;
  for (int S=1;S<F;S++) for (int v=0;v<n;v++) if (f[S][v]){ u64 x=f[S][v]; uint16_t o=adj[v]&~S; while(o){ int u=__builtin_ctz(o); o&=o-1; f[S|1<<u][u]+=x; } }
  u64 h=0; for (int v=0;v<n;v++) h+=f[F-1][v]; return h; }
static int nu_of(const uint16_t *adj){   /* max number of vertex-disjoint cyclic triangles (<= 4 for n <= 13) */
  uint16_t tri[400]; int t=0;
  for (int a=0;a<n;a++) for (int b=a+1;b<n;b++) for (int c=b+1;c<n;c++){
    int ab=adj[a]>>b&1, bc=adj[b]>>c&1, ca=adj[c]>>a&1; if (ab==bc && bc==ca) tri[t++]=(1u<<a)|(1u<<b)|(1u<<c); }
  if (!t) return 0; int best=1;
  for (int i=0;i<t;i++) for (int j=i+1;j<t;j++) if (!(tri[i]&tri[j])){ if (best<2) best=2;
    for (int k=j+1;k<t;k++) if (!(tri[k]&(tri[i]|tri[j]))){ if (best<3) best=3;
      for (int l=k+1;l<t;l++) if (!(tri[l]&(tri[i]|tri[j]|tri[k]))){ return 4; } } }
  return best; }
static void build(u64 t, uint16_t *adj){ memset(adj,0,sizeof(uint16_t)*MAXN);
  for (int i=0;i+1<n;i++) adj[i] |= 1u<<(i+1);
  for (int b=0;b<m;b++){ int i=ti[b], j=tj[b]; if (t>>b & 1) adj[j] |= 1u<<i; else adj[i] |= 1u<<j; } }
int main(int argc, char **argv){
  n = atoi(argv[1]); int mode = atoi(argv[2]); long samples = argc>3 ? atol(argv[3]) : 0; unsigned seed = argc>4 ? atoi(argv[4]) : 1;
  int idx[MAXN][MAXN]; m=0; for (int i=0;i<n;i++) for (int j=i+2;j<n;j++){ ti[m]=i; tj[m]=j; idx[i][j]=m; m++; }
  for (int b=0;b<m;b++) sg[b] = idx[n-1-tj[b]][n-1-ti[b]];
  int done[80]={0}; norb=0; for (int b=0;b<m;b++) if (!done[b]){ orb[norb] = (1ULL<<b) | (1ULL<<sg[b]); done[b]=done[sg[b]]=1; norb++; }
  u64 all = (m==64)?~0ULL:((1ULL<<m)-1);
  static long hm8[8][8], nun[5][5], h4[4][4]; long lines=0;
  /* lines inside Fix(sigma): the all-ones vector is a sum of all orbits; pair t with t^all; enumerate t with the
     LAST orbit bit clear (exactly one of t, t^all has it) */
  long total = 1L << (norb-1);
  srand(seed);
  long iters = (mode==0) ? total : samples;
  for (long it=0; it<iters; it++){
    long code = (mode==0) ? it : (((long)rand()<<15) ^ rand()) % total;
    u64 t = 0; for (int k=0;k<norb-1;k++) if (code>>k & 1) t |= orb[k];
    uint16_t a1[MAXN], a2[MAXN]; build(t, a1); build(t^all, a2);
    u64 H1 = hpaths(a1), H2 = hpaths(a2); int v1 = nu_of(a1), v2 = nu_of(a2);
    int x = H1 & 7, y = H2 & 7; if (x > y){ int z=x; x=y; y=z; } hm8[x][y]++;
    int p = v1, q = v2; if (p > q){ int z=p; p=q; q=z; } nun[p][q]++;
    lines++;
  }
  printf("n=%d %s: blue lines %ld (of %ld)\n", n, mode==0?"all":"sampled", lines, total);
  printf("  (H mod 8 pairs): "); for (int x=0;x<8;x++) for (int y=x;y<8;y++) if (hm8[x][y]) printf("(%d,%d):%ld ", x, y, hm8[x][y]); printf("\n");
  long a=0,b=0,c=0; for (int x=0;x<8;x++) for (int y=x;y<8;y++){ int u=x&3, v=y&3; if (u>v){int z=u;u=v;v=z;} if (u==1&&v==1) a+=hm8[x][y]; else if (u==1&&v==3) b+=hm8[x][y]; else if (u==3&&v==3) c+=hm8[x][y]; }
  printf("  (H mod 4 pairs): (1,1):%ld (1,3):%ld (3,3):%ld\n", a, b, c);
  printf("  (nu pairs): "); for (int x=0;x<5;x++) for (int y=x;y<5;y++) if (nun[x][y]) printf("(%d,%d):%ld ", x, y, nun[x][y]); printf("\n");
  return 0;
}
