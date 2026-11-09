/* lam.c - kind-pasteur-2026-11-08-S17.  Arc-strong connectivity lambda(T) = min over nonempty proper X of the
 * number of arcs from V\X into X (the fewest arcs whose reversal makes T non-strong), against H and first-creator
 * status.  Prints, for every lambda, the number of classes, of first creators, and min/max H overall and among
 * non-first-creators.  usage: lam <n> <txt|bin> <classfile> <S(n-1)> */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
#define MAXN 12
#define HCAP 100000
static unsigned char inP[HCAP];
static void decode(u64 c, int m, uint16_t *adj){ memset(adj,0,sizeof(uint16_t)*m); int b=0;
  for (int i=0;i<m;i++) for (int j=i+1;j<m;j++){ if (c>>b & 1) adj[i] |= 1u<<j; else adj[j] |= 1u<<i; b++; } }
static u64 hpaths(const uint16_t *adj, int m){
  static u64 F2[1<<MAXN][MAXN]; int F = 1<<m;
  for (int S=0;S<F;S++) for (int v=0;v<m;v++) F2[S][v]=0;
  for (int v=0;v<m;v++) F2[1<<v][v]=1;
  for (int S=1;S<F;S++) for (int v=0;v<m;v++) if (F2[S][v]){ u64 x=F2[S][v]; uint16_t o = adj[v] & ~S; while(o){ int u=__builtin_ctz(o); o&=o-1; F2[S|1<<u][u]+=x; } }
  u64 h=0; for (int v=0;v<m;v++) h+=F2[F-1][v]; return h; }
static int lambdaq(const uint16_t *adj, int m){ unsigned full=(1u<<m)-1; int best=1<<30;
  int mind = m; for (int v=0;v<m;v++){ int od=__builtin_popcount(adj[v]); if (od<mind) mind=od; if (m-1-od<mind) mind=m-1-od; }
  best = mind;   /* singletons */
  for (unsigned X=1; X<full && best>0; X++){ int k=__builtin_popcount(X); if (k<2 || k>m-2) continue; int b=0;
    for (int y=0;y<m && b<best;y++) if (!(X>>y&1)) b += __builtin_popcount(adj[y] & X);
    if (b < best) best=b; }
  return best; }
int main(int argc, char **argv){
  int n = atoi(argv[1]); int bin = !strcmp(argv[2],"bin"); FILE *fp = fopen(argv[3], bin?"rb":"r");
  { FILE *s=fopen(argv[4],"r"); long v; while (fscanf(s,"%ld",&v)==1) if (v>=0 && v<HCAP) inP[v]=1; fclose(s); }
  long long cnt[MAXN]={0}, fc[MAXN]={0}; long long hmin[MAXN], hmax[MAXN], nmin[MAXN], nmax[MAXN]; long long minfcH=1LL<<40, minfcHlam=-1;
  for (int l=0;l<MAXN;l++){ hmin[l]=nmin[l]=1LL<<40; hmax[l]=nmax[l]=-1; }
  char line[8192];
  for (;;){ u64 code; long long H; uint16_t adj[MAXN];
    if (bin){ unsigned char rec[34]; if (fread(rec,34,1,fp)!=1) break; memcpy(&code,rec,8); uint32_t h32; memcpy(&h32,rec+8,4); H=h32; decode(code,n,adj); }
    else { if (!fgets(line,sizeof line,fp)) break; unsigned long long c,h; if (sscanf(line,"%llu %llu",&c,&h)!=2) continue; code=c; decode(code,n,adj); H=(long long)hpaths(adj,n); }
    int l = lambdaq(adj,n); int isfc = !inP[H];
    cnt[l]++; if (isfc) fc[l]++;
    if (H<hmin[l]) hmin[l]=H; if (H>hmax[l]) hmax[l]=H;
    if (!isfc){ if (H<nmin[l]) nmin[l]=H; if (H>nmax[l]) nmax[l]=H; }
  }
  printf("order %d (lambda: classes, first creators, H range, H range of non-first-creators)\n", n);
  for (int l=0;l<MAXN;l++) if (cnt[l]) printf("  lambda=%d: %ld classes, %ld first creators (%.4f), H in [%ld,%ld], non-FC H in [%ld,%ld]\n",
     l, cnt[l], fc[l], (double)fc[l]/cnt[l], hmin[l], hmax[l], nmin[l]<(1L<<40)?nmin[l]:-1, nmax[l]);
  return 0; }
