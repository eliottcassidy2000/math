/* lam11.c - kind-pasteur-2026-11-08-S17.  Exhaustive search of the 11-tournaments with arc-strong connectivity
 * lambda >= 3 and H <= HMAXQ (= max H(10) = 15745).
 * Every such T' has minimum in- and out-degree >= 3, so T = T' - v has minimum in- and out-degree >= 2 for every v;
 * extending every order-10 class with that property by every pattern therefore reaches every such T'.
 * H(T + v_x) comes from the quadratic attachment form (as in ext.c); lambda(T') >= 3 is then tested directly.
 * usage: lam11 <n> <txt|bin> <classfile> <lo> <hi> <S(n) values> <hmax> <outprefix>  -> <outprefix>.cnt (u64 per H) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
#define MAXN 12
#define HCAPQ 100000
static int n = 10; static int HMAXQ = 15745;
static uint32_t f[1<<11][MAXN], g[1<<11][MAXN];
static u64 cnt[HCAPQ];
static unsigned char inS10[100000];
static void decode(u64 c, int m, uint16_t *adj){ memset(adj,0,sizeof(uint16_t)*m); int b=0;
  for (int i=0;i<m;i++) for (int j=i+1;j<m;j++){ if (c>>b & 1) adj[i] |= 1u<<j; else adj[j] |= 1u<<i; b++; } }
static void dps(const uint16_t *adj){
  int F = 1<<n;
  for (int S=0;S<F;S++) for (int v=0;v<n;v++){ f[S][v]=0; g[S][v]=0; }
  for (int v=0;v<n;v++){ f[1<<v][v]=1; g[1<<v][v]=1; }
  uint16_t in[MAXN]; for (int v=0;v<n;v++){ in[v]=0; for (int u=0;u<n;u++) if (adj[u]>>v&1) in[v]|=1u<<u; }
  for (int S=1;S<F;S++) for (int v=0;v<n;v++){
    if (f[S][v]){ uint32_t x=f[S][v]; uint16_t o = adj[v] & ~S; while(o){ int u=__builtin_ctz(o); o&=o-1; f[S|1<<u][u]+=x; } }
    if (g[S][v]){ uint32_t x=g[S][v]; uint16_t o = in[v] & ~S; while(o){ int u=__builtin_ctz(o); o&=o-1; g[S|1<<u][u]+=x; } }
  }
}
static int32_t form(int U, int32_t *L, int32_t M[MAXN][MAXN]){
  int32_t W[MAXN][MAXN]; memset(W,0,sizeof W);
  for (int S=(U-1)&U; S>0; S=(S-1)&U){ int R = U & ~S;
    for (int a=0;a<n;a++) if ((S>>a&1) && f[S][a]){ int32_t fa=f[S][a];
      for (int b=0;b<n;b++) if ((R>>b&1) && g[R][b]) W[a][b] += fa * (int32_t)g[R][b]; } }
  int32_t h0=0; for (int b=0;b<n;b++) if (U>>b&1) h0 += g[U][b];
  for (int a=0;a<n;a++) if (U>>a&1){ int32_t s=(int32_t)f[U][a]-(int32_t)g[U][a];
    for (int b=0;b<n;b++) if ((U>>b&1) && b!=a){ s += W[a][b]; M[a][b] = W[a][b]+W[b][a]; } L[a]=s; }
  return h0;
}
/* lambda(T') >= 3 ?  (m = 11; early exit on a cut with at most 2 back arcs) */
static int lam3(const uint16_t *adj, int m){ unsigned full=(1u<<m)-1;
  for (int v=0;v<m;v++){ int od=__builtin_popcount(adj[v]); if (od<3 || m-1-od<3) return 0; }
  for (unsigned X=1; X<full; X++){ int k=__builtin_popcount(X); if (k<2 || k>m-2) continue; int b=0;
    for (int y=0;y<m;y++) if (!(X>>y&1)){ b += __builtin_popcount(adj[y] & X); if (b>2) break; }
    if (b<=2) return 0; }
  return 1; }
int main(int argc, char **argv){
  n = atoi(argv[1]); int bin = !strcmp(argv[2],"bin"); FILE *fp = fopen(argv[3], bin?"rb":"r"); long lo=atol(argv[4]), hi=atol(argv[5]);
  { FILE *s=fopen(argv[6],"r"); long v; while (fscanf(s,"%ld",&v)==1) if (v>=0 && v<100000) inS10[v]=1; fclose(s); }
  HMAXQ = atoi(argv[7]); const char *outp = argv[8]; char line[8192];
  long idx=-1, used=0, tested=0, found=0, nwit=0; long long minH = 1LL<<40; u64 minCode=0; unsigned minX=0;
  int full=(1<<n)-1;
  if (bin && lo>0) fseek(fp, (long)lo*34, SEEK_SET), idx = lo-1;
  for (;;){ u64 code;
    if (bin){ unsigned char rec[34]; if (fread(rec,34,1,fp)!=1) break; memcpy(&code,rec,8); }
    else { if (!fgets(line,sizeof line,fp)) break; unsigned long long c,hh; if (sscanf(line,"%llu %llu",&c,&hh)!=2) continue; code=c; }
    idx++; if (idx<lo) continue; if (idx>=hi) break;
    uint16_t adj[MAXN]; decode(code,n,adj);
    int od[MAXN], ok=1; for (int v=0;v<n;v++){ od[v]=__builtin_popcount(adj[v]); if (od[v]<2 || od[v]>n-3) ok=0; }
    if (!ok) continue;
    used++; dps(adj); int32_t L[MAXN], M[MAXN][MAXN]; int32_t h = form(full, L, M);
    int32_t D[MAXN]; for (int a=0;a<n;a++) D[a]=L[a]; unsigned x=0;
    for (unsigned i=0;; i++){
      int px = __builtin_popcount(x);
      if (px>=3 && px<=n-3 && h<=HMAXQ){
        int dok=1; for (int a=0;a<n;a++){ int o = od[a] + (x>>a&1); if (o<3 || o>n-3){ dok=0; break; } }
        if (dok){ uint16_t a2[MAXN]; for (int v=0;v<n;v++) a2[v]=adj[v]; a2[n]=0;
          for (int v=0;v<n;v++){ if (x>>v&1) a2[v] |= 1u<<n; else a2[n] |= 1u<<v; }
          tested++;
          if (lam3(a2, n+1)){ found++; cnt[h]++; if (h<minH){ minH=h; minCode=code; minX=x; }
            if (nwit < 1000){ nwit++; printf("W %llu %u %d %s\n", code, x, h, inS10[h] ? "inS10" : "new"); } } } }
      if (i==(unsigned)full) break;
      int c=__builtin_ctz(i+1);
      if (!(x>>c&1)){ h += D[c]; x |= 1u<<c; for (int b=0;b<n;b++) if (b!=c) D[b] -= M[b][c]; }
      else { x &= ~(1u<<c); for (int b=0;b<n;b++) if (b!=c) D[b] += M[b][c]; h -= D[c]; }
    }
  }
  char fn[1024]; snprintf(fn,sizeof fn,"%s.cnt",outp); FILE *o=fopen(fn,"wb"); fwrite(cnt,sizeof(u64),HCAPQ,o); fclose(o);
  printf("slice [%ld,%ld): parents used %ld, extensions tested %ld, lambda>=3 with H<=%d: %ld; min H %lld (parent %llu, pattern %u)\n",
         lo, hi, used, tested, HMAXQ, found, found?minH:-1, minCode, minX);
  return 0; }
