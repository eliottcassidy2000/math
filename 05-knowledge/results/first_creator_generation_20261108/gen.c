/* gen.c - kind-pasteur-2026-11-08-S17.  Generation census: how every tournament class of order n arises from
 * first creators.  For each class T (order n) with H = H(T):
 *   FC   : H not in S(n-1)                          (a first creator: the minimal set)
 *   EXT  : some deletion T - v is a first creator of order n-1 (H(T-v) not in S(n-2)): T = (first creator) + v
 *   SUM  : T is not strong (an ordinal sum of smaller tournaments; H multiplies)
 *   SUB  : T has a nontrivial module (a substitution Q[T_1..T_k] of smaller tournaments)
 *   REST : none of the above (strong, prime, not a first creator, no first-creator deletion)
 * usage: gen <n> <txt|bin> <classfile> <S(n-2)> <S(n-1)> [maxprint] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
#define MAXN 12
#define HCAP 100000
static int n;
static unsigned char inPP[HCAP], inP[HCAP];
static void decode(u64 c, int m, uint16_t *adj){ memset(adj,0,sizeof(uint16_t)*m); int b=0;
  for (int i=0;i<m;i++) for (int j=i+1;j<m;j++){ if (c>>b & 1) adj[i] |= 1u<<j; else adj[j] |= 1u<<i; b++; } }
static u64 hall(const uint16_t *adj, int m, u64 *dh){
  static u64 F2[1<<MAXN][MAXN]; int F = 1<<m;
  for (int S=0;S<F;S++) for (int v=0;v<m;v++) F2[S][v]=0;
  for (int v=0;v<m;v++) F2[1<<v][v]=1;
  for (int S=1;S<F;S++) for (int v=0;v<m;v++) if (F2[S][v]){ u64 x=F2[S][v]; uint16_t o = adj[v] & ~S; while(o){ int u=__builtin_ctz(o); o&=o-1; F2[S|1<<u][u]+=x; } }
  u64 h=0; for (int v=0;v<m;v++) h+=F2[F-1][v];
  for (int x=0;x<m;x++){ int S=(F-1)&~(1<<x); u64 s=0; for(int v=0;v<m;v++) s+=F2[S][v]; dh[x]=s; }
  return h;
}
static int strongq(const uint16_t *adj, int m){ if(m<=1) return 1; uint16_t full=(1u<<m)-1, r=1, prev=0;
  while(r!=prev){ prev=r; for(int v=0;v<m;v++) if(r>>v&1) r|=adj[v]; } if(r!=full) return 0;
  uint16_t in[MAXN]={0}; for(int v=0;v<m;v++) for(int u=0;u<m;u++) if(adj[v]>>u&1) in[u]|=1u<<v;
  r=1; prev=0; while(r!=prev){ prev=r; for(int v=0;v<m;v++) if(r>>v&1) r|=in[v]; } return r==full; }
/* smallest module containing a and b; returns its mask */
static unsigned modclose(const uint16_t *adj, int m, int a, int b){
  unsigned M = (1u<<a)|(1u<<b); int changed = 1;
  while (changed){ changed = 0;
    for (int z=0; z<m; z++) if (!(M>>z&1)){ unsigned outz = adj[z] & M; if (outz && outz != M){ M |= 1u<<z; changed = 1; } }
  }
  return M;
}
static int primeq(const uint16_t *adj, int m){ unsigned full=(1u<<m)-1;
  for (int a=0;a<m;a++) for (int b=a+1;b<m;b++) if (modclose(adj,m,a,b) != full) return 0;
  return 1; }
/* arc-strong connectivity: min over nonempty proper X of the number of arcs from V\X into X; *bestX = a minimiser */
static int lambdaq(const uint16_t *adj, int m, unsigned *bestX){ unsigned full=(1u<<m)-1; int best=1<<30;
  for (unsigned X=1; X<full; X++){ int b=0; for (int y=0;y<m;y++) if (!(X>>y&1)) b += __builtin_popcount(adj[y] & X);
    if (b < best){ best=b; *bestX=X; } }
  return best; }
static u64 hsub(const uint16_t *adj, int m, unsigned X){ uint16_t a2[MAXN]; int idx[MAXN], k=0; for (int v=0;v<m;v++) if (X>>v&1) idx[k++]=v;
  for (int i=0;i<k;i++){ a2[i]=0; for (int j=0;j<k;j++) if (adj[idx[i]]>>idx[j]&1) a2[i]|=1u<<j; }
  u64 dh[MAXN]; return k ? hall(a2,k,dh) : 1; }
static unsigned char inS[13][HCAP];   /* S(k) for k <= n, loaded on demand */
static void loadset(const char *fn, unsigned char *a);
static int fcsub(const uint16_t *adj, int m, unsigned X){ int k=__builtin_popcount(X); u64 h=hsub(adj,m,X); if (k<=1) return 1; return !inS[k-1][h]; }
static void loadset(const char *fn, unsigned char *a){ FILE *fp=fopen(fn,"r"); long v; if(!fp){perror(fn); exit(1);} while (fscanf(fp,"%ld",&v)==1) if (v>=0 && v<HCAP) a[v]=1; fclose(fp); }
int main(int argc, char **argv){
  n = atoi(argv[1]); int bin = !strcmp(argv[2],"bin"); FILE *fp = fopen(argv[3], bin?"rb":"r");
  loadset(argv[4], inPP); loadset(argv[5], inP); int maxprint = argc>6 ? atoi(argv[6]) : 20;
  long tot=0, cFC=0, cEXT=0, cSUM=0, cSUB=0, cREST=0, nprime=0, nprimeFC=0, nstrong=0;
  long nonFC=0, nonFCstrong=0, nonFCstrongEXT=0, nonFCprime=0, nonFCprimeEXT=0;
  long strongNoFCdel=0, strongNoFCdelPrime=0;
  static long restH[HCAP]; char line[8192]; long lamHist[12]={0}, graftFF[4]={0}, restPend=0, restGraft=0, restOther=0;
  for (int k=1;k<=n-1 && k<=12;k++){ char fn[64]; snprintf(fn,sizeof fn,"S%d.txt",k); loadset(fn, inS[k]); }
  for (;;){
    u64 code; uint32_t H; u64 dh[MAXN]; int strong;
    uint16_t adj[MAXN];
    if (bin){ unsigned char rec[34]; if (fread(rec,34,1,fp)!=1) break; memcpy(&code,rec,8); memcpy(&H,rec+8,4);
      for (int v=0;v<n;v++){ uint16_t d; memcpy(&d, rec+12+2*v, 2); dh[v]=d; } strong = rec[32]; decode(code,n,adj); }
    else { if (!fgets(line,sizeof line,fp)) break; unsigned long long c,h; if (sscanf(line,"%llu %llu",&c,&h)!=2) continue; code=c;
      decode(code,n,adj); H=(uint32_t)hall(adj,n,dh); if (H!=h){ fprintf(stderr,"H mismatch\n"); return 2; } strong = strongq(adj,n); }
    tot++;
    int fc = !inP[H];
    int k = 0; for (int v=0;v<n;v++) if (!inPP[dh[v]]) k++;
    int pr = primeq(adj,n);
    if (strong) nstrong++;
    if (pr){ nprime++; if (fc) nprimeFC++; }
    if (fc) cFC++;
    else { nonFC++;
      if (strong){ nonFCstrong++; if (k) nonFCstrongEXT++; else { strongNoFCdel++; if (pr) strongNoFCdelPrime++; } }
      if (pr){ nonFCprime++; if (k) nonFCprimeEXT++; }
      if (k) cEXT++; else if (!strong) cSUM++; else if (!pr) cSUB++;
      else { cREST++; restH[H]++;
        unsigned X=0; int lam = lambdaq(adj,n,&X); if (lam < 12) lamHist[lam]++;
        { int pend=0; for (int v=0;v<n;v++){ int od=__builtin_popcount(adj[v]); if (od==1 || od==n-2) pend=1; }
          if (pend) restPend++; else if (lam==1) restGraft++; else { restOther++; if (restOther<=10) printf("RESTOTHER code=%llu H=%u lambda=%d\n", code, H, lam); } }
        if (lam == 1){ unsigned Y=((1u<<n)-1)&~X; int fa=fcsub(adj,n,X), fb=fcsub(adj,n,Y); graftFF[fa*2+fb]++;
          if (cREST <= maxprint) printf("REST code=%llu H=%u lambda=1 A(|%d|,H=%llu,%s) => B(|%d|,H=%llu,%s)\n", code, H, __builtin_popcount(X), hsub(adj,n,X), fa?"FC":"old", __builtin_popcount(Y), hsub(adj,n,Y), fb?"FC":"old"); }
        else if (cREST <= maxprint) printf("REST code=%llu H=%u lambda=%d\n", code, H, lam); }
    }
  }
  printf("order %d: classes %ld (strong %ld, prime %ld, prime first creators %ld)\n", n, tot, nstrong, nprime, nprimeFC);
  printf("  first creators %ld; non-first-creators %ld = EXT %ld + SUM %ld + SUB %ld + REST %ld\n", cFC, nonFC, cEXT, cSUM, cSUB, cREST);
  printf("  strong non-first-creators %ld (with a first-creator deletion %ld; without %ld, of which prime %ld)\n", nonFCstrong, nonFCstrongEXT, strongNoFCdel, strongNoFCdelPrime);
  printf("  prime non-first-creators %ld (with a first-creator deletion %ld)\n", nonFCprime, nonFCprimeEXT);
  printf("  REST values:"); for (int h=0;h<HCAP;h++) if (restH[h]) printf(" %d[%ld]", h, restH[h]); printf("\n");
  printf("  REST arc-strong connectivity histogram:"); for (int l=0;l<12;l++) if (lamHist[l]) printf(" lambda=%d:%ld", l, lamHist[l]); printf("\n");
  printf("  REST: pendant (a vertex of in- or out-degree 1) %ld, other single-back-arc graft %ld, neither %ld\n", restPend, restGraft, restOther);
  printf("  REST single-back-arc grafts A=>B: (old,old) %ld (old,FC) %ld (FC,old) %ld (FC,FC) %ld\n", graftFF[0], graftFF[1], graftFF[2], graftFF[3]);
  return 0;
}
