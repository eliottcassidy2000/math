// collatz_procgen_20260922_mirror_chain.c
// Q1 mirror of half_chain.c.  For integers n = 7 mod 8 (the -1 thread at precision m = v_2(n+1) >= 3),
// find a shortest forward E-path (moves x -> 3x+1 always, x -> x/2 for even x) from n to a value < n,
// "shortest" = fewest halvings D (iterative deepening), values capped at CAP*n.
// Pruning (exact): a path with a multiplications can end below n after D halvings only if 3^a < 2^D, and
// from a state x with b halvings done, at least f_k(x mod 2^k) more multiplications occur before the next
// k halvings (k = min(D-b, P0)); f is the dimension lane's value function (oracle table, P0 bits).
// Recorded: D, the precision m, the multiplications A at the (m+1)-th halving (exit) versus A*(m),
// the largest excursion value/n, and "hostile re-landings": states after a halving b >= m+1 (after the
// exit) that are = 7 mod 8 (-1 thread, precision >= 3), resp. = 3 mod 4 (precision >= 2).
// Usage: mirror_chain N CAP DMAX [P0]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
static int P0; static uint8_t **T;
static inline int fP(int p, uint64_t y){ if(p==0) return 0; uint64_t m=(1ULL<<p)-1; y&=m;
  if(y&1) return 1+T[p][((3*y+1)&m)>>1]; return T[p][y>>1]; }
static uint64_t N0; static double CAP; static int D, AMAX, M; static uint64_t path[256]; static int pb[256], pa[256], plen;
static uint64_t best[256]; static int bestb[256], besta[256], bestlen;
static int THD[512];
static int dfs(uint64_t x,int b,int a,int len){
  path[len]=x; pb[len]=b; pa[len]=a;
  if(len>0 && pb[len]>pb[len-1] && x<N0){ bestlen=len; for(int i=0;i<=len;i++){best[i]=path[i];bestb[i]=pb[i];besta[i]=pa[i];} return 1; }
  if(b==D) return 0;
  int k=D-b; if(k>P0) k=P0;
  if(a+fP(k,x)>AMAX) return 0;
  if((double)x>CAP*(double)N0) return 0;
  if(x&1) return dfs(3*x+1,b,a+1,len+1);
  if(dfs(x>>1,b+1,a,len+1)) return 1;
  return dfs(3*x+1,b,a+1,len+1);
}
int main(int argc,char**argv){
  uint64_t N=strtoull(argv[1],0,10); CAP=atof(argv[2]); int DMAX=atoi(argv[3]); P0=argc>4?atoi(argv[4]):20;
  T=calloc(P0+1,sizeof(uint8_t*)); T[1]=malloc(1); T[1][0]=0;
  for(int p=2;p<=P0;p++){ uint64_t n=1ULL<<(p-1); T[p]=malloc(n); uint64_t mp=(1ULL<<p)-1;
    for(uint64_t i=0;i<n;i++){ uint64_t y=2*i; int bst=255; for(int j=0;2*j<bst;j++){ int v=fP(p-1,y>>1)+2*j; if(v<bst) bst=v; y=(9*y+4)&mp; } T[p][i]=bst; } }
  double TH=log(2.0)/log(3.0);
  for(int d=0;d<512;d++){ int a=(int)floor(d*TH); while(pow(3.0,a)>=pow(2.0,d)) a--; THD[d]=a; } // largest a with 3^a<2^d
  long cnt=0, fail=0, rev3=0, rev2=0, canon=0; long hist[256]={0}; long histm[64][64]={{0}};
  double maxexc=0; uint64_t maxexc_n=0; int maxd=0; uint64_t maxd_n=0; long revcount_hist[16]={0};
  for(uint64_t n=7;n<=N;n+=8){
    N0=n; M=__builtin_ctzll(n+1); int ok=0;
    for(D=M+1; D<=DMAX; D++){ AMAX=THD[D]; if(dfs(n,0,0,0)){ ok=1; break; } }
    cnt++;
    if(!ok){ fail++; continue; }
    hist[D]++; if(M<64 && D<64) histm[M][D]++;
    if(D>maxd){maxd=D;maxd_n=n;}
    double exc=0; int r3=0, r2=0, nre=0; int aexit=-1;
    for(int i=1;i<=bestlen;i++){ double r=(double)best[i]/(double)n; if(r>exc) exc=r;
      if(bestb[i]>bestb[i-1]){ // state right after a halving
        if(bestb[i]==M+1) aexit=besta[i];
        if(bestb[i]>=M+1 && i<bestlen){ if((best[i]&7)==7){r3=1;nre++;} if((best[i]&3)==3) r2=1; } } }
    if(exc>maxexc){maxexc=exc;maxexc_n=n;}
    rev3+=r3; rev2+=r2; if(nre>15) nre=15; revcount_hist[nre]++;
  }
  printf("n = 7 mod 8, n <= %llu: count=%ld, no descent within D<=%d halvings or cap %.0f: %ld\n",(unsigned long long)N,cnt,DMAX,CAP,fail);
  printf("max shortest D=%d (n=%llu); max excursion value/n on the shortest path found = %.3f (n=%llu)\n",maxd,(unsigned long long)maxd_n,maxexc,(unsigned long long)maxexc_n);
  printf("paths with a post-exit state = 7 mod 8 (hostile re-landing, precision >= 3): %ld (%.4f); = 3 mod 4 (precision >= 2): %ld (%.4f)\n",
         rev3,(double)rev3/(cnt-fail),rev2,(double)rev2/(cnt-fail));
  printf("number of post-exit precision>=3 re-landings per path:"); for(int i=0;i<16;i++) if(revcount_hist[i]) printf(" %d:%ld",i,revcount_hist[i]); printf("\n");
  printf("D histogram:"); for(int d=0;d<256;d++) if(hist[d]) printf(" %d:%ld",d,hist[d]); printf("\n");
  printf("D - (m+1) by precision m (rows m, entries excess:count):\n");
  for(int m=3;m<64;m++){ long s=0; for(int d=0;d<64;d++) s+=histm[m][d]; if(!s) continue; printf("  m=%2d (%8ld):",m,s);
    for(int d=0;d<64;d++) if(histm[m][d]) printf(" %d:%ld",d-(m+1),histm[m][d]); printf("\n"); }
  return 0;
}
