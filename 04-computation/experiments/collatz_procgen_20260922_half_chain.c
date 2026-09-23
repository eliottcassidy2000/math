// For integers m = 14 mod 27 (hostile 1/2-neighbourhood), find a shortest reverse-E path to a value < m
// (iterative deepening DFS, moves x -> (2^k x - 1)/3, k>=0, 2^k x=1 mod 3, result not = 0 mod 3, values <= CAP*m),
// and record: depth, max excursion ratio, and whether the path revisits a hostile class (1 or 14 mod 27, depth>=3).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static uint64_t M0; static double CAPR; static int best_depth; static uint64_t path[64]; static uint64_t bestpath[64]; static int found;
static int v3(uint64_t x){ int v=0; while(x%3==0 && x){x/=3;v++;} return v; }
static int dfs(uint64_t x,int d,int D){
  path[d]=x;
  if(d>0 && x<M0){ for(int i=0;i<=d;i++) bestpath[i]=path[i]; best_depth=d; return 1; }
  if(d==D) return 0;
  int k0=(x%3==1)?0:1; u128 t=(u128)x<<k0;
  for(int k=k0;k<=40;k+=2){ if(t>(u128)(CAPR*(double)M0)*3) break;
    uint64_t y=(uint64_t)((t-1)/3); if(y%3!=0 && y>=1){ if(dfs(y,d+1,D)) return 1; }
    t<<=2; }
  return 0;
}
int main(int argc,char**argv){
  uint64_t N=atoll(argv[1]); CAPR=atof(argv[2]); int DMAX=atoi(argv[3]);
  long cnt=0, revisit=0, fail=0; int maxd=0; uint64_t maxd_m=0; double maxexc=0; uint64_t maxexc_m=0;
  long hist[64]={0};
  for(uint64_t m=14; m<=N; m+=27){
    M0=m; found=0; int ok=0;
    for(int D=1; D<=DMAX; D++){ if(dfs(m,0,D)){ ok=1; break; } }
    cnt++;
    if(!ok){ fail++; continue; }
    hist[best_depth]++;
    if(best_depth>maxd){maxd=best_depth; maxd_m=m;}
    double exc=0; int rv=0;
    for(int i=1;i<best_depth;i++){ double r=(double)bestpath[i]/m; if(r>exc) exc=r; uint64_t z=bestpath[i]; if((z%27==14 && v3(2*z-1)>=3) || (z%27==1 && z!=1 && v3(z-1)>=3)) rv=1; }
    if(exc>maxexc){maxexc=exc; maxexc_m=m;}
    revisit+=rv;
  }
  printf("m=14 mod 27, m<=%llu: count=%ld fail(depth>%d or cap)=%ld; max shortest depth=%d (m=%llu); max excursion=%.3f (m=%llu); paths revisiting a hostile class: %ld\n",
    (unsigned long long)N,cnt,DMAX,fail,maxd,(unsigned long long)maxd_m,maxexc,(unsigned long long)maxexc_m,revisit);
  printf("depth histogram:"); for(int d=1;d<64;d++) if(hist[d]) printf(" %d:%ld",d,hist[d]); printf("\n");
  return 0;
}
