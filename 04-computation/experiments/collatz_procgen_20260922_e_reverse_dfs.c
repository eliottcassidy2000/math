// Thread-deepening for Q2 (backward E, 3-adic): lift bad classes mod 3^(r+1) to mod 3^(r+2), test by DFS.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static const double L2=0.6931471805599453, L3=1.0986122886681098;
static uint64_t P3[45]; static long long nodes;
// x known mod 3^prec (prec>=1). returns 1 if a descending path exists.
static int dfs(uint64_t x, int prec, double cost, int depth){
  nodes++;
  if(depth>0 && cost < -1e-12) return 1;
  if(prec<=1) return 0;            // need x mod 9 to judge legality of any move
  if(cost - (prec-1)*L3 >= -1e-12) return 0;   // at most prec-1 more moves, each >= -log3
  uint64_t M=P3[prec], M1=P3[prec-1];
  int k0 = (x%3==1)?0:1;
  u128 pw = (k0==0)?1:2;
  for(int k=k0;;k+=2){
    double step=k*L2-L3;
    if(cost+step - (prec-2)*L3 >= -1e-12) break;  // even best future can't bring below 0
    uint64_t t=(uint64_t)((pw*(u128)x)%M);
    uint64_t y=((t+M-1)%M)/3; // mod 3^(prec-1)
    if(y%3!=0){ if(dfs(y%M1,prec-1,cost+step,depth+1)) return 1; }
    pw=(pw*4)%M;
  }
  return 0;
}
int main(int argc,char**argv){
  int r0=atoi(argv[1]), RMAX=atoi(argv[2]); const char*fn=argv[3];
  P3[0]=1; for(int i=1;i<45;i++) P3[i]=P3[i-1]*3;
  FILE*f=fopen(fn,"r"); size_t cap=10000000; uint64_t *cur=malloc(8*cap),*nxt=malloc(8*cap); size_t n=0; unsigned long long c; double w;
  while(fscanf(f,"%llu %lf",&c,&w)==2) cur[n++]=c; fclose(f);
  printf("r=%d bad=%zu\n",r0,n);
  for(int r=r0+1;r<=RMAX;r++){
    size_t n2=0; nodes=0; uint64_t Mold=P3[r]; // old modulus 3^(r0+1) -> classes mod 3^r
    for(size_t i=0;i<n;i++) for(int t=0;t<3;t++){
      uint64_t cc=cur[i]+t*Mold;  // lift to mod 3^(r+1)
      if(!dfs(cc,r+1,0.0,0)) nxt[n2++]=cc;
    }
    uint64_t*tmp=cur;cur=nxt;nxt=tmp;n=n2;
    printf("r=%d mod 3^%d: bad=%zu nodes=%lld\n",r,r+1,n,nodes); fflush(stdout);
    char fnm[64]; sprintf(fnm,"q2deep_r%d.txt",r); FILE*g=fopen(fnm,"w"); for(size_t i=0;i<n;i++) fprintf(g,"%llu\n",(unsigned long long)cur[i]); fclose(g);
  }
  return 0;
}
