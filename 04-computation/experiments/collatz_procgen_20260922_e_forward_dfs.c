// Thread-deepening for Q1 bad classes via DFS with precision tracking.
// Input: bad classes at level m0 (from q1bad_dump.txt). For m=m0+1..MMAX: lift, test with DFS, keep bad.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
static const double L2=0.6931471805599453, L3=1.0986122886681098;
static long long nodes;
// returns 1 if a multiplicatively descending path exists
static int dfs(uint64_t xr, int prec, int a, int b){
  nodes++;
  double cost=a*L3-b*L2;
  if((a+b)>0 && cost< -1e-12) return 1;
  if(prec<=0) return 0;
  if(cost - prec*L2 >= -1e-12) return 0; // cannot get below 0 even halving every remaining bit
  uint64_t mask = (prec>=64)?~0ULL:((1ULL<<prec)-1);
  if(xr&1){ return dfs((3*xr+1)&mask,prec,a+1,b); }
  // even: halve first (cheap), then excursion
  if(dfs((xr>>1)&(mask>>1),prec-1,a,b+1)) return 1;
  return dfs((3*xr+1)&mask,prec,a+1,b);
}
int main(int argc,char**argv){
  int m0=atoi(argv[1]), MMAX=atoi(argv[2]); const char*fn=argv[3];
  FILE*f=fopen(fn,"r"); uint64_t *cur=malloc(sizeof(uint64_t)*5000000), *nxt=malloc(sizeof(uint64_t)*5000000); size_t n=0; unsigned long long c; double w;
  while(fscanf(f,"%llu %lf",&c,&w)==2) cur[n++]=c; fclose(f);
  printf("m=%d bad=%zu\n",m0,n);
  for(int m=m0+1;m<=MMAX;m++){
    size_t n2=0; nodes=0;
    for(size_t i=0;i<n;i++){
      for(int t=0;t<2;t++){
        uint64_t cc = cur[i] + ((uint64_t)t<<(m-1));
        if(!dfs(cc,m,0,0)) nxt[n2++]=cc;
      }
    }
    uint64_t *tmp=cur; cur=nxt; nxt=tmp; n=n2;
    printf("m=%d bad=%zu nodes=%lld\n",m,n,nodes); fflush(stdout);
    { char fnm[64]; sprintf(fnm,"q1deep_m%d.txt",m); FILE*g=fopen(fnm,"w"); for(size_t i=0;i<n;i++) fprintf(g,"%llu\n",(unsigned long long)cur[i]); fclose(g); }
  }
  return 0;
}
