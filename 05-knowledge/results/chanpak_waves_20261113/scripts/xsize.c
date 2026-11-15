// distribution of |X| for unit out-sets (2W = H +- 1), all classes of order n
#include "armed_lib.h"
int main(int argc,char**argv){
  int n=atoi(argv[1]); CAPN=n+1; size_t NN=(size_t)1<<CAPN;
  he=malloc(sizeof(i64)*NN*CAPN); g=malloc(sizeof(i64)*NN*CAPN); Hc=malloc(sizeof(i64)*NN);
  char line[4096]; long long sz[2][MAXN]={{0}}; long long leafish[2]={0}, tot[2]={0};
  while(fgets(line,sizeof line,stdin)){
    if(line[0]=='>') continue;
    int out[MAXN]; parse(line,n,out); i64 Pi[MAXN*MAXN],H; build(n,out,Pi,&H);
    static int p[1<<16],m[1<<16]; int a,b; enum_units(n,Pi,H,p,&a,m,&b,1<<16);
    for(int s=0;s<2;s++){ int*L=s?m:p; int c=s?b:a; for(int i=0;i<c;i++){ int k=__builtin_popcount(L[i]); sz[s][k]++; tot[s]++;
        // leaf-like: X = {q} or V\{p}; also 'near-ordinal': X closed under... just record size
        if(k==1||k==n-1) leafish[s]++; } }
  }
  for(int s=0;s<2;s++){ printf("order %d %s units: total %lld, leaf-like %lld; |X| distribution:",n,s?"minus":"plus",tot[s],leafish[s]); for(int k=0;k<=n;k++) printf(" %lld",sz[s][k]); printf("\n"); }
  return 0;}
