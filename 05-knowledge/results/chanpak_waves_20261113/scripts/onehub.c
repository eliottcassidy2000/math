// One-hub tournaments TT_k + v (v beats exactly the spine vertices in S). Count spectrum.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
int main(int argc,char**argv){
  int K=atoi(argv[1]); u64 LIM=strtoull(argv[2],0,10);
  unsigned char *cov=calloc(LIM/2+2,1);
  for(int k=1;k<=K;k++){
    for(u64 S=0;S<(1ULL<<k);S++){
      // bit i-1 of S <-> vertex i in S (v beats i)
      u64 H=0; int first=(S&1), last=!((S>>(k-1))&1);
      H=first+last;
      // transitions i notin S, i+1 in S
      for(int i=1;i<k;i++) if(!((S>>(i-1))&1) && ((S>>i)&1)) H++;
      // pairs b<a, b in S, a notin S: 2^(a-b-1)
      // accumulate: run over a, keep T = sum_{b in S, b<a} 2^(a-1-b)
      u64 Tacc=0; for(int a=1;a<=k;a++){ if(!((S>>(a-1))&1)) H+=Tacc; Tacc=Tacc*2+(((S>>(a-1))&1)?1:0); if(H>LIM) break; }
      if(H<=LIM && (H&1)) cov[H>>1]=1;
    }
  }
  int shown=0; u64 cnt=0; printf("one-hub (TT_k + v, k<=%d) misses among odd <= %llu:",K,LIM);
  for(u64 n=1;n<=LIM;n+=2) if(!cov[n>>1]){ cnt++; if(shown<60){printf(" %llu",n);shown++;} }
  printf("\n  total misses %llu of %llu\n",cnt,(LIM+1)/2);
  return 0;
}
