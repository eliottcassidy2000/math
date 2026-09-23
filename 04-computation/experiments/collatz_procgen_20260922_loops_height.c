// collatz_procgen_20260922_loops_height.c
// Lexicographic DP (min K, then min height) for reverse E-loops through 1 with all values <= XMAX.
// Valid because every loop with K = K0(s) is K-minimal at every intermediate state (a smaller prefix K
// would give a loop with K < K0(s), impossible by the lower bound 2^(K+1) > 3^(s+1), s >= 2).
// Output per s: minK, K0(s), and H(s) = least possible max value on a loop with K = minK (exact if <= XMAX).
// Usage: loops_height SMAX XMAX     memory 12*(XMAX+1) bytes (16-bit K)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
int main(int argc,char**argv){
  int SMAX=atoi(argv[1]); uint64_t XMAX=strtoull(argv[2],0,10);
  uint16_t *cK=malloc(2*(XMAX+1)),*nK=malloc(2*(XMAX+1)); uint32_t *cH=malloc(4*(XMAX+1)),*nH=malloc(4*(XMAX+1));
  if(!cK||!nK||!cH||!nH){fprintf(stderr,"oom\n");return 1;}
  memset(cK,255,2*(XMAX+1)); cK[1]=0; cH[1]=1;
  for(int s=1;s<=SMAX;s++){
    memset(nK,255,2*(XMAX+1));
    for(uint64_t x=1;x<=XMAX;x++){
      if(cK[x]==0xFFFF) continue;
      int k0=(x%3==1)?0:1; unsigned __int128 p=((unsigned __int128)1)<<k0;
      for(int k=k0;k<=100;k+=2){
        unsigned __int128 t=p*x; if(t>(unsigned __int128)3*XMAX+1) break;
        uint64_t y=(uint64_t)((t-1)/3);
        if(y>=1 && y%3!=0){
          int K=cK[x]+k; uint32_t h=cH[x]; if(y>h) h=(uint32_t)y;
          if(K<0xFFFF && (K<nK[y] || (K==nK[y] && h<nH[y]))){ nK[y]=(uint16_t)K; nH[y]=h; }
        }
        p<<=2;
      }
    }
    int k0=(int)floorl((long double)(s+1)*log2l(3.0L));
    if(nK[1]==0xFFFF) printf("s=%4d K0=%4d none<=X\n",s,k0);
    else printf("s=%4d K0=%4d minK=%4d c=%.6f H=%u %s\n",s,k0,nK[1],pow(2.0,nK[1]-s*log2(3.0)),nH[1],nK[1]==k0?"=K0":"!=K0");
    fflush(stdout);
    uint16_t*t1=cK;cK=nK;nK=t1; uint32_t*t2=cH;cH=nH;nH=t2;
  }
  return 0;
}
