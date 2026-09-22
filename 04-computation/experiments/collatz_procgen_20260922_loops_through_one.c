// Reverse E-loops at 1: x0=1, x_i=(2^{k_i} x_{i-1}-1)/3 legal (x_i unit, positive integer), x_s=1.
// For each s find minimal K=sum k_i (values bounded by XMAX). Report ratio 2^K/3^s.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#define XMAX 20000000
int main(int argc,char**argv){
  int SMAX=atoi(argv[1]);
  // dist[x] = min K to reach x from 1 in exactly current s steps (layered), 255=inf
  uint8_t *cur=malloc(XMAX+1), *nxt=malloc(XMAX+1);
  memset(cur,255,XMAX+1); cur[1]=0;
  for(int s=1;s<=SMAX;s++){
    memset(nxt,255,XMAX+1);
    for(uint64_t x=1;x<=XMAX;x++){ if(cur[x]==255) continue; if(x%3==0) continue;
      int k0=(x%3==1)?0:1; uint64_t p=((uint64_t)1)<<k0;
      for(int k=k0;k<=60;k+=2){ unsigned __int128 t=(unsigned __int128)p*x; if(t>(unsigned __int128)3*XMAX+1) break;
        uint64_t y=(uint64_t)((t-1)/3); if(y>=1 && y%3!=0){ int K=cur[x]+k; if(K<255 && K<nxt[y]) nxt[y]=K; }
        p<<=2; }
    }
    uint8_t *tmp=cur; cur=nxt; nxt=tmp;
    if(cur[1]!=255){ double r=pow(2.0,cur[1])/pow(3.0,s); printf("s=%2d minK=%3d ratio 2^K/3^s=%.4f  %s\n",s,cur[1],r, r<3?"<3 OK":">=3"); }
    else printf("s=%2d no loop within XMAX\n",s);
    fflush(stdout);
  }
  return 0;
}
