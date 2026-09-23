// For each depth k (v_3(m-1)=k>=3), check exit A (loop len k-1, ratio<3) and exit B (loop len k-2, ratio<9/4),
// using minimal-K loops through 1 (reverse E-loops, values <= XMAX).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#define XMAX 20000000
int main(int argc,char**argv){
  int SMAX=atoi(argv[1]);
  uint8_t *cur=malloc(XMAX+1), *nxt=malloc(XMAX+1); memset(cur,255,XMAX+1); cur[1]=0;
  int minK[200]; for(int i=0;i<200;i++) minK[i]=-1;
  for(int s=1;s<=SMAX;s++){
    memset(nxt,255,XMAX+1);
    for(uint64_t x=1;x<=XMAX;x++){ if(cur[x]==255||x%3==0) continue; int k0=(x%3==1)?0:1; uint64_t p=((uint64_t)1)<<k0;
      for(int k=k0;k<=60;k+=2){ unsigned __int128 t=(unsigned __int128)p*x; if(t>(unsigned __int128)3*XMAX+1) break; uint64_t y=(uint64_t)((t-1)/3);
        if(y>=1 && y%3!=0){ int K=cur[x]+k; if(K<255 && K<nxt[y]) nxt[y]=K; } p<<=2; } }
    uint8_t*tmp=cur;cur=nxt;nxt=tmp; minK[s]=(cur[1]==255)?-1:cur[1];
  }
  int allok=1;
  for(int k=3;k<=SMAX+1;k++){
    int sA=k-1, sB=k-2; double rA= (sA>=1&&minK[sA]>=0)? pow(2,minK[sA])/pow(3,sA):1e9; double rB=(sB>=1&&minK[sB]>=0)? pow(2,minK[sB])/pow(3,sB):1e9;
    int A=rA<3, B=rB<2.25; if(!(A||B)) allok=0;
    printf("k=%2d  A: loop len %2d ratio %.4f %s   B: loop len %2d ratio %.4f %s  => %s\n",k,sA,rA,A?"ok":"--",sB,rB,B?"ok":"--",(A||B)?"ESCAPES":"FAIL");
  }
  printf("all depths 3..%d escape: %s\n",SMAX+1, allok?"YES":"NO");
  return 0;
}
