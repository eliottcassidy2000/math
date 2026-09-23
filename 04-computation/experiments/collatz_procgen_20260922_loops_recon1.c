// collatz_procgen_20260922_loops_recon1.c
// Memory-light reconstruction of ONE reverse E-loop through 1 of length S with minimal K, values <= XMAX,
// optionally through a hub (start/end point h instead of 1).  Checkpoints every B layers; when backtracking,
// the layers of each block are recomputed from its checkpoint.  Memory ~ (S/B + B) * 2 * (XMAX+1) bytes.
// Usage: loops_recon1 S XMAX [B] [h]     prints "s=S K0=.. minK=.. ratio=.. max=.. k:k_1,...,k_S"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
static uint64_t X;
static void step(const uint16_t *cur, uint16_t *nxt){
  memset(nxt,255,2*(X+1));
  for(uint64_t x=1;x<=X;x++){
    if(cur[x]==0xFFFF) continue;
    int k0=(x%3==1)?0:1; unsigned __int128 p=((unsigned __int128)1)<<k0;
    for(int k=k0;k<=100;k+=2){
      unsigned __int128 t=p*x; if(t>(unsigned __int128)3*X+1) break;
      uint64_t y=(uint64_t)((t-1)/3);
      if(y>=1 && y%3!=0){ int K=cur[x]+k; if(K<0xFFFF && K<nxt[y]) nxt[y]=(uint16_t)K; }
      p<<=2;
    }
  }
}
int main(int argc,char**argv){
  int S=atoi(argv[1]); X=strtoull(argv[2],0,10); int B=(argc>3)?atoi(argv[3]):32; uint64_t h=(argc>4)?strtoull(argv[4],0,10):1;
  int nck=S/B+1;
  uint16_t **ck=malloc(sizeof(uint16_t*)*(nck+1));
  uint16_t *a=malloc(2*(X+1)),*b=malloc(2*(X+1));
  memset(a,255,2*(X+1)); a[h]=0;
  ck[0]=malloc(2*(X+1)); memcpy(ck[0],a,2*(X+1));
  for(int i=1;i<=S;i++){ step(a,b); uint16_t*t=a;a=b;b=t; if(i%B==0 && i/B<=nck){ ck[i/B]=malloc(2*(X+1)); memcpy(ck[i/B],a,2*(X+1)); } }
  int mk=a[h]; int k0=(int)floorl((long double)(S+1)*log2l(3.0L));
  if(mk==0xFFFF){ printf("s=%d no loop with values <= %llu\n",S,(unsigned long long)X); return 0; }
  // backtrack block by block
  int *ks=malloc(sizeof(int)*(S+1)); uint64_t x=h; int K=mk; uint64_t vmax=h;
  uint16_t **blk=malloc(sizeof(uint16_t*)*(B+1)); for(int j=0;j<=B;j++) blk[j]=malloc(2*(X+1));
  int i=S;
  while(i>=1){
    int base=((i-1)/B)*B;            // layers base..i-1 needed; base is a checkpoint index multiple of B
    memcpy(blk[0],ck[base/B],2*(X+1));
    for(int j=1;j<i-base;j++) step(blk[j-1],blk[j]);   // blk[j] = layer base+j
    for(;i>base;i--){
      uint16_t *prev=blk[i-1-base];
      uint64_t t=3*x+1; int found=0;
      for(int k=0;;k++){ if(k>0){ if(t&1) break; t>>=1; } uint64_t xp=t; if(xp>X) continue;
        if(prev[xp]!=0xFFFF && prev[xp]+k==K){ ks[i]=k; K-=k; x=xp; if(x>vmax) vmax=x; found=1; break; } }
      if(!found){ printf("RECON FAILED at layer %d\n",i); return 1; }
    }
  }
  if(x!=h){ printf("RECON FAILED (end %llu)\n",(unsigned long long)x); return 1; }
  printf("s=%d K0=%d minK=%d ratio=%.6f %s max=%llu k:",S,k0,mk,pow(2.0,mk-S*log2(3.0)),mk==k0?"=K0":"!=K0",(unsigned long long)vmax);
  for(int j=1;j<=S;j++) printf("%s%d",j>1?",":"",ks[j]);
  printf("\n");
  return 0;
}
