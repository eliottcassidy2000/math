// collatz_procgen_20260922_loops_dp.c
// Layered DP for reverse E-loops through 1 with all values <= XMAX, with loop reconstruction.
//   Reverse move: x -> y = (2^k x - 1)/3, k >= 0, legal iff y is a positive integer, 3 !| y.
//   Loop of length s: x_0 = 1 -> x_1 -> ... -> x_s = 1, K = sum k_i  (forward: an E-cycle through 1
//   with s multiplications 3n+1 and K halvings).
// Layer i stores L_i[x] = min K_i over legal i-move paths 1 -> x with all values <= XMAX (0xFFFF = none).
// Usage:
//   loops_dp SMAX XMAX            : print minK(s), K0(s)=floor((s+1)log2 3), ratio, flag
//   loops_dp SMAX XMAX recon [SMIN]: additionally store every layer and print one optimal loop per s
//                                   (list of k_i and the max value on the loop), only for s >= SMIN.
//                                   With XMAX = H(s) (from loops_height) the loop printed has minimal height.
// Memory: 4*(XMAX+1) bytes, or 2*(SMAX+1)*(XMAX+1) bytes in recon mode (16-bit K).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

static int K0(int s){ // largest K with 2^K < 3^(s+1) = floor((s+1) log2 3); long double is exact here
  // (for s+1 <= 10^4 the distance of (s+1)log2(3) to the nearest integer exceeds 1e-6 >> rounding error)
  long double v=(long double)(s+1)*log2l(3.0L);
  int K=(int)floorl(v);
  return K;
}

int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage: %s SMAX XMAX [recon]\n",argv[0]);return 1;}
  int SMAX=atoi(argv[1]); uint64_t XMAX=strtoull(argv[2],0,10); int recon=(argc>3); int SMIN=(argc>4)?atoi(argv[4]):1;
  uint16_t **L=NULL; uint16_t *cur,*nxt;
  if(recon){
    L=malloc(sizeof(uint16_t*)*(SMAX+1));
    for(int i=0;i<=SMAX;i++){ L[i]=malloc(2*(XMAX+1)); if(!L[i]){fprintf(stderr,"oom\n");return 1;} }
    cur=L[0];
  } else { cur=malloc(2*(XMAX+1)); nxt=malloc(2*(XMAX+1)); }
  memset(cur,255,2*(XMAX+1)); cur[1]=0;
  for(int s=1;s<=SMAX;s++){
    if(recon) nxt=L[s];
    memset(nxt,255,2*(XMAX+1));
    for(uint64_t x=1;x<=XMAX;x++){
      if(cur[x]==0xFFFF) continue;
      int k0=(x%3==1)?0:1; unsigned __int128 p=((unsigned __int128)1)<<k0;
      for(int k=k0;k<=100;k+=2){
        unsigned __int128 t=p*x; if(t>(unsigned __int128)3*XMAX+1) break;
        uint64_t y=(uint64_t)((t-1)/3);
        if(y>=1 && y%3!=0){ int K=cur[x]+k; if(K<0xFFFF && K<nxt[y]) nxt[y]=(uint16_t)K; }
        p<<=2;
      }
    }
    int mk=nxt[1];
    int k0=K0(s);
    if(s<SMIN){ if(!recon){ uint16_t *tmp=cur; cur=nxt; nxt=tmp; } else cur=nxt; continue; }
    if(mk==0xFFFF) printf("s=%3d K0=%3d  no loop with values <= %llu\n",s,k0,(unsigned long long)XMAX);
    else {
      double r=pow(2.0,mk)/pow(3.0,s);
      printf("s=%3d K0=%3d minK=%3d ratio=%.6f %s",s,k0,mk,r,(mk==k0)?"=K0":(mk<k0?"<K0":">K0"));
      if(recon){
        // backtrack from (s,1)
        int *ks=malloc(sizeof(int)*(s+1)); uint64_t x=1; int K=mk; uint64_t vmax=1; int ok=1;
        for(int i=s;i>=1;i--){
          // predecessor xp with x=(2^k xp -1)/3, i.e. xp=(3x+1)/2^k
          uint64_t t=3*x+1; int found=0;
          for(int k=0;;k++){
            if(k>0){ if(t&1) break; t>>=1; }
            uint64_t xp=t; if(xp>XMAX) continue;
            if(L[i-1][xp]!=0xFFFF && L[i-1][xp]+k==K){ ks[i]=k; K-=k; x=xp; if(x>vmax) vmax=x; found=1; break; }
          }
          if(!found){ ok=0; break; }
        }
        if(!ok||x!=1) printf("  RECON FAILED");
        else { printf("  max=%llu k:",(unsigned long long)vmax); for(int i=1;i<=s;i++) printf("%s%d",i>1?",":"",ks[i]); }
        free(ks);
      }
      printf("\n");
    }
    fflush(stdout);
    if(!recon){ uint16_t *tmp=cur; cur=nxt; nxt=tmp; } else cur=nxt;
  }
  return 0;
}
