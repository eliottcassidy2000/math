// collatz_procgen_20260922_loops_hubs.c
// Minimal cycles through a hub h in the E-graph (reverse-move DP from h back to h, values <= XMAX).
// A cycle through h with q multiplications and p halvings has ratio 2^p/3^q = prod over its multiplication
// points v of (1 + 1/(3v)) > 1.  Inserting it into any loop through 1 that passes through h changes
// (s,K) by (q,p) and multiplies the loop ratio by 2^p/3^q.
// Usage: loops_hubs h QMAX XMAX [recon]   prints for each q the minimal p (and, with recon, one cycle)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
int main(int argc,char**argv){
  uint64_t h=strtoull(argv[1],0,10); int QMAX=atoi(argv[2]); uint64_t XMAX=strtoull(argv[3],0,10); int recon=argc>4;
  uint16_t **L=malloc(sizeof(uint16_t*)*(QMAX+1));
  int nl = recon? QMAX+1 : 2;
  for(int i=0;i<nl;i++){ L[i]=malloc(2*(XMAX+1)); if(!L[i]){fprintf(stderr,"oom\n");return 1;} }
  uint16_t *cur=L[0],*nxt; memset(cur,255,2*(XMAX+1)); cur[h]=0;
  for(int q=1;q<=QMAX;q++){
    nxt = recon? L[q] : L[q&1]; memset(nxt,255,2*(XMAX+1));
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
    if(nxt[h]!=0xFFFF){
      int pp=nxt[h]; long double eps=pp-q*log2l(3.0L);
      printf("h=%llu q=%3d p=%3d ratio=%.6Lf eps=p-q*log2(3)=%.6Lf",(unsigned long long)h,q,pp,powl(2.0L,eps),eps);
      if(recon){
        int K=pp; uint64_t x=h; int ok=1; int *ks=malloc(sizeof(int)*(q+1));
        for(int i=q;i>=1;i--){ uint64_t t=3*x+1; int found=0;
          for(int k=0;;k++){ if(k>0){ if(t&1) break; t>>=1; } uint64_t xp=t; if(xp>XMAX) continue;
            if(L[i-1][xp]!=0xFFFF && L[i-1][xp]+k==K){ ks[i]=k; K-=k; x=xp; found=1; break; } }
          if(!found){ok=0;break;} }
        if(ok && x==h){ printf(" k:"); for(int i=1;i<=q;i++) printf("%s%d",i>1?",":"",ks[i]); }
        else printf(" RECON FAILED");
        free(ks);
      }
      printf("\n");
    }
    fflush(stdout);
    cur=nxt;
  }
  return 0;
}
