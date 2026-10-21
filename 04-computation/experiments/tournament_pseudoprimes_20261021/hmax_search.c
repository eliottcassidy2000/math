/* hmax_search.c -- steepest-ascent search for tournaments with many Hamiltonian paths, moves = reversal of a
 * cyclic triangle (preserves scores) and single arcs; restarts from the rotational circulant and random regular
 * tournaments.  Reports the best H found.  mac-mini-2026-10-21-fanopsp.  Usage: ./hmax_search n restarts seed */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n; static uint32_t out[32];
static uint64_t *f;
static uint64_t H(void){ int full=(1<<n)-1; memset(f,0,sizeof(uint64_t)*(size_t)(full+1)*n);
  for(int v=0;v<n;v++) f[((size_t)1<<v)*n+v]=1;
  for(int m=1;m<=full;m++){ uint64_t*row=f+(size_t)m*n; for(int v=0;v<n;v++){ uint64_t c=row[v]; if(!c) continue;
      uint32_t cand=out[v]&~m; while(cand){int w=__builtin_ctz(cand); cand&=cand-1; f[((size_t)(m|(1<<w)))*n+w]+=c;} } }
  uint64_t s=0; for(int v=0;v<n;v++) s+=f[(size_t)full*n+v]; return s; }
static void flip(int a,int b){ out[a]^=1u<<b; out[b]^=1u<<a; }
static int arc(int a,int b){ return (out[a]>>b)&1; }
static uint64_t rng=88172645463325252ULL; static uint64_t xr(void){ rng^=rng<<13; rng^=rng>>7; rng^=rng<<17; return rng; }
int main(int argc,char**argv){ n=atoi(argv[1]); int R=atoi(argv[2]); rng^=atoll(argv[3]);
  f=malloc(sizeof(uint64_t)*((size_t)1<<n)*n); uint64_t best=0; uint32_t bestT[32];
  for(int r=0;r<R;r++){
    memset(out,0,sizeof out); int m=(n-1)/2;
    for(int i=0;i<n;i++) for(int k=1;k<=m;k++){ int j=(i+k)%n; out[i]|=1u<<j; } /* rotational */
    if(r>0){ for(int t=0;t<20*n;t++){ int a=xr()%n,b=xr()%n,c=xr()%n; if(a==b||b==c||a==c) continue;
        if(arc(a,b)&&arc(b,c)&&arc(c,a)){ flip(a,b);flip(b,c);flip(c,a);} } }
    uint64_t h=H(); int improved=1;
    while(improved){ improved=0; uint64_t bh=h; int ba=-1,bb=-1,bc=-1;
      for(int a=0;a<n;a++) for(int b=0;b<n;b++) if(arc(a,b)) for(int c=0;c<n;c++) if(c!=a&&c!=b&&arc(b,c)&&arc(c,a)&&a<b&&a<c){
            flip(a,b);flip(b,c);flip(c,a); uint64_t hh=H(); if(hh>bh){bh=hh;ba=a;bb=b;bc=c;} flip(c,a);flip(b,c);flip(a,b); }
      if(ba>=0){ flip(ba,bb);flip(bb,bc);flip(bc,ba); h=bh; improved=1; } }
    if(h>best){ best=h; memcpy(bestT,out,sizeof out); }
    printf("restart %d: local max H=%llu (best %llu)\n",r,(unsigned long long)h,(unsigned long long)best); fflush(stdout); }
  printf("BEST n=%d H=%llu rows:",n,(unsigned long long)best); for(int i=0;i<n;i++) printf(" %x",bestT[i]); printf("\n"); return 0; }
