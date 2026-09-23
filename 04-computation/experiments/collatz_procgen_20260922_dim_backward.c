// Backward E-game (Q2) exceptional threads on Z_3 units, exact integer formulation + oracle table.
//
// Move from a 3-adic unit x: k>=0 with 2^k x = 1 mod 3, y=(2^k x-1)/3, legal iff 3 !| y; cost k log2 - log3.
// A path of s moves has multiplier 2^K/3^s (K = sum k); it descends iff K <= C_s := floor(s log_2 3).
//   g_s(x) = min K over legal s-move paths from x (decidable from x mod 3^(s+1)).
//   g_s(x mod 3^(s+1)) is nondecreasing in s (prefix argument), so
//   Bad_r = { x mod 3^(r+1) : x mod 3^r in Bad_{r-1} and g_r(x) >= C_r + 1 }.
// Recursion: g_0 = 0;  g_s(x) = min_{legal k} k + g_{s-1}(y_k mod 3^s),  y_{k+2} = 4 y_k + 1.
// Phase 1: full DP table g_P0 on units mod 3^(P0+1) (uint8), and Bad_r for r <= P0.
// Phase 2: for r > P0 test the three lifts of each bad class by DFS over the first r-P0 moves:
//   exact tail  g_P0(y) at precision 3^(P0+1) (table),
//   lower bound g_n(y) >= g_P0(y mod 3^(P0+1)) for n >= P0 (prefix argument) -> prune.
//   Optional transposition table (ttbits>0) of exact bounds on g, persistent across levels.
// usage: dim_backward P0 RMAX [ttbits] [dumpdir]     (RMAX <= 39 unless compiled -DW128)
// fast build: -DP0C=17 -DP0C_M0=387420489ULL  (M0 = 3^(P0+1))
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef W128
typedef unsigned __int128 word;
#define NMAX 78
#else
typedef uint64_t word;
#define NMAX 40
#endif
typedef struct { word y; uint8_t n, lo, hi, used; } Entry;
#ifdef P0C
#define P0 P0C                 // compile-time oracle depth (faster modular reductions)
#else
static int P0;                 // table covers g_P0 on units mod 3^(P0+1)
#endif
static uint8_t *G;             // G[idx(x)] = g_P0(x), idx(x)=2*(x/3)+(x%3)-1
static word P3[NMAX+2];
static long long nodes, tthits;
static int C[400];
static Entry *TT=NULL; static int TTBITS=0;
static inline size_t idx(uint64_t x){ return 2*(x/3)+(x%3)-1; }
#ifdef P0C
static const uint64_t M0 = (uint64_t)P0C_M0;   // 3^(P0+1), compile-time
#else
static uint64_t M0;
#endif
static inline uint64_t mod0(word y){
#ifdef W128
  uint64_t hi=(uint64_t)(y>>64), lo=(uint64_t)y;
  if(!hi) return lo % M0;
  static uint64_t two64=0; if(!two64){ two64 = (uint64_t)((((unsigned __int128)1)<<64) % M0); }
  return (uint64_t)(((unsigned __int128)(hi % M0) * two64 + lo % M0) % M0);
#else
  return y % M0;
#endif
}
static inline int gP0(word y){ return G[idx(mod0(y))]; }
static inline uint64_t hsh(word y,int n){
  uint64_t lo=(uint64_t)y, hi=0;
#ifdef W128
  hi=(uint64_t)(y>>64);
#endif
  uint64_t h = lo*0x9E3779B97F4A7C15ULL ^ (hi+0x632BE59BD9B4E019ULL)*0xC2B2AE3D27D4EB4FULL ^ (uint64_t)n*0x94D049BB133111EBULL;
  h ^= h>>31; return (h*0xBF58476D1CE4E5B9ULL)>>(64-TTBITS);
}
#define BUCKET 4
static inline Entry* tt_find(word y,int n){
  uint64_t h=hsh(y,n)&~(uint64_t)(BUCKET-1);
  for(int i=0;i<BUCKET;i++){ Entry*e=&TT[h+i]; if(e->used && e->n==n && e->y==y) return e; }
  return NULL;
}
static inline void tt_store(word y,int n,int lo,int hi){
  uint64_t h=hsh(y,n)&~(uint64_t)(BUCKET-1); Entry*victim=NULL;
  if(lo>255) lo=255; if(hi>255) hi=255;
  for(int i=0;i<BUCKET;i++){ Entry*e=&TT[h+i];
    if(e->used && e->n==n && e->y==y){ if(lo>e->lo) e->lo=lo; if(hi<e->hi) e->hi=hi; return; }
    if(!e->used){ if(!victim||victim->used) victim=e; }
    else if(!victim || (victim->used && e->n < victim->n)) victim=e; }
  victim->y=y; victim->n=n; victim->lo=lo; victim->hi=hi; victim->used=1;
}
static int TTMIN;   // use TT only at precision n >= TTMIN
// y unit known mod 3^n (n > P0+1), must make exactly n-1 more moves.  Returns v:
//   v <= r -> g_{n-1}(y) <= v (path found);  v > r -> g_{n-1}(y) >= v (exact lower bound).
static int search(word y,int n,int r){
  nodes++;
  int lb=gP0(y);
  if(lb>r) return lb;
  Entry*e=NULL;
  if(TT && n>=TTMIN){ e=tt_find(y,n); if(e){ tthits++; if(e->lo>r) return e->lo; if(e->hi<=r) return e->hi; if(e->lo>lb) lb=e->lo; } }
  word M=P3[n-1];                       // child modulus 3^(n-1)
  int k=((uint64_t)(y%3)==1)?0:1;
  word yk = (k==0) ? (y-1)/3 : (2*y-1)/3;  yk %= M;   // y_k mod 3^(n-1)
  int best=1000;                         // lower bound over the not-yet-found region
  for(;k<=r;k+=2){
    if((uint64_t)(yk%3)!=0){
      int v;
      if(n-1==P0+1) v = k + gP0(yk);
      else v = k + search(yk,n-1,r-k);
      if(v<=r){ if(TT && n>=TTMIN) tt_store(y,n,lb,v); return v; }
      if(v<best) best=v;
    }
    { word t=4*yk+1; while(t>=M) t-=M; yk=t; }
  }
  if(k<best) best=k;                     // unexplored k' >= k cost at least k' > r
  if(best<lb) best=lb;
  if(TT && n>=TTMIN) tt_store(y,n,best,255);
  return best;
}
static void print_word(FILE*f, word x){
  char buf[64]; int n=0; if(x==0){ fputc('0',f); return; }
  while(x){ buf[n++]='0'+(int)(x%10); x/=10; } while(n) fputc(buf[--n],f);
}
int main(int argc,char**argv){
  if(argc<3){ fprintf(stderr,"usage: %s P0 RMAX [ttbits] [dumpdir]\n",argv[0]); return 1; }
#ifdef P0C
  if(atoi(argv[1])!=P0C){ fprintf(stderr,"compiled for P0=%d\n",P0C); return 1; }
#else
  P0=atoi(argv[1]);
#endif
  int RMAX=atoi(argv[2]); TTBITS=argc>3?atoi(argv[3]):0; const char*dump=argc>4?argv[4]:NULL;
  TTMIN = P0+3; if(getenv("TTMIN")) TTMIN=atoi(getenv("TTMIN"));
  P3[0]=1; for(int i=1;i<=NMAX+1;i++) P3[i]=P3[i-1]*3;
#ifndef P0C
  M0=(uint64_t)P3[P0+1];
#endif
  if(P0<1||P0>19||RMAX>=NMAX||RMAX<P0){ fprintf(stderr,"bad args\n"); return 1; }
  const double l23=log(3.0)/log(2.0);
  for(int s=0;s<400;s++){ double v=s*l23; C[s]=(int)floor(v); if(s>0 && fabs(v-floor(v+0.5))<1e-9){fprintf(stderr,"precision\n");return 1;} }
  clock_t t0=clock();
  // ---------- phase 1 ----------
  uint8_t *prev=malloc(2), *cur=NULL; prev[0]=0; prev[1]=0;       // g_0 on units mod 3
  size_t cap=1<<22; word *bad=malloc(sizeof(word)*cap), *nb=malloc(sizeof(word)*cap);
  uint8_t *L=malloc(cap), *nL=malloc(cap); size_t nbad=2; bad[0]=1; bad[1]=2; L[0]=0; L[1]=0;
  printf("# backward E-game: r C_r bad tested min_nontriv_rep nodes cpu_s   (classes mod 3^(r+1); 'nontriv' = not the class of 1)\n");
  for(int s=1;s<=P0;s++){
    uint64_t Mn=(uint64_t)P3[s+1], Mc=(uint64_t)P3[s];     // x mod 3^(s+1), child mod 3^s
    size_t nu=2*(size_t)(Mn/3); cur=malloc(nu);
    for(uint64_t x=1;x<Mn;x++){
      if(x%3==0) continue;
      int k=(x%3==1)?0:1; uint64_t yk=((k==0)?(x-1)/3:(2*x-1)/3)%Mc; int best=255;
      for(;k<best;k+=2){
        if(yk%3){ int v=k+prev[idx(yk)]; if(v<best) best=v; }
        yk=(4*yk+1)%Mc;
      }
      cur[idx(x)]=(uint8_t)best;
    }
    free(prev); prev=cur;
    size_t n2=0;
    for(size_t i=0;i<nbad;i++) for(int t=0;t<3;t++){
      uint64_t x=(uint64_t)bad[i]+t*Mc;
      int g=cur[idx(x)];
      if(g>=C[s]+1){ nb[n2]=x; nL[n2]=g; n2++; }
    }
    word*tmp=bad; bad=nb; nb=tmp; uint8_t*tl=L; L=nL; nL=tl; nbad=n2;
    uint64_t minrep=~0ULL; for(size_t i=0;i<nbad;i++){ uint64_t x=(uint64_t)bad[i]; if(x!=1 && x<minrep) minrep=x; }
    printf("r=%d C=%d bad=%zu tested=%zu minrep=%llu nodes=0 t=%.1f\n",s,C[s],nbad,3*nbad,(unsigned long long)minrep,(double)(clock()-t0)/CLOCKS_PER_SEC);
    fflush(stdout);
  }
  G=cur;
  if(TTBITS>0){ TT=calloc((size_t)1<<TTBITS,sizeof(Entry)); if(!TT){fprintf(stderr,"TT alloc\n");return 1;} }
  // ---------- phase 2 ----------
  for(int r=P0+1;r<=RMAX;r++){
    size_t n2=0, ntest=0; nodes=0; tthits=0; word Mold=P3[r];
    if(3*nbad>cap){ cap=4*nbad; nb=realloc(nb,sizeof(word)*cap); bad=realloc(bad,sizeof(word)*cap); nL=realloc(nL,cap); L=realloc(L,cap); }
    for(size_t i=0;i<nbad;i++) for(int t=0;t<3;t++){
      word x=bad[i]+(word)t*Mold; int lb=L[i];
      if(lb<=C[r]){ ntest++;
        int v=search(x,r+1,C[r]);
        if(v<=C[r]) continue;
        lb=v; }
      nb[n2]=x; nL[n2]=(uint8_t)(lb>255?255:lb); n2++;
    }
    word*tmp=bad; bad=nb; nb=tmp; uint8_t*tl=L; L=nL; nL=tl; nbad=n2;
    word minrep=~(word)0; for(size_t i=0;i<nbad;i++){ word x=bad[i]; if(x!=1 && x<minrep) minrep=x; }
    printf("r=%d C=%d bad=%zu tested=%zu minrep=",r,C[r],nbad,ntest); print_word(stdout,minrep);
    printf(" nodes=%lld tthits=%lld t=%.1f\n",nodes,tthits,(double)(clock()-t0)/CLOCKS_PER_SEC); fflush(stdout);
    if(dump){ char fn[512]; snprintf(fn,sizeof fn,"%s/bwd_bad_r%d.txt",dump,r); FILE*g=fopen(fn,"w");
      for(size_t i=0;i<nbad;i++){ print_word(g,bad[i]); fprintf(g," %d\n",L[i]); } fclose(g); }
  }
  return 0;
}
