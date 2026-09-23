// Forward E-game (Q1) exceptional threads: exact integer formulation, low-precision oracle table,
// and a transposition table of exact bounds on f_p(y) that persists across levels.
//
// Definitions (see collatz_procgen_20260922_dim_forward.c for the plain DFS version):
//   f_p(y) = min #(x3-moves) over E-paths from y (known mod 2^p) making exactly p halvings.
//   x mod 2^m is exceptional iff f_b(x mod 2^b) >= A_b + 1 for all b <= m, A_b = floor(b log_3 2).
//   f_b(x mod 2^b) is nondecreasing in b.  If A_m == A_{m-1}, every lift of a bad class is bad.
// Search: search(p,y,r) for even y, p > P0, returns v with
//     v <= r  ->  f_p(y) <= v   (a path of cost v was found)
//     v >  r  ->  f_p(y) >= v   (exact lower bound; IDA*-style minimum over the pruned frontier)
//   lower bounds used:  f_p(y) >= f_P0(y mod 2^P0)  (prefix argument; table)  and TT entries.
//   Every bad class carries a lower bound L(x) <= f_m(x); a lift is tested only if L <= A_m.
// usage: dim_forward_tt P0 MMAX log2(TT entries) [dumpdir]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef W128
typedef unsigned __int128 word;
#define WBITS 128
#else
typedef uint64_t word;
#define WBITS 64
#endif
typedef struct { word y; uint8_t p, lo, hi, used; } Entry;
static int P0;
static uint8_t *T;               // T[i] = f_P0(2i)
static uint64_t MASK0;
static long long nodes, tthits, ttstores;
static int A[400];
static Entry *TT; static uint64_t TTMASK; static int TTBITS; static int TTMIN;
static inline word MASKP(int p){ return p>=WBITS ? ~(word)0 : (((word)1)<<p)-1; }
static inline int fP0(uint64_t y){ y &= MASK0; if(y&1) return 1 + T[((3*y+1)&MASK0)>>1]; return T[y>>1]; }
static inline uint64_t hsh(word y,int p){
  uint64_t lo=(uint64_t)y, hi=(uint64_t)(y>>(WBITS>64?64:0)); if(WBITS==64) hi=0;
  uint64_t h = lo*0x9E3779B97F4A7C15ULL ^ (hi+0x632BE59BD9B4E019ULL)*0xC2B2AE3D27D4EB4FULL ^ (uint64_t)p*0x94D049BB133111EBULL;
  h ^= h>>31; return (h*0xBF58476D1CE4E5B9ULL)>>(64-TTBITS);
}
#define BUCKET 4
static inline Entry* tt_find(word y,int p){
  uint64_t h=hsh(y,p)&~(uint64_t)(BUCKET-1);
  for(int i=0;i<BUCKET;i++){ Entry*e=&TT[h+i]; if(e->used && e->p==p && e->y==y) return e; }
  return NULL;
}
static inline void tt_store(word y,int p,int lo,int hi){
  uint64_t h=hsh(y,p)&~(uint64_t)(BUCKET-1); Entry*victim=NULL;
  for(int i=0;i<BUCKET;i++){ Entry*e=&TT[h+i];
    if(e->used && e->p==p && e->y==y){ if(lo>e->lo) e->lo=lo; if(hi<e->hi) e->hi=hi; return; }
    if(!e->used){ if(!victim||victim->used) victim=e; }
    else if(!victim || (victim->used && e->p < victim->p)) victim=e; }
  victim->y=y; victim->p=p; victim->lo=lo; victim->hi=hi; victim->used=1; ttstores++;
}
static int search(word y,int p,int r);
// cost to finish from state z (any parity) at precision q (q >= P0), budget r; same return convention
static inline int child(word z,int q,int r){
  if(q==P0) return fP0((uint64_t)z);
  if(z&1) return 1 + search((3*z+1)&MASKP(q), q, r-1);
  return search(z,q,r);
}
static int search(word y,int p,int r){        // y even, p > P0
  nodes++;
  int lb = T[((uint64_t)y&MASK0)>>1];
  if(lb>r) return lb;
  int use = (p>=TTMIN);
  if(use){ Entry*e=tt_find(y,p);
    if(e){ tthits++; if(e->lo>r) return e->lo; if(e->hi<=r) return e->hi; if(e->lo>lb) lb=e->lo; } }
  int v1 = child(y>>1,p-1,r);
  if(v1<=r){ if(use) tt_store(y,p,lb,v1); return v1; }
  int v2 = 2 + search((9*y+4)&MASKP(p),p,r-2);
  if(v2<=r){ if(use) tt_store(y,p,lb,v2); return v2; }
  int v = v1<v2?v1:v2; if(v<lb) v=lb;           // failure: exact lower bound > r
  if(use) tt_store(y,p,v,255); return v;
}
static void print_word(FILE*f, word x){
  char buf[64]; int n=0; if(x==0){ fputc('0',f); return; }
  while(x){ buf[n++]='0'+(int)(x%10); x/=10; } while(n) fputc(buf[--n],f);
}
int main(int argc,char**argv){
  if(argc<4){ fprintf(stderr,"usage: %s P0 MMAX ttbits [dumpdir]\n",argv[0]); return 1; }
  P0=atoi(argv[1]); int MMAX=atoi(argv[2]); TTBITS=atoi(argv[3]); const char*dump=argc>4?argv[4]:NULL;
  if(P0<2||P0>31||MMAX>WBITS||MMAX<P0||TTBITS<10||TTBITS>30){ fprintf(stderr,"bad args\n"); return 1; }
  const double l32=log(2.0)/log(3.0);
  for(int b=0;b<400;b++){ double v=b*l32; A[b]=(int)floor(v); if(b>0 && fabs(v-floor(v+0.5))<1e-9){fprintf(stderr,"precision\n");return 1;} }
  clock_t t0=clock();
  uint8_t *prev=NULL, *cur=NULL;
  size_t cap=1<<22; word *bad=malloc(sizeof(word)*cap), *nb=malloc(sizeof(word)*cap);
  uint8_t *L=malloc(cap), *nL=malloc(cap); size_t nbad=0;
  cur=malloc(1); cur[0]=0; bad[0]=1; L[0]=1; nbad=1;
  printf("# forward E-game (TT): m A_m beatty bad tested min_nontriv_rep min_abs_signed nodes cpu_s\n");
  printf("m=1 A=0 beatty=1 bad=1 tested=0 minrep=- minabs=- nodes=0 t=0.0\n");
  for(int p=2;p<=P0;p++){
    prev=cur; uint64_t n=1ULL<<(p-1); cur=malloc(n);
    uint64_t mp=(1ULL<<p)-1, mq=(1ULL<<(p-1))-1;
    for(uint64_t i=0;i<n;i++){
      uint64_t y=2*i; int best=255;
      for(int j=0;2*j<best;j++){
        uint64_t z=y>>1; int v;
        if(z&1) v=1+prev[((3*z+1)&mq)>>1]; else v=prev[z>>1];
        v+=2*j; if(v<best) best=v; y=(9*y+4)&mp;
      }
      cur[i]=(uint8_t)best;
    }
    free(prev);
    size_t n2=0; uint64_t half=1ULL<<(p-1);
    for(size_t i=0;i<nbad;i++) for(int t=0;t<2;t++){
      uint64_t x=(uint64_t)bad[i]+(t?half:0);
      int f = 1 + cur[((3*x+1)&mp)>>1];
      if(f>=A[p]+1){ nb[n2]=x; nL[n2]=f; n2++; }
    }
    word*tmp=bad; bad=nb; nb=tmp; uint8_t*tl=L; L=nL; nL=tl; nbad=n2;
    uint64_t minrep=~0ULL, minabs=~0ULL; uint64_t full=(1ULL<<p);
    for(size_t i=0;i<nbad;i++){ uint64_t x=(uint64_t)bad[i]; if(x==full-1) continue; if(x<minrep) minrep=x; uint64_t ab = x < full/2 ? x : full-x; if(ab<minabs) minabs=ab; }
    printf("m=%d A=%d beatty=%d bad=%zu tested=%zu minrep=%llu minabs=%llu nodes=0 t=%.1f\n",p,A[p],A[p]>A[p-1],nbad,2*nbad,
      (unsigned long long)minrep,(unsigned long long)minabs,(double)(clock()-t0)/CLOCKS_PER_SEC);
    fflush(stdout);
  }
  T=cur; MASK0=(1ULL<<P0)-1;
  TTMIN=P0+4; if(getenv("TTMIN")) TTMIN=atoi(getenv("TTMIN"));
  TT=calloc((size_t)1<<TTBITS,sizeof(Entry)); TTMASK=((uint64_t)1<<TTBITS)-1;
  if(!TT){ fprintf(stderr,"TT alloc failed\n"); return 1; }
  for(int m=P0+1;m<=MMAX;m++){
    size_t n2=0, ntest=0; nodes=0; tthits=0; ttstores=0; word half=((word)1)<<(m-1);
    if(2*nbad>cap){ cap=4*nbad; nb=realloc(nb,sizeof(word)*cap); bad=realloc(bad,sizeof(word)*cap); nL=realloc(nL,cap); L=realloc(L,cap); }
    int Bm=A[m];
    for(size_t i=0;i<nbad;i++) for(int t=0;t<2;t++){
      word x=bad[i]+(t?half:0); int lb=L[i];
      if(lb<=Bm){ ntest++;
        int v = 1 + search((3*x+1)&MASKP(m), m, Bm-1);   // x odd: forced first move
        if(v<=Bm) continue;                                 // descends: good
        lb=v; }
      nb[n2]=x; nL[n2]=(uint8_t)(lb>255?255:lb); n2++;
    }
    word*tmp=bad; bad=nb; nb=tmp; uint8_t*tl=L; L=nL; nL=tl; nbad=n2;
    word minus1 = MASKP(m);
    word minrep=~(word)0, minabs=~(word)0;
    for(size_t i=0;i<nbad;i++){ word x=bad[i]; if(x==minus1) continue; if(x<minrep) minrep=x; word ab = (x < half) ? x : (minus1-x)+1; if(ab<minabs) minabs=ab; }
    int slack[8]={0}; for(size_t i=0;i<nbad;i++){ int s=L[i]-A[m]-1; if(s<0) s=0; if(s>7) s=7; slack[s]++; }
    printf("m=%d A=%d beatty=%d bad=%zu tested=%zu minrep=",m,A[m],A[m]>A[m-1],nbad,ntest); print_word(stdout,minrep);
    printf(" minabs="); print_word(stdout,minabs);
    printf(" nodes=%lld tthits=%lld ttstores=%lld slackLB=[%d %d %d %d %d] t=%.1f\n",nodes,tthits,ttstores,slack[0],slack[1],slack[2],slack[3],slack[4],(double)(clock()-t0)/CLOCKS_PER_SEC); fflush(stdout);
    if(dump){ char fn[512]; snprintf(fn,sizeof fn,"%s/fwdtt_bad_m%d.txt",dump,m); FILE*g=fopen(fn,"w");
      for(size_t i=0;i<nbad;i++){ print_word(g,bad[i]); fprintf(g," %d\n",L[i]);} fclose(g); }
  }
  return 0;
}
