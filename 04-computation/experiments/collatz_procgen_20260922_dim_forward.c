// Forward E-game (Q1) exceptional threads, exact integer formulation + low-precision oracle.
//
// E-moves on Z_2: odd x -> 3x+1 (forced); even x -> x/2 (consumes one bit) or x -> 3x+1 -> 9x+4.
// Cost only drops at halvings, so a class x mod 2^m has a descending certificate iff
//     exists b <= m :  f_b(x mod 2^b) <= A_b := floor(b log_3 2)      (3^a < 2^b  <=>  a <= A_b)
// where f_b(y) = min number of x3-moves over E-paths from y that make exactly b halvings
// (decidable from y mod 2^b).  f_b(x mod 2^b) is nondecreasing in b (prefix argument), so
//     Bad_m = { x : x mod 2^{m-1} in Bad_{m-1}  and  f_m(x) >= A_m + 1 },
// and if A_m == A_{m-1} every lift of a bad class is bad (no test needed).
//
// Recursion (exact, integer):  f_0 = 0;  odd y: f_p(y) = 1 + f_p(3y+1);
//   even y: f_p(y) = min_{j>=0} 2j + f_{p-1}(y_j / 2),  y_0 = y, y_{j+1} = 9 y_j + 4 (mod 2^p).
// Phase 1: full DP tables f_p for p <= P0 (even entries only, uint8), and Bad_p.
// Phase 2: for m > P0 test lifts by DFS over the top m-P0 bits with
//   exact tail:   state at precision P0 contributes f_P0(y) exactly (table lookup);
//   lower bound:  f_p(y) >= f_P0(y mod 2^P0) for p >= P0 (prefix argument) -> prune.
// Default: one combined search decides both lifts (dfs2); SINGLE_LIFT=1 tests them separately (dfs).
// usage: dim_forward P0 MMAX [dumpdir]     (P0 <= 31; MMAX <= 64 unless compiled -DW128)
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
static int P0, B;               // oracle precision, current budget A_m
static uint8_t *T;              // T[i] = f_P0(2i), i < 2^(P0-1)
static uint64_t MASK0;          // 2^P0 - 1
static long long nodes;
static int A[400];
static inline word MASKP(int p){ return p>=WBITS ? ~(word)0 : (((word)1)<<p)-1; }
static inline int fP0(uint64_t y){ // y mod 2^P0
  y &= MASK0;
  if(y&1) return 1 + T[((3*y+1)&MASK0)>>1];
  return T[y>>1];
}
// y known mod 2^p, p > P0, a = x3-moves so far; returns 1 iff some continuation consumes all p bits
// with total x3-moves <= B.
static int dfs(word y, int p, int a){
  nodes++;
  word mk = MASKP(p);
  if(y&1){ y=(3*y+1)&mk; a++; }
  for(;;){
    if(a + T[((uint64_t)y & MASK0)>>1] > B) return 0;      // f_p(y) >= f_P0(y mod 2^P0)
    word z = y>>1;
    if(p-1==P0){ if(a + fP0((uint64_t)z) <= B) return 1; }
    else if(dfs(z,p-1,a)) return 1;
    y = (9*y+4)&mk; a += 2; nodes++;                          // excursion (stays even)
  }
}
// Combined search for both lifts x, x+2^(m-1) of a parent class (y = state of the lift with top bit 0).
// The two lifts have identical parities and identical pruning lookups at every internal node
// (their states differ by 3^a 2^(p-1), i.e. only in the top bit); at the leaf (precision P0)
// the residues differ exactly in bit P0-1.  Returns the subset of `need` (bit t = lift t) that
// has a continuation with total x3-moves <= B.
static int dfs2(word y, int p, int a, int need){
  nodes++;
  word mk = MASKP(p); int got=0;
  if(y&1){ y=(3*y+1)&mk; a++; }
  for(;;){
    if(a + T[((uint64_t)y & MASK0)>>1] > B) return got;
    word z = y>>1;
    if(p-1==P0){ uint64_t z0=(uint64_t)z&MASK0;
      if((need&1) && a + fP0(z0) <= B) got|=1;
      if((need&2) && a + fP0(z0^(1ULL<<(P0-1))) <= B) got|=2; }
    else got |= dfs2(z,p-1,a,need&~got);
    if(got==need) return got;
    y = (9*y+4)&mk; a += 2; nodes++;
  }
}
static void print_word(FILE*f, word x){
  char buf[64]; int n=0; if(x==0){ fputc('0',f); return; }
  while(x){ buf[n++]='0'+(int)(x%10); x/=10; } while(n) fputc(buf[--n],f);
}
int main(int argc,char**argv){
  if(argc<3){ fprintf(stderr,"usage: %s P0 MMAX [dumpdir]\n",argv[0]); return 1; }
  P0=atoi(argv[1]); int MMAX=atoi(argv[2]); const char*dump=argc>3?argv[3]:NULL;
  int single = getenv("SINGLE_LIFT")!=NULL;   // independent check: test the two lifts separately
  if(P0<2||P0>31||MMAX>WBITS||MMAX<P0){ fprintf(stderr,"bad args\n"); return 1; }
  const double l32=log(2.0)/log(3.0);
  for(int b=0;b<400;b++){ double v=b*l32; A[b]=(int)floor(v); if(b>0 && fabs(v-floor(v+0.5))<1e-9){fprintf(stderr,"precision\n");return 1;} }
  clock_t t0=clock();
  // ---------- phase 1: DP ----------
  uint8_t *prev=NULL, *cur=NULL;
  size_t cap=1<<22; word *bad=malloc(sizeof(word)*cap), *nb=malloc(sizeof(word)*cap); size_t nbad=0;
  // level 1
  cur=malloc(1); cur[0]=0;  // f_1(0)=0 ; f_1(1)=1+f_1(0)=1
  bad[0]=1; nbad=1;
  printf("# forward E-game: m A_m beatty count min_nontriv_rep min_abs_signed_nontriv nodes cpu_s\n");
  printf("m=1 A=0 beatty=1 bad=1 minrep=- minabs=- nodes=0 t=0.0\n");
  for(int p=2;p<=P0;p++){
    prev=cur; uint64_t n=1ULL<<(p-1); cur=malloc(n);
    uint64_t mp=(1ULL<<p)-1, mq=(1ULL<<(p-1))-1;
    for(uint64_t i=0;i<n;i++){
      uint64_t y=2*i; int best=255;
      for(int j=0;2*j<best;j++){
        uint64_t z=y>>1; int v;  // z mod 2^(p-1)
        if(z&1) v=1+prev[((3*z+1)&mq)>>1]; else v=prev[z>>1];
        v+=2*j; if(v<best) best=v;
        y=(9*y+4)&mp;
      }
      cur[i]=(uint8_t)best;
    }
    free(prev);
    // lift bad list
    size_t n2=0; uint64_t half=1ULL<<(p-1);
    for(size_t i=0;i<nbad;i++) for(int t=0;t<2;t++){
      uint64_t x=(uint64_t)bad[i]+(t?half:0);
      int f = 1 + cur[((3*x+1)&mp)>>1];   // x odd
      if(f>=A[p]+1) nb[n2++]=x;
    }
    word*tmp=bad; bad=nb; nb=tmp; nbad=n2;
    uint64_t minrep=~0ULL, minabs=~0ULL; uint64_t full=(1ULL<<p);
    for(size_t i=0;i<nbad;i++){ uint64_t x=(uint64_t)bad[i]; if(x==full-1) continue; if(x<minrep) minrep=x; uint64_t ab = x < full/2 ? x : full-x; if(ab<minabs) minabs=ab; }
    printf("m=%d A=%d beatty=%d bad=%zu minrep=%llu minabs=%llu nodes=0 t=%.1f\n",p,A[p],A[p]>A[p-1],nbad,
      (unsigned long long)minrep,(unsigned long long)minabs,(double)(clock()-t0)/CLOCKS_PER_SEC);
    fflush(stdout);
  }
  T=cur; MASK0=(P0>=64)?~0ULL:((1ULL<<P0)-1);
  // sanity: recompute f_P0 of the bad list via fP0
  if(dump){ char fn[512]; snprintf(fn,sizeof fn,"%s/fwd_bad_m%d.txt",dump,P0); FILE*g=fopen(fn,"w");
    for(size_t i=0;i<nbad;i++){ print_word(g,bad[i]); fprintf(g," %d\n",fP0((uint64_t)bad[i])); } fclose(g); }
  // ---------- phase 2: thread DFS ----------
  for(int m=P0+1;m<=MMAX;m++){
    size_t n2=0; nodes=0; word half=((word)1)<<(m-1);
    if(2*nbad>cap){ cap=4*nbad; nb=realloc(nb,sizeof(word)*cap); bad=realloc(bad,sizeof(word)*cap); }
    B=A[m]; int beatty = A[m]>A[m-1];
    for(size_t i=0;i<nbad;i++){
      if(!beatty){ nb[n2++]=bad[i]; nb[n2++]=bad[i]+half; continue; }
      if(single){ for(int t=0;t<2;t++){ word x=bad[i]+(t?half:0); if(!dfs(x,m,0)) nb[n2++]=x; } continue; }
      int got=dfs2(bad[i],m,0,3);
      if(!(got&1)) nb[n2++]=bad[i];
      if(!(got&2)) nb[n2++]=bad[i]+half;
    }
    word*tmp=bad; bad=nb; nb=tmp; nbad=n2;
    word minus1 = MASKP(m);                     // 2^m - 1, the class of -1
    word minrep=~(word)0, minabs=~(word)0;
    for(size_t i=0;i<nbad;i++){ word x=bad[i]; if(x==minus1) continue; if(x<minrep) minrep=x; word ab = (x < half) ? x : (minus1-x)+1; if(ab<minabs) minabs=ab; }
    printf("m=%d A=%d beatty=%d bad=%zu minrep=",m,A[m],beatty,nbad); print_word(stdout,minrep);
    printf(" minabs="); print_word(stdout,minabs);
    printf(" nodes=%lld t=%.1f\n",nodes,(double)(clock()-t0)/CLOCKS_PER_SEC); fflush(stdout);
    if(dump){ char fn[512]; snprintf(fn,sizeof fn,"%s/fwd_bad_m%d.txt",dump,m); FILE*g=fopen(fn,"w");
      for(size_t i=0;i<nbad;i++){ print_word(g,bad[i]); fputc('\n',g);} fclose(g); }
  }
  return 0;
}
