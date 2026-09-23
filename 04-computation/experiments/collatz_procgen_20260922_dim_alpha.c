// alpha_h(i) = f_i(h mod 2^i): minimum number of x3-moves an E-path from the 2-adic point h has made
// at its i-th halving (exact).  Used by the perturbation lemma (Theorem P in
// collatz_procgen_20260922_exceptional_dimension.md): h - 2^i c/3^j is hostile when h is hostile
// (with negative orbit), j <= alpha_h(i), j >= j_h, and (1-3^-j)/2 > (|h|-1) + 2^i c/3^j.
// usage: dim_alpha P0 IMAX num den [num den ...]   (h = num/den, den odd; prints alpha_h(i), i=1..IMAX)
// For i <= P0 the value is a table lookup; for i > P0 an exact DFS with the table as tail/lower bound,
// run with budgets B = lower bound, lower bound+1, ... until a path is found.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128;
static int P0, B; static uint8_t *T; static uint64_t MASK0;
static uint8_t **TAB;   // TAB[p][i] = f_p(2i) for p <= P0
static inline u128 MASKP(int p){ return p>=128 ? ~(u128)0 : (((u128)1)<<p)-1; }
static inline int fP(int p, uint64_t y){ // f_p(y mod 2^p), p <= P0
  if(p==0) return 0; uint64_t m=(p>=64)?~0ULL:((1ULL<<p)-1); y&=m;
  if(y&1) return 1 + TAB[p][((3*y+1)&m)>>1]; return TAB[p][y>>1];
}
static inline int fP0(uint64_t y){ return fP(P0,y); }
static int dfs(u128 y,int p,int a){               // y mod 2^p, p > P0
  u128 mk=MASKP(p);
  if(y&1){ y=(3*y+1)&mk; a++; }
  for(;;){
    if(a + T[((uint64_t)y & MASK0)>>1] > B) return 0;
    u128 z=y>>1;
    if(p-1==P0){ if(a + fP0((uint64_t)z) <= B) return 1; }
    else if(dfs(z,p-1,a)) return 1;
    y=(9*y+4)&mk; a+=2;
  }
}
static u128 inv_mod2k(u128 d){ u128 x=1; for(int i=0;i<7;i++) x*=2-d*x; return x; } // d odd, inverse mod 2^128
int main(int argc,char**argv){
  P0=atoi(argv[1]); int IMAX=atoi(argv[2]);
  TAB=calloc(P0+1,sizeof(uint8_t*));
  TAB[1]=malloc(1); TAB[1][0]=0;
  for(int p=2;p<=P0;p++){
    uint64_t n=1ULL<<(p-1); TAB[p]=malloc(n); uint64_t mp=(1ULL<<p)-1;
    for(uint64_t i=0;i<n;i++){ uint64_t y=2*i; int best=255;
      for(int j=0;2*j<best;j++){ int v=fP(p-1,y>>1)+2*j; if(v<best) best=v; y=(9*y+4)&mp; }
      TAB[p][i]=best; }
  }
  T=TAB[P0]; MASK0=(1ULL<<P0)-1;
  for(int arg=3; arg+1<argc; arg+=2){
    long long num=atoll(argv[arg]), den=atoll(argv[arg+1]);
    u128 d=(u128)den; u128 inv=inv_mod2k(d);
    u128 h=((u128)(__int128)num)*inv;             // h mod 2^128
    printf("h=%lld/%lld alpha:",num,den);
    int prev=0;
    for(int i=1;i<=IMAX;i++){
      int val;
      if(i<=P0) val=fP(i,(uint64_t)h);
      else {
        u128 y=h&MASKP(i); int lb=prev>fP0((uint64_t)y)?prev:fP0((uint64_t)y);
        for(B=lb;;B++){ if(dfs(y,i,0)) break; }
        val=B;
      }
      printf(" %d",val); prev=val; fflush(stdout);
    }
    printf("\n");
  }
  return 0;
}
