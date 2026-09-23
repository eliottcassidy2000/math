// collatz_procgen_20260923_sweep_rn_cover.c  (HYP-9122, section 3.5 of the sweep note)
// R_j = { sum_{i<j} 3^i 2^{f_i} : f nonincreasing >= 0 } = U_a 2^a (3^{j-1} + R_{j-1}); HYP-9122(s) <=> 2^{K0(s)} - 3^s in R_s.
// Prints, for each level j, the number of uncovered x (x%3!=0) in ratio bins of x/3^j,
// and the max ratio of an uncovered x in [0.5,C) (excluding near-lower-end), plus min ratio of uncovered above 0.6.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef uint64_t u64;
#define GET(b,x) ((b[(x)>>6]>>((x)&63))&1ULL)
#define SET(b,x) (b[(x)>>6]|=(1ULL<<((x)&63)))
int main(int argc,char**argv){
  int n=atoi(argv[1]); double C=argc>2?atof(argv[2]):3.0;
  u64 p3[40]; p3[0]=1; for(int i=1;i<40;i++) p3[i]=p3[i-1]*3;
  u64 X=(u64)(C*(double)p3[n]); u64 W=X/64+2;
  u64 *A=calloc(W,8),*B=calloc(W,8);
  if(!A||!B){printf("alloc fail\n");return 1;}
  SET(A,0);
  for(int j=1;j<=n;j++){
    memset(B,0,W*8);
    u64 sh=p3[j-1];
    for(u64 x=1;x<=X;x++){
      int t = (x>=sh) ? (int)GET(A,x-sh) : 0;
      if(t || ((x&1)==0 && GET(B,x>>1))) SET(B,x);
    }
    u64 *tmp=A;A=B;B=tmp;
    if(j<8) continue;
    // bins of width 0.05 from 0.5 up to min(C, X/3^j)
    double top = (double)X/p3[j]; if(top> 4.0) top=4.0;
    printf("j=%d uncovered per bin(0.1):",j);
    double maxr=0; u64 total=0;
    for(double lam=0.5; lam<top-1e-9; lam+=0.1){
      u64 lo=(u64)ceil(lam*p3[j]), hi=(u64)floor((lam+0.1)*p3[j]); if(hi>X)hi=X;
      u64 unc=0;
      for(u64 x=lo;x<=hi;x++){ if(x%3==0) continue; if(!GET(A,x)){unc++; double r=(double)x/p3[j]; if(r>maxr) maxr=r;} }
      total+=unc;
      printf(" %llu",(unsigned long long)unc);
    }
    printf(" | total %llu maxratio %.4f\n",(unsigned long long)total,maxr); fflush(stdout);
  }
  return 0;
}

