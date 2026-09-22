// Deterministic analogues: (1) Collatz T on classes mod 2^m: no multiplicative descent within m halvings.
// (2) greedy G (3-adic) on classes mod 3^{J+1}: no descent within J steps.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
int main(){
  const double L2=log(2.0), L3=log(3.0);
  printf("Collatz T (shortcut n->n/2, (3n+1)/2): classes mod 2^m never descending within m steps\n");
  for(int m=1;m<=26;m++){
    uint64_t N=1ULL<<m, bad=0;
    for(uint64_t c=0;c<N;c++){
      // simulate parity of class c for m steps using exact integer c (valid for whole class)
      unsigned __int128 x=c; double cost=0; int desc=0;
      for(int s=0;s<m;s++){
        if(x&1){ x=(3*x+1)/2; cost+=L3-L2; } else { x/=2; cost-=L2; }
        if(cost<0){desc=1;break;}
      }
      if(!desc) bad++;
    }
    printf("m=%2d bad=%llu frac=%.4g  log2(bad)/m=%.4f\n",m,(unsigned long long)bad,(double)bad/N,log2((double)bad)/m);
  }
  printf("Greedy G: classes mod 3^{J+1} (units) with no descent within J steps\n");
  uint64_t P[30]; P[0]=1; for(int i=1;i<30;i++) P[i]=P[i-1]*3;
  for(int J=1;J<=15;J++){
    uint64_t M=P[J+1], bad=0, units=0;
    for(uint64_t c=1;c<M;c++){ if(c%3==0) continue; units++;
      unsigned __int128 x=c + M; // representative > small to avoid 1-cycle degeneracy? use c+M
      double cost=0; int desc=0;
      for(int s=0;s<J;s++){
        int k=0; unsigned __int128 t=x; // minimal k with 2^k x = 4 or 7 mod 9
        while(1){ uint64_t r=(uint64_t)(t%9); if(r==4||r==7) break; t*=2; k++; }
        x=(t-1)/3; cost+=k*L2-L3; if(cost<0){desc=1;break;}
      }
      if(!desc) bad++;
    }
    printf("J=%2d bad=%llu frac=%.4g  log3(bad)/J=%.4f\n",J,(unsigned long long)bad,(double)bad/units,log((double)bad)/log(3.0)/J);
  }
  return 0;
}
