// Wave 23, P48: an Apollonian-style reciprocity test for the F1 family.
// F1 values: N = g(h1) g(h2) + 2^s (2h1+1)(2h2+1) - [s=0],  g = 2 psi + 1,
//   psi(0)=0, psi(2h)=psi(h)+h, psi(2h+1)=psi(h)+floor(h/2)+1   (HYP-9657 / THM-4885 coordinates).
// r(N) = number of triples (h1,h2,s) (saturating at 255).
// Test: miss rate of F1 on odd squares m^2 (and t*m^2) in [2^lo, 2^K) against locally matched controls
// N + j*M, M = 8*3*5*7*11*13 = 120120, |j| <= J (same residues mod 8,3,5,7,11,13).
// Usage: f1sq K lo J
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
int main(int argc,char**argv){
  int K=argc>1?atoi(argv[1]):27, lo=argc>2?atoi(argv[2]):22, J=argc>3?atoi(argv[3]):20;
  i64 X=1LL<<K;
  i64 Hmax=X/2+2;
  i64 *psi=malloc(Hmax*sizeof(i64));
  psi[0]=0;
  for(i64 h=1;h<Hmax;h++){ i64 q=h>>1; psi[h]= (h&1)? psi[q]+(q>>1)+1 : psi[q]+q; }
  uint8_t *r=calloc(X+1,1);
  u64 triples=0;
  for(int s=0;s<K;s++){
    for(i64 h1=0;h1<Hmax;h1++){
      i64 y1=(2*h1+1)<<s; if(y1>X) break;
      i64 g1=2*psi[h1]+1;
      i64 lim=X/y1; i64 h2max=(lim-1)/2; if(h2max<0) break;
      for(i64 h2=0;h2<=h2max;h2++){
        i64 N=g1*(2*psi[h2]+1)+y1*(2*h2+1)-(s==0);
        if(N<=X){ if(r[N]<255) r[N]++; triples++; }
      }
    }
  }
  i64 L=1LL<<lo;
  u64 nall=0,mall=0;
  for(i64 N=L|1;N<X;N+=2){ nall++; if(!r[N]) mall++; }
  printf("F1 census: K=%d, triples=%llu; odd N in [2^%d,2^%d): n=%llu misses=%llu rate=%.5f\n",K,triples,lo,K,nall,mall,(double)mall/nall);
  const i64 M=120120;
  int ts[]={1,3,5,7,11,15,21,35,105};
  for(int ti=0;ti<9;ti++){
    i64 t=ts[ti];
    u64 nsq=0,msq=0,nc=0,mc=0; double rs=0,rc=0;
    for(i64 m=1;;m+=2){
      i64 N=t*m*m; if(N>=X) break; if(N<L) continue;
      if(t>1){ // t*m^2 with t squarefree: require gcd(m, t)=1 for a clean square class
        i64 a=m,b=t; while(b){i64 c=a%b;a=b;b=c;} if(a!=1) continue; }
      nsq++; rs+=r[N]; if(!r[N]) msq++;
      for(int j=-J;j<=J;j++){ if(!j) continue; i64 C=N+j*M; if(C<L||C>=X) continue;
        // exclude controls that are themselves in the square class t*k^2
        i64 q=C/t; if(q*t==C){ i64 k=(i64)(__builtin_sqrt((double)q)+0.5); if(k*k==q) continue; }
        nc++; rc+=r[C]; if(!r[C]) mc++; }
    }
    double ps=nsq?(double)msq/nsq:0, pc=nc?(double)mc/nc:0;
    printf("t=%3lld: class t*m^2: n=%7llu misses=%6llu rate=%.5f mean r=%.3f | controls: n=%9llu misses=%7llu rate=%.5f mean r=%.3f | miss ratio %.3f  mean ratio %.3f\n",
      t,nsq,msq,ps,nsq?rs/nsq:0,nc,mc,pc,nc?rc/nc:0, pc>0?ps/pc:-1.0, (nc&&nsq)?(rs/nsq)/(rc/nc):-1.0);
  }
  // small squares listing
  printf("odd squares below 2^14 missed by F1:");
  for(i64 m=1;m*m<(1<<14);m+=2) if(!r[m*m]) printf(" %lld",m*m);
  printf("\n");
  return 0;
}
