// Wave 23: F1 misses on primes versus rough composites (no prime factor <= 97), stratified by
// the number of primes p <= 97 dividing N+1 and by N mod 16.  Usage: f1prime K lo
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef long long i64; typedef unsigned long long u64;
int main(int argc,char**argv){
  int K=atoi(argv[1]), lo=atoi(argv[2]); i64 X=1LL<<K, L=1LL<<lo, Hmax=X/2+2;
  uint32_t *psi=malloc(Hmax*4); psi[0]=0;
  for(i64 h=1;h<Hmax;h++){ i64 q=h>>1; psi[h]=(h&1)? psi[q]+(uint32_t)(q>>1)+1 : psi[q]+(uint32_t)q; }
  uint8_t *r=calloc(X+1,1);
  for(int s=0;s<K;s++) for(i64 h1=0;h1<Hmax;h1++){
    i64 y1=(2*h1+1)<<s; if(y1>X) break; i64 g1=2*(i64)psi[h1]+1; i64 lim=X/y1; i64 h2max=(lim-1)/2; if(h2max<0) break;
    for(i64 h2=0;h2<=h2max;h2++){ i64 N=g1*(2*(i64)psi[h2]+1)+y1*(2*h2+1)-(s==0); if(N<=X && r[N]<255) r[N]++; } }
  free(psi);
  // sieve: smallest prime factor flag for odd numbers up to X+1
  uint8_t *comp=calloc(X/2+2,1);           // comp[n>>1]=1 if odd n composite
  for(i64 p=3;p*p<=X+1;p+=2) if(!comp[p>>1]) for(i64 q=p*p;q<=X+1;q+=2*p) comp[q>>1]=1;
  int SP[30],ns=0; for(int p=3;p<=97;p+=2){ int ok=1; for(int d=3;d*d<=p;d+=2) if(p%d==0) ok=0; if(ok) SP[ns++]=p; }
  // strata: w = #small primes dividing N+1 (cap 4), c = N mod 16
  u64 pn[5][16]={{0}}, pm[5][16]={{0}}, cn[5][16]={{0}}, cm[5][16]={{0}};
  for(i64 N=L|1;N<X;N+=2){
    int rough=1; for(int i=0;i<ns;i++) if(N%SP[i]==0){ rough=0; break; } if(!rough) continue;
    int w=0; for(int i=0;i<ns;i++) if((N+1)%SP[i]==0) w++; if(w>4) w=4; int c=N&15;
    int isp=!comp[N>>1];
    if(isp){ pn[w][c]++; if(!r[N]) pm[w][c]++; } else { cn[w][c]++; if(!r[N]) cm[w][c]++; } }
  printf("F1, odd N in [2^%d,2^%d) with no prime factor <= 97: primes vs rough composites\n",lo,K);
  u64 TP=0,TPm=0,TC=0,TCm=0; double expP=0;
  for(int w=0;w<5;w++){ u64 a=0,b=0,c_=0,d=0; double e=0;
    for(int c=0;c<16;c++){ a+=pn[w][c]; b+=pm[w][c]; c_+=cn[w][c]; d+=cm[w][c]; if(cn[w][c]) e+= (double)pn[w][c]*cm[w][c]/cn[w][c]; }
    printf("  #small p | N+1 = %d%s: primes n=%9llu miss %.5f | rough composites n=%9llu miss %.5f | prime misses observed %llu vs composite-rate expectation (matched mod 16) %.1f\n",
      w,w==4?"+":"",a,a?(double)b/a:0,c_,c_?(double)d/c_:0,b,e);
    TP+=a;TPm+=b;TC+=c_;TCm+=d;expP+=e; }
  printf("  total: primes %llu misses %llu (rate %.5f); rough composites %llu misses %llu (rate %.5f); matched ratio primes/composites = %.3f\n",
    TP,TPm,(double)TPm/TP,TC,TCm,(double)TCm/TC,TPm/expP);
  return 0; }
