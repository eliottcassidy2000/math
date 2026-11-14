// Wave 23 (P48/P49): square-class misses of F1 against an empirical local model with strata
// (w0 = #primes 11<=p<=97 dividing N, w1 = #primes 11<=p<=97 dividing N+1, N mod 16, N mod 105).
// Usage: f1strata105 K lo   (in-sample: 27 22; disjoint out-of-sample window: 29 27)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef long long i64; typedef unsigned long long u64;
static int SP[30],ns=0;
static int stratum(i64 N){ int w0=0,w1=0; for(int i=0;i<ns;i++){ if(N%SP[i]==0) w0++; if((N+1)%SP[i]==0) w1++; }
  if(w0>4) w0=4; if(w1>5) w1=5; return ((w0*6+w1)*16+(int)(N&15))*105 + (int)(N%105); }
int main(int argc,char**argv){
  int K=atoi(argv[1]), lo=atoi(argv[2]); i64 X=1LL<<K, L=1LL<<lo, Hmax=X/2+2;
  for(int p=11;p<=97;p+=2){ int ok=1; for(int d=3;d*d<=p;d+=2) if(p%d==0) ok=0; if(ok) SP[ns++]=p; }
  uint32_t *psi=malloc(Hmax*4); psi[0]=0;
  for(i64 h=1;h<Hmax;h++){ i64 q=h>>1; psi[h]=(h&1)? psi[q]+(uint32_t)(q>>1)+1 : psi[q]+(uint32_t)q; }
  uint8_t *r=calloc(X+1,1);
  for(int s=0;s<K;s++) for(i64 h1=0;h1<Hmax;h1++){
    i64 y1=(2*h1+1)<<s; if(y1>X) break; i64 g1=2*(i64)psi[h1]+1; i64 lim=X/y1; i64 h2max=(lim-1)/2; if(h2max<0) break;
    for(i64 h2=0;h2<=h2max;h2++){ i64 N=g1*(2*(i64)psi[h2]+1)+y1*(2*h2+1)-(s==0); if(N<=X && r[N]<255) r[N]++; } }
  free(psi);
  const int NS=5*6*16*105; u64 *sn=calloc(NS,8), *sm=calloc(NS,8);
  for(i64 N=L|1;N<X;N+=2){ int s=stratum(N); sn[s]++; if(!r[N]) sm[s]++; }
  printf("F1 odd N in [2^%d,2^%d): expected misses from %d empirical strata (w0,w1 over primes 11..97, N mod 16, N mod 105)\n",lo,K,NS);
  int ts[]={1,3,5,7,11,13,15,19,21,35}; 
  for(int ti=0;ti<10;ti++){ i64 t=ts[ti]; u64 n=0,obs=0; double e=0;
    for(i64 m=1;;m+=2){ i64 N=t*m*m; if(N>=X) break; if(N<L) continue; i64 a=m,b=t; while(b){i64 c=a%b;a=b;b=c;} if(t>1&&a!=1) continue;
      int s=stratum(N); n++; e+=(double)sm[s]/sn[s]; if(!r[N]) obs++; }
    double z=(obs-e)/__builtin_sqrt(e>0?e:1);
    printf("  t=%3lld: n=%6llu observed %4llu expected %7.1f ratio %.2f  z=%+.2f\n",t,n,obs,e,e>0?obs/e:-1.0,z); }
  return 0; }
