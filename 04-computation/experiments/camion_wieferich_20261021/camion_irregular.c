/* camion_irregular.c -- Camion zero classes for every prime p in [P0,P1]: r in [4, p-1] with p | D(r) = s(r) - c(r),
 * using s(r) = -r! * u_r (mod p), u_r = [x^r] 1/T(x) (r < p).  Prints per-prime counts of odd / even sporadic classes
 * and a summary (fraction of primes with no odd sporadic class, the analogue of Siegel's regular-prime density).
 * mac-mini-2026-10-21-fanopsp.  Usage: ./camion_irregular P0 P1 [list] */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 pw(u64 b,u64 e,u64 m){u64 r=1%m;b%=m;while(e){if(e&1)r=(u128)r*b%m;b=(u128)b*b%m;e>>=1;}return r;}
int main(int argc,char**argv){ u64 P0=atoll(argv[1]),P1=atoll(argv[2]); int list=(argc>3);
  uint32_t *t=malloc((P1+2)*4),*rev=malloc((P1+2)*4); u64 *fact=malloc((P1+2)*8),*invfact=malloc((P1+2)*8);
  char *comp=calloc(P1+2,1); for(u64 i=2;i*i<=P1;i++) if(!comp[i]) for(u64 j=i*i;j<=P1;j+=i) comp[j]=1;
  long np=0, oddreg=0, evenreg=0, bothreg=0; long hodd[16]={0}, heven[16]={0};
  for(u64 p=(P0<7?7:P0);p<=P1;p++){ if(comp[p]||p%2==0) continue;
    fact[0]=1; for(u64 k=1;k<p;k++) fact[k]=fact[k-1]*k%p; invfact[p-1]=pw(fact[p-1],p-2,p); for(u64 k=p-1;k>=1;k--) invfact[k-1]=invfact[k]*k%p;
    u64 two=1; for(u64 k=0;k<p;k++){ t[k]=(uint32_t)(two*invfact[k]%p); two=two*pw(2,k,p)%p; }
    u64 M=p; rev[M]=1;
    int nodd=0,neven=0; char buf[4096]; int bl=0; buf[0]=0;
    for(u64 nn=1;nn<p;nn++){ u64 acc=0; const uint32_t *tt=t+1,*rr=rev+(M-nn+1); for(u64 k=0;k<nn;k++) acc+=(u64)tt[k]*rr[k];
      u64 un=(p-acc%p)%p; rev[M-nn]=(uint32_t)un;
      if(nn>=4){ u64 s=(p - fact[nn]*un%p)%p; u64 c=(nn&1)?pw(2,(nn-1)/2,p):0; if(s==c){ if(nn&1) nodd++; else neven++; if(list && bl<4000) bl+=sprintf(buf+bl," %llu",(unsigned long long)nn); } } }
    np++; if(!nodd) oddreg++; if(!neven) evenreg++; if(!nodd&&!neven) bothreg++; hodd[nodd<15?nodd:15]++; heven[neven<15?neven:15]++;
    if(list) printf("p=%llu odd=%d even=%d classes:%s\n",(unsigned long long)p,nodd,neven,buf); }
  printf("primes %llu..%llu: %ld; no odd sporadic class: %ld (%.4f; Siegel-type e^-1/2 = 0.6065); no even: %ld (%.4f); neither: %ld (%.4f; e^-1 = 0.3679)\n",
    (unsigned long long)P0,(unsigned long long)P1,np,oddreg,(double)oddreg/np,evenreg,(double)evenreg/np,bothreg,(double)bothreg/np);
  printf("odd-class count histogram:"); for(int i=0;i<8;i++) printf(" %d:%ld",i,hodd[i]); printf("\neven-class count histogram:"); for(int i=0;i<8;i++) printf(" %d:%ld",i,heven[i]); printf("\n");
  return 0; }
