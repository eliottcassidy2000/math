/* cw_fast.c -- Camion quotient q_C(p) = D(p)/p mod p for primes p in [P0,P1], using
 *   s(p) = eps + p*u_p (mod p^2),  u_p = [x^p] 1/T_{<p}(x),  T_{<p} = sum_{k<p} 2^C(k,2) x^k/k!   (arithmetic mod p only),
 *   c(p) = 2^((p-1)/2) = eps + p*e' (mod p^2),  so  q_C(p) = u_p - e' (mod p);  Camion-Wieferich iff q_C(p) = 0.
 * Inner loop: u_n = -sum_{k=1}^{n} t_k u_{n-k}, 64-bit accumulation without reduction (valid for p < 2.6e6).
 * mac-mini-2026-10-21-fanopsp.  Usage: ./cw_fast P0 P1 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 pw(u64 b,u64 e,u64 m){u64 r=1%m;b%=m;while(e){if(e&1)r=(u128)r*b%m;b=(u128)b*b%m;e>>=1;}return r;}
int main(int argc,char**argv){ u64 P0=atoll(argv[1]),P1=atoll(argv[2]);
  uint32_t *t=malloc((P1+2)*4),*rev=malloc((P1+2)*4); u64 *fact=malloc((P1+2)*8);
  char *comp=calloc(P1+2,1); for(u64 i=2;i*i<=P1;i++) if(!comp[i]) for(u64 j=i*i;j<=P1;j+=i) comp[j]=1;
  for(u64 p=(P0<3?3:P0);p<=P1;p++){ if(comp[p]||p%2==0) continue;
    /* t_k = 2^C(k,2) / k!  mod p, k = 0..p-1 */
    fact[0]=1; for(u64 k=1;k<p;k++) fact[k]=fact[k-1]*k%p;
    u64 invf=pw(fact[p-1],p-2,p); /* inverse of (p-1)! */
    u64 two=1; /* 2^C(k,2) */
    u64 *invfact=malloc(p*8); invfact[p-1]=invf; for(u64 k=p-1;k>=1;k--) invfact[k-1]=invfact[k]*k%p;
    for(u64 k=0;k<p;k++){ t[k]=(uint32_t)(two*invfact[k]%p); two=two*pw(2,k,p)%p; }
    free(invfact);
    /* u_0 = 1; rev[p - j] = u_j */
    u64 M=p; rev[M]=1;
    for(u64 nn=1;nn<=p;nn++){ u64 acc=0; u64 kmax=(nn<p)?nn:p-1; /* t_p = 0 for T_{<p} */
      const uint32_t *tt=t+1, *rr=rev+(M-nn+1);
      for(u64 k=0;k<kmax;k++) acc+=(u64)tt[k]*rr[k];
      u64 un=(p-acc%p)%p; rev[M-nn]=(uint32_t)un; }
    u64 up=rev[0];
    int eps=(p%8==1||p%8==7)?1:-1; u64 p2=p*p; u64 h=pw(2,(p-1)/2,p2); u64 hm=(eps==1)?(h+p2-1)%p2:(h+1)%p2; u64 e1=(hm/p)%p;
    u64 qc=(up+p-e1)%p;
    if(qc==0){ printf("CAMION-WIEFERICH p=%llu\n",(unsigned long long)p); fflush(stdout); }
    if(p<=40) { printf("  p=%llu q_C=%llu\n",(unsigned long long)p,(unsigned long long)qc); }
  }
  printf("searched [%llu,%llu]\n",(unsigned long long)P0,(unsigned long long)P1); return 0; }
