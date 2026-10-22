/* korselt_full.c -- THM-4748 checks.
 * (A) Full Korselt-type criterion for EVERY composite n in [4, N]: for each odd prime power p^a || n,
 *     s(n) = c(n) (mod p^a)  iff  p^a | D(r), r = the integer in [2a-1, 2a-2+phi(p^a)] with r = n (mod phi(p^a));
 *     the 2-part is automatic.  Compared with direct s(n) mod n.
 * (B) Tower: D(p^a) = (2/p) D(p^(a-1)) (mod p^a) for odd prime powers p^a <= TMAX.
 * mac-mini-2026-10-21-fanopsp.  Usage: ./korselt_full N TMAX */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 pw(u64 b,u64 e,u64 m){u64 r=1%m;b%=m;while(e){if(e&1)r=(u128)r*b%m;b=(u128)b*b%m;e>>=1;}return r;}
static u64 *row,*T;
static void smod(long M,u64 m,u64*s){ for(long j=0;j<=M;j++) T[j]=pw(2,(u64)j*(j-1)/2,m); row[0]=1%m;
  for(long k=1;k<=M;k++){ for(long i=k-1;i>=1;i--){ row[i]+=row[i-1]; if(row[i]>=m) row[i]-=m; } row[k]=1%m;
    u128 acc=0; for(long i=1;i<k;i++){ acc+=(u128)((u128)row[i]*s[i]%m)*T[k-i]; if((i&7)==0) acc%=m; } acc%=m; s[k]=(T[k]+m-(u64)acc)%m; } }
static u64 cc(long m,u64 q){ return (m&1)? pw(2,(m-1)/2,q) : (m==0?1%q:0); }
int main(int argc,char**argv){ long N=atol(argv[1]), TMAX=atol(argv[2]); long BIG=(N>TMAX?N:TMAX)*3+64;
  row=malloc((BIG+2)*8); T=malloc((BIG+2)*8); u64 *s=malloc((BIG+2)*8);
  /* (A) */
  long agree=0,dis=0,pass=0,nonsqf=0,nonsqfpass=0;
  for(long n=4;n<=N;n++){ long x=n; int prime=1; for(long d=2;d*d<=x;d++) if(x%d==0){prime=0;break;} if(prime) continue;
    int pred=1, sqf=1; long y=n; while(y%2==0) y/=2;
    for(long p=3;p<=y;p+=2){ if(y%p) continue; int a=0; u64 q=1; while(y%p==0){ y/=p; a++; q*=p; } if(a>1) sqf=0;
      u64 L=q/p*(p-1); long lo=2*a-1; long r=lo+((n-lo)%(long)L+(long)L)%(long)L; /* r in [lo, lo+L-1], r = n mod L */
      smod(r,q,s); u64 D=(s[r]+q-cc(r,q))%q; if(D!=0) pred=0; }
    smod(n,n,s); int direct=(s[n]==cc(n,n));
    if(pred==direct) agree++; else { dis++; printf("DISAGREE n=%ld pred=%d direct=%d\n",n,pred,direct); }
    pass+=direct; if(!sqf){ nonsqf++; nonsqfpass+=direct; } }
  printf("(A) composites 4..%ld: agree %ld, disagree %ld; passing %ld; with a squared odd prime factor: %ld (passing %ld)\n",N,agree,dis,pass,nonsqf,nonsqfpass);
  /* (B) tower */
  long tb=0,tt=0;
  for(u64 p=3;p<=(u64)TMAX;p+=2){ int pr=1; for(u64 d=3;d*d<=p;d+=2) if(p%d==0){pr=0;break;} if(!pr) continue;
    int eps=(p%8==1||p%8==7)?1:-1; u64 q=p, prev=1; int a=1;
    while(q<=(u64)TMAX){ smod((long)q,q,s); u64 Dq=(s[q]+q-cc((long)q,q))%q; u64 Dprev=(s[prev]+q-cc((long)prev,q))%q;
      u64 rhs=(eps==1)?Dprev:(q-Dprev)%q; tt++; if(Dq!=rhs){ tb++; printf("TOWER FAIL p=%llu a=%d\n",p,a); }
      if(a>=2 && Dq==0) printf("  p^%d = %llu passes (p^a | D(p^(a-1)))\n",a,q);
      prev=q; q*=p; a++; } }
  printf("(B) tower checks: %ld, failures %ld\n",tt,tb); return 0; }
