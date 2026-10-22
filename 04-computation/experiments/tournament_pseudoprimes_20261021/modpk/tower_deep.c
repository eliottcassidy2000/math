/* tower_deep.c -- D(p^a) = (2/p) D(p^(a-1)) (mod p^a) for given p and a = 2..A, plus p^a | D(p^(a-1)) report.
 * mac-mini-2026-10-21-fanopsp.  Usage: ./tower_deep p A */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 pw(u64 b,u64 e,u64 m){u64 r=1%m;b%=m;while(e){if(e&1)r=(u128)r*b%m;b=(u128)b*b%m;e>>=1;}return r;}
int main(int argc,char**argv){ u64 p=atoll(argv[1]); int A=atoi(argv[2]); int eps=(p%8==1||p%8==7)?1:-1;
  u64 q=1; for(int i=0;i<A;i++) q*=p; long M=(long)q; /* need s up to p^A mod p^a for each a: compute mod p^A once */
  u64 *T=malloc((M+1)*8),*s=malloc((M+1)*8),*row=malloc((M+2)*8);
  for(long j=0;j<=M;j++) T[j]=pw(2,(u64)j*(j-1)/2,q); row[0]=1;
  for(long k=1;k<=M;k++){ for(long i=k-1;i>=1;i--){ row[i]+=row[i-1]; if(row[i]>=q) row[i]-=q; } row[k]=1;
    u128 acc=0; for(long i=1;i<k;i++){ acc+=(u128)((u128)row[i]*s[i]%q)*T[k-i]; if((i&7)==0) acc%=q; } acc%=q; s[k]=(T[k]+q-(u64)acc)%q; }
  u64 pa=1; long prev=1;
  for(int a=1;a<=A;a++){ pa*=p; long n=(long)pa; u64 cn=pw(2,(n-1)/2,pa), cp=(prev&1)?pw(2,(prev-1)/2,pa):0;
    u64 Dn=(s[n]%pa+pa-cn)%pa, Dp=(s[prev]%pa+pa-cp)%pa; u64 rhs=(eps==1)?Dp:(pa-Dp)%pa;
    /* p-adic valuation of D(p^(a-1)) and of D(p^a), computed mod p^A */
    u64 Dpa=(s[prev]+q-((prev&1)?pw(2,(prev-1)/2,q):0))%q, Dna=(s[n]+q-pw(2,(n-1)/2,q))%q;
    int vp=0; u64 t=Dpa; if(t==0) vp=A; else while(t%p==0){t/=p;vp++;}
    int vn=0; t=Dna; if(t==0) vn=A; else while(t%p==0){t/=p;vn++;}
    printf("p=%llu a=%d: tower %s;  v_p(D(p^%d)) = %s%d,  v_p(D(p^%d)) = %s%d;  p^%d is %sa Camion pseudoprime\n",p,a,(Dn==rhs)?"OK":"FAIL",a-1,(vp==A)?">=":"",vp,a,(vn==A)?">=":"",vn,a,(Dn==0)?"":"not ");
    prev=n; }
  return 0; }
