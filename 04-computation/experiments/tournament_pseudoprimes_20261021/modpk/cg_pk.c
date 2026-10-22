/* cg_pk.c -- connected labelled graphs cg(n) (A001187) mod p^k: test cg(n + phi(p^k)) = (2/p) cg(n) (mod p^k)
 * for 1 <= n <= M - L; THM-4748's argument (shift-closed algebra) predicts failures only at n <= 2k-1.
 * Recurrence: 2^C(n,2) = sum_{j=1}^{n} C(n-1,j-1) cg(j) 2^C(n-j,2).  mac-mini-2026-10-21-fanopsp.  Usage: ./cg_pk p k */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 pw(u64 b,u64 e,u64 m){u64 r=1%m;b%=m;while(e){if(e&1)r=(u128)r*b%m;b=(u128)b*b%m;e>>=1;}return r;}
int main(int argc,char**argv){ u64 p=atoll(argv[1]); int k=atoi(argv[2]); u64 q=1; for(int i=0;i<k;i++) q*=p;
  u64 L=q/p*(p-1); long M=2*L+4*k+20; int eps=(p%8==1||p%8==7)?1:-1;
  u64 *G=malloc((M+1)*8),*c=malloc((M+1)*8),*row=malloc((M+2)*8);
  for(long j=0;j<=M;j++) G[j]=pw(2,(u64)j*(j-1)/2,q);
  row[0]=1; /* row = C(n-1, .) */
  for(long n=1;n<=M;n++){ if(n>=2){ for(long i=n-2;i>=1;i--){ row[i]+=row[i-1]; if(row[i]>=q) row[i]-=q; } row[n-1]=1; }
    u128 acc=0; for(long j=1;j<n;j++){ acc+=(u128)((u128)row[j-1]*c[j]%q)*G[n-j]; if((j&7)==0) acc%=q; } acc%=q; c[n]=(G[n]+q-(u64)acc)%q; }
  printf("p=%llu k=%d L=%llu: cg-law fails at n:",p,k,L); int bad=0;
  for(long n=1;n+(long)L<=M;n++){ u64 rhs=(eps==1)?c[n]:(q-c[n])%q; if(c[n+L]!=rhs){ printf(" %ld",n); if(n>=2*k) bad++; } }
  printf("   failures at n >= 2k: %d %s\n",bad,bad?"*** CONTRADICTS ***":"OK"); return 0; }
