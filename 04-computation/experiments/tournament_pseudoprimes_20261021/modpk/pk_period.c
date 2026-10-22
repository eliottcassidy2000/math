/* pk_period.c -- test s(m + phi(p^k)) = (2/p) s(m) (mod p^k), s = strong labelled tournaments (A054946),
 * and the same for D = s - c (c = circulant count), for all 1 <= m <= M - L, L = phi(p^k), M = 2L + 4k + 20.
 * Prints the full list of failing m (expected: only m <= 2k-2, by THM-4748).  mac-mini-2026-10-21-fanopsp
 * Usage: ./pk_period p k */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 pw(u64 b,u64 e,u64 m){u64 r=1%m;b%=m;while(e){if(e&1)r=(u128)r*b%m;b=(u128)b*b%m;e>>=1;}return r;}
int main(int argc,char**argv){ u64 p=atoll(argv[1]); int k=atoi(argv[2]);
  u64 q=1; for(int i=0;i<k;i++) q*=p; u64 L=q/p*(p-1); long M=2*L+4*k+20;
  u64 *T=malloc((M+1)*8),*s=malloc((M+1)*8),*row=malloc((M+2)*8);
  for(long j=0;j<=M;j++){ u64 e=((u64)j*(j-1)/2); T[j]=pw(2,e,q); }
  row[0]=1%q;
  for(long m=1;m<=M;m++){ for(long i=m-1;i>=1;i--){ row[i]+=row[i-1]; if(row[i]>=q) row[i]-=q; } row[m]=1%q;
    u128 acc=0; for(long i=1;i<m;i++){ acc+=(u128)((u128)row[i]*s[i]%q)*T[m-i]; if((i&7)==0) acc%=q; } acc%=q;
    s[m]=(T[m]+q-(u64)acc)%q; }
  int eps=(p%8==1||p%8==7)?1:-1;
  printf("p=%llu k=%d q=%llu L=%llu M=%ld eps=%+d\n",p,k,q,L,M,eps);
  int nf=0,nfD=0; printf("  s-law fails at m:"); 
  for(long m=1;m+(long)L<=M;m++){ u64 rhs=(eps==1)?s[m]:(q-s[m])%q; if(s[m+L]!=rhs){ printf(" %ld",m); nf++; if(nf>40){printf(" ...");break;} } }
  printf("  (count %d)\n  D-law fails at m:",nf);
  for(long m=1;m+(long)L<=M;m++){ u64 c=(m&1)?pw(2,(m-1)/2,q):0, c2=((m+L)&1)?pw(2,(m+L-1)/2,q):0;
    u64 D1=(s[m]+q-c)%q, D2=(s[m+L]+q-c2)%q; u64 rhs=(eps==1)?D1:(q-D1)%q; if(D2!=rhs){ printf(" %ld",m); nfD++; if(nfD>40){printf(" ...");break;} } }
  printf("  (count %d)\n  threshold 2k-1=%d : %s\n",nfD,2*k-1, (nf==0||1)?"":"");
  /* verify that no failure occurs at m >= 2k-1 */
  int bad=0; for(long m=2*k-1;m+(long)L<=M;m++){ u64 rhs=(eps==1)?s[m]:(q-s[m])%q; if(s[m+L]!=rhs) bad++; }
  printf("  failures at m >= 2k-1: %d  %s\n",bad, bad?"*** CONTRADICTS THM ***":"OK");
  return bad?1:0; }
