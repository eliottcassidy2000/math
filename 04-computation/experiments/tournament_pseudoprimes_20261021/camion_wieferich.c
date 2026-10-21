/* camion_wieferich.c -- primes p in [P0,P1] with s(p) = 2^((p-1)/2) (mod p^2), s = strong labelled tournaments.
 * O(p^2) recurrence mod p^2 per prime.  mac-mini-2026-10-21-fanopsp */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 pw(u64 b,u64 e,u64 m){u64 r=1%m;b%=m;while(e){if(e&1)r=(u128)r*b%m;b=(u128)b*b%m;e>>=1;}return r;}
int main(int argc,char**argv){ int P0=atoi(argv[1]),P1=atoi(argv[2]);
  u64 *T=malloc((P1+2)*8),*s=malloc((P1+2)*8),*row=malloc((P1+2)*8);
  for(int p=P0;p<=P1;p++){ if(p<3||p%2==0) continue; int ok=1; for(int d=3;d*d<=p;d+=2) if(p%d==0){ok=0;break;} if(!ok) continue;
    u64 m=(u64)p*p; for(int j=0;j<=p;j++) T[j]=pw(2,(u64)j*(j-1)/2,m); row[0]=1;
    for(int k=1;k<=p;k++){ for(int i=k-1;i>=1;i--){ row[i]+=row[i-1]; if(row[i]>=m) row[i]-=m; } row[k]=1;
      u128 acc=0; for(int i=1;i<k;i++){ acc+=(u128)((u128)row[i]*s[i]%m)*T[k-i]; if((i&15)==0) acc%=m; } acc%=m;
      s[k]=(T[k]+m-(u64)acc)%m; }
    u64 c=pw(2,(p-1)/2,m); if(s[p]==c) { printf("CAMION-WIEFERICH %d\n",p); fflush(stdout); }
    if(p%1000<2) { fprintf(stderr,"done p=%d\n",p); } }
  printf("searched %d..%d\n",P0,P1); return 0; }
