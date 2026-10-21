/* camion_list.c -- (A) squarefree Camion pseudoprimes n <= N via the Korselt-type criterion (zero classes Z_p);
 *                  (B) Camion-Wieferich primes: D(p) = s(p)-2^((p-1)/2) = 0 mod p^2, p <= PW.
 * mac-mini-2026-10-21-fanopsp */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 pw(u64 b,u64 e,u64 m){u64 r=1%m;b%=m;while(e){if(e&1)r=(u128)r*b%m;b=(u128)b*b%m;e>>=1;}return r;}
static void smod(int M,u64 m,u64*s){ u64*T=malloc((M+1)*8),*row=malloc((M+2)*8);
  for(int j=0;j<=M;j++) T[j]=pw(2,(u64)j*(j-1)/2,m); row[0]=1%m;
  for(int k=1;k<=M;k++){ for(int i=k-1;i>=1;i--) row[i]=(row[i]+row[i-1])%m; row[k]=1%m; /* row = C(k, .) in place */
    u64 acc=0; for(int i=1;i<k;i++) acc=(acc+(u128)row[i]*s[i]%m*T[k-i])%m; s[k]=(T[k]+m-acc)%m; }
  free(T);free(row); }
int main(int argc,char**argv){ int N=atoi(argv[1]), PW=atoi(argv[2]); int PM=N/2+2;
  char *isc=calloc(PM+1,1); for(int i=2;i<=PM;i++) if(!isc[i]) for(long j=(long)i*i;j<=PM;j+=i) isc[j]=1;
  char **Z=calloc(PM+1,sizeof(char*));
  for(int p=3;p<=PM;p+=2){ if(isc[p]) continue; u64*s=malloc((p+1)*8); smod(p,p,s); Z[p]=calloc(p,1);
    for(int r=1;r<=p-1;r++){ u64 c=(r&1)?pw(2,(r-1)/2,p):0; Z[p][r]=((s[r]+p-c)%p==0); } free(s); }
  /* (A) */
  FILE*fo=fopen("camion_squarefree_list.txt","w"); int cnt=0, carm=0, carmpass=0;
  for(int n=4;n<=N;n++){ int x=n, sqf=1, np=0, ps[32], prime=1; if(x%4==0) sqf=0; if(x%2==0) x/=2;
    for(int d=3;(long)d*d<=x && sqf;d+=2) if(x%d==0){ ps[np++]=d; x/=d; if(x%d==0) sqf=0; }
    if(!sqf) continue; if(x>1) ps[np++]=x; if(n%2==0) prime=0; else prime=(np==1 && ps[0]==n);
    if(prime || n<4) continue; if(np==0) continue; /* n=2 */
    int ok=1; for(int i=0;i<np;i++){ int p=ps[i]; int r=((n-1)%(p-1))+1; if(!Z[p][r]) {ok=0;break;} }
    /* Carmichael? odd squarefree, >=3 primes, p-1 | n-1 */
    int isCarm = (n%2==1 && np>=2); for(int i=0;i<np && isCarm;i++) if((n-1)%(ps[i]-1)) isCarm=0;
    if(isCarm){carm++; carmpass+=ok;}
    if(ok){ cnt++; fprintf(fo,"%d%s\n",n,(n%2)?"":" even"); } }
  fclose(fo); printf("(A) squarefree Camion pseudoprimes <= %d: %d ; Carmichael numbers %d of which pass %d\n",N,cnt,carm,carmpass);
  /* (B) */
  printf("(B) Camion-Wieferich primes p<=%d (D(p)=0 mod p^2):",PW); fflush(stdout);
  for(int p=3;p<=PW;p+=2){ int q=1; for(int d=3;d*d<=p;d+=2) if(p%d==0){q=0;break;} if(!q) continue;
    u64 m=(u64)p*p; u64*s=malloc((p+1)*8); smod(p,m,s); if(s[p]==pw(2,(p-1)/2,m)) {printf(" %d",p); fflush(stdout);} free(s); }
  printf("\n"); return 0; }
