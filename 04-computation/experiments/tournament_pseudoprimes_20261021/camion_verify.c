/* camion_verify.c -- (1) periodicity s(m+p-1) = (2/p) s(m) mod p; (2) zero classes Z_p of D(r)=s(r)-2^((r-1)/2);
 * (3) Korselt-type prediction vs direct test for squarefree odd composite n <= NMAX;
 * (4) Camion-Wieferich primes: p^2 passes iff F(p) = 2^((p^2-1)/2) mod p^2, F = strong tournaments on Z_{p^2}
 *     invariant under x -> x+p (recurrence over cosets, t_k = 2^(k(kp-1)/2)).   mac-mini-2026-10-21-fanopsp */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 pw(u64 b,u64 e,u64 m){u64 r=1%m;b%=m;while(e){if(e&1)r=(u128)r*b%m;b=(u128)b*b%m;e>>=1;}return r;}
static int isprime(u64 n){if(n<2)return 0;for(u64 d=2;d*d<=n;d++)if(n%d==0)return 0;return 1;}
/* s(1..M) mod m */
static void smod(int M,u64 m,u64*s){ u64*T=malloc((M+1)*8),*row=malloc((M+2)*8),*nr=malloc((M+2)*8);
  for(int j=0;j<=M;j++) T[j]=pw(2,(u64)j*(j-1)/2,m); row[0]=1%m;
  for(int k=1;k<=M;k++){ nr[0]=1%m; for(int i=1;i<k;i++) nr[i]=(row[i-1]+row[i])%m; nr[k]=1%m; memcpy(row,nr,(k+1)*8);
    u64 acc=0; for(int i=1;i<k;i++) acc=(acc+(u128)row[i]*s[i]%m*T[k-i])%m; s[k]=(T[k]+m-acc)%m; }
  free(T);free(row);free(nr); }
int main(int argc,char**argv){ int P=atoi(argv[1]), NMAX=atoi(argv[2]), PW=atoi(argv[3]);
  static char Z[4000][4000]; /* Z[p][r] */
  int perfail=0;
  for(int p=3;p<=P;p+=2){ if(!isprime(p)) continue; int M=3*p+5; u64*s=malloc((M+1)*8); smod(M,p,s);
    int eps=(p%8==1||p%8==7)?1:-1;
    for(int m=1;m+p-1<=M;m++){ u64 lhs=s[m+p-1], rhs=(eps==1)?s[m]:(p-s[m])%p; if(lhs!=rhs){perfail++; printf("PERIOD FAIL p=%d m=%d\n",p,m);} }
    int nz=0; printf("p=%d eps=%+d zero classes r mod %d:",p,eps,p-1);
    for(int r=1;r<=p-2;r+=2){ u64 D=(s[r]+p-pw(2,(r-1)/2,p))%p; Z[p][r]=(D==0); if(D==0){printf(" %d",r);nz++;} }
    printf("  (%d of %d)\n",nz,(p-1)/2); free(s); }
  printf("periodicity failures: %d\n",perfail);
  /* (3) Korselt prediction vs direct */
  int agree=0,disagree=0,npass=0;
  for(int n=9;n<=NMAX;n+=2){ if(isprime(n)) continue; int x=n,sqf=1; int ps[20],np=0;
    for(int d=3;d*d<=x;d+=2) if(x%d==0){ ps[np++]=d; x/=d; if(x%d==0){sqf=0;break;} }
    if(!sqf) continue; if(x>1) ps[np++]=x;
    int pred=1; for(int i=0;i<np;i++){ int p=ps[i]; if(p>P){pred=-1;break;} int r=n%(p-1); if(!Z[p][r]) pred=0; }
    if(pred<0) continue;
    u64*s=malloc((n+1)*8); smod(n,n,s); int direct=(s[n]==pw(2,(n-1)/2,n)); free(s);
    if(direct==pred) agree++; else {disagree++; printf("DISAGREE n=%d pred=%d direct=%d\n",n,pred,direct);}
    npass+=direct; }
  printf("Korselt check (squarefree odd composite n<=%d with all p<=%d): agree %d disagree %d, passing %d\n",NMAX,P,agree,disagree,npass);
  /* (4) Camion-Wieferich */
  printf("Camion-Wieferich primes p<=%d (p^2 passes):",PW);
  for(int p=3;p<=PW;p+=2){ if(!isprime(p)) continue; u64 m=(u64)p*p; int K=p;
    u64*t=malloc((K+1)*8),*F=malloc((K+1)*8),*row=malloc((K+2)*8),*nr=malloc((K+2)*8);
    for(int k=0;k<=K;k++) t[k]=pw(2,(u64)k*((u64)k*p-1)/2,m); row[0]=1;
    for(int k=1;k<=K;k++){ nr[0]=1; for(int i=1;i<k;i++) nr[i]=(row[i-1]+row[i])%m; nr[k]=1; memcpy(row,nr,(k+1)*8);
      u64 acc=0; for(int i=1;i<k;i++) acc=(acc+(u128)row[i]*F[i]%m*t[k-i])%m; F[k]=(t[k]+m-acc)%m; }
    if(F[K]==pw(2,(m-1)/2,m)) printf(" %d",p); fflush(stdout); free(t);free(F);free(row);free(nr); }
  printf("\n"); return 0; }
