/* camion_psp.c -- "Camion pseudoprimes": composite n for which the number s(n) of strong
 * (equivalently, by Camion's theorem, Hamiltonian) labeled tournaments satisfies the
 * prime-order congruence s(n) = 2^((n-1)/2) (mod n), together with the tournament-count
 * test T(n) = 2^C(n,2) = 2^((n-1)/2) (mod n).  mac-mini-2026-10-21-fanopsp.
 * Recurrence: T(n) = sum_{k=1}^{n} C(n,k) s(k) T(n-k), T(0)=1.  Everything mod n.
 * Usage: ./camion_psp NMAX
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 pw(u64 b,u64 e,u64 m){ u64 r=1%m; b%=m; while(e){ if(e&1) r=(u128)r*b%m; b=(u128)b*b%m; e>>=1;} return r; }
static int isprime(u64 n){ if(n<2) return 0; for(u64 d=2;d*d<=n;d++) if(n%d==0) return 0; return 1; }
int main(int argc,char**argv){ int NMAX=atoi(argv[1]);
  u64 *s=malloc((NMAX+1)*sizeof(u64)), *T=malloc((NMAX+1)*sizeof(u64)), *row=malloc((NMAX+2)*sizeof(u64)), *nr=malloc((NMAX+2)*sizeof(u64));
  for(int n=3;n<=NMAX;n++){ u64 m=n;
    /* T(j) mod m */
    for(int j=0;j<=n;j++) T[j]=pw(2,(u64)j*(j-1)/2,m);
    /* s(k) mod m for k=1..n, using binomials C(k,i) mod m built row by row */
    row[0]=1%m; /* row = C(k-1, .) */
    for(int k=1;k<=n;k++){ /* build row for C(k, i) */
      nr[0]=1%m; for(int i=1;i<k;i++) nr[i]=(row[i-1]+row[i])%m; nr[k]=1%m;
      for(int i=0;i<=k;i++) row[i]=nr[i];
      u64 acc=0; for(int i=1;i<k;i++) acc=(acc+(u128)row[i]*s[i]%m*T[k-i])%m;
      s[k]=(T[k]+m-acc)%m; }
    if(n%2==1){ u64 circ=pw(2,(n-1)/2,m); int pr=isprime(n);
      int camion=(s[n]==circ), tfp=(T[n]==circ);
      /* Jacobi symbol (2/n) */
      int r8=n%8; int jac=(r8==1||r8==7)?1:-1; u64 jm=(jac==1)?1%m:m-1;
      int euler=(circ==jm);
      if(!pr && (camion||tfp)) printf("n=%d composite: camion=%d tfp=%d eulerjacobi2=%d s mod n=%llu circ=%llu\n",n,camion,tfp,euler,s[n],circ);
      if(pr && !camion) printf("PRIME FAILURE n=%d\n",n); }
    else { if(s[n]%m==0) printf("even n=%d : n | s(n)\n",n); }
  }
  return 0; }
