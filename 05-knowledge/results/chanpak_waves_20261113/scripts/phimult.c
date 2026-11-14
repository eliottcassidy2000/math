// Max multiplicity of Phi_V(W) = 917 W - 9 Psi(W) - 9 V rev_m(W) - 9 (v0 * w_{m-1}) - 9 w_0 over m-bit words W,
// for V in [0, VMAX] (Psi = binary pair energy: g(W) = 1 + W - w0 - Psi(W)). Also the level-m energy sum_V E_m(V).
#include <stdio.h>
#include <stdlib.h>
typedef long long i64; typedef unsigned long long u64;
static i64 g(i64 U){ if(U==0) return 1; int k=64-__builtin_clzll((u64)U); i64 s=1;
  for(int t=0;t<k;t++){ int bt=(U>>t)&1; if(!bt) s+=U>>(t+1); if(t+1<k && ((U>>(t+1))&1) && !bt) s+=1; } return s; }
static int cmp(const void*a,const void*b){ i64 x=*(const i64*)a,y=*(const i64*)b; return x<y?-1:x>y; }
int main(int argc,char**argv){
  int m=atoi(argv[1]); i64 VMAX=atoll(argv[2]); i64 VSTEP=argc>3?atoll(argv[3]):1;
  i64 M=1LL<<m; i64 *psi=malloc(sizeof(i64)*M), *rv=malloc(sizeof(i64)*M), *phi=malloc(sizeof(i64)*M);
  for(i64 W=0;W<M;W++){ i64 w0=W&1; psi[W]=1+W-w0-g(W); i64 r=0; for(int t=0;t<m;t++) if((W>>t)&1) r|=1LL<<(m-1-t); rv[W]=r; }
  i64 worst=0, worstV=0; double energy=0; i64 nV=0;
  for(i64 V=0;V<=VMAX;V+=VSTEP){ nV++;
    for(i64 W=0;W<M;W++){ i64 w0=W&1, wtop=(W>>(m-1))&1; phi[W]=917*W-9*psi[W]-9*V*rv[W]-9*((V&1)*wtop)-9*w0; }
    qsort(phi,M,sizeof(i64),cmp);
    i64 run=1, mx=1; double e=0;
    for(i64 i=1;i<=M;i++){ if(i<M && phi[i]==phi[i-1]) run++; else { if(run>mx) mx=run; e+=(double)run*(run-1); run=1; } }
    energy+=e; if(mx>worst){ worst=mx; worstV=V; } }
  printf("m=%d, V in [0,%lld] step %lld (%lld values): max multiplicity of Phi_V %lld (at V=%lld); mean off-diagonal energy per V %.2f (vs 2^m=%lld)\n",m,VMAX,VSTEP,nV,worst,worstV,energy/nV,M);
  return 0; }
