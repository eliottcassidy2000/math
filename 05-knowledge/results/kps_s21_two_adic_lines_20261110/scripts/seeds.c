/* seeds.c - unit-offset roots. For each tournament class K (canonical code, order m) and each root v compute
 *   H = H(K), d = H(K - v), E = # Hamiltonian paths ending at v, U = # covers (P1 possibly empty, P2 starting at v),
 * and print the roots with c0 = H - E + 2d - U = 1:   "m code v H d E U M"  with M = U - d.
 * For these roots the tail law gives H(K o 0^t) = 1 + M 2^t for every t >= 1 (new chain vertices beat v and lose to the rest).
 * usage: seeds m < cls_m.txt        (with a second argument "all": print every root, c0 included) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXV 11
static u64 en[1<<MAXV][MAXV], st[1<<MAXV][MAXV], hs[1<<MAXV];
int main(int argc,char**argv){ int m=atoi(argv[1]); int all=(argc>2); char line[1<<16]; int full=(1<<m)-1;
  while(fgets(line,sizeof line,stdin)){ u64 code=strtoull(line,NULL,10); unsigned out[MAXV]={0},inn[MAXV]={0}; int b=0;
    for(int i=0;i<m;i++) for(int j=i+1;j<m;j++){ if(code>>b&1){out[i]|=1u<<j;inn[j]|=1u<<i;} else {out[j]|=1u<<i;inn[i]|=1u<<j;} b++; }
    for(int s=0;s<=full;s++) for(int v=0;v<m;v++){en[s][v]=0;st[s][v]=0;}
    for(int v=0;v<m;v++){en[1<<v][v]=1;st[1<<v][v]=1;}
    for(int s=1;s<=full;s++) for(int v=0;v<m;v++){ u64 c=en[s][v]; if(c){unsigned nx=out[v]&~s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;en[s|1u<<u][u]+=c;}}
      c=st[s][v]; if(c){unsigned nx=inn[v]&~s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;st[s|1u<<u][u]+=c;}} }
    hs[0]=1; for(int s=1;s<=full;s++){u64 t=0; for(int v=0;v<m;v++) t+=en[s][v]; hs[s]=t;}
    i64 H=hs[full];
    for(int v=0;v<m;v++){ i64 d=hs[full^(1<<v)], E=en[full][v], U=0;
      for(int X=1;X<=full;X++) if((X>>v&1)&&st[X][v]) U+=(i64)st[X][v]*(i64)hs[full^X];
      i64 c0=H-E+2*d-U;
      if(all) printf("%d %llu %d %lld %lld %lld %lld %lld %lld\n",m,code,v,H,d,E,U,U-d,c0);
      else if(c0==1) printf("%d %llu %d %lld %lld %lld %lld %lld\n",m,code,v,H,d,E,U,U-d); } }
  return 0; }
