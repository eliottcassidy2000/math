/* singles.c - one-vertex out-sets: for every tournament class B (order m, code = first field of each input line) and vertex o,
 * h = H(B), U = U_B(o) (covers (P1 possibly empty, P2 starting at o)), E = paths ending at o.
 * Prints "m code o h U E" when h - 1 is a power of two and h >= 3 (sources with a Fermat-type offset).   usage: singles m < cls_m.txt */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXV 10
static u64 en[1<<MAXV][MAXV], st[1<<MAXV][MAXV], hs[1<<MAXV];
int main(int argc,char**argv){ int m=atoi(argv[1]); char line[1<<16]; int full=(1<<m)-1;
  while(fgets(line,sizeof line,stdin)){ u64 code=strtoull(line,NULL,10); unsigned out[MAXV]={0},inn[MAXV]={0}; int b=0;
    for(int i=0;i<m;i++) for(int j=i+1;j<m;j++){ if(code>>b&1){out[i]|=1u<<j;inn[j]|=1u<<i;} else {out[j]|=1u<<i;inn[i]|=1u<<j;} b++; }
    for(int s=0;s<=full;s++) for(int v=0;v<m;v++){en[s][v]=0;st[s][v]=0;}
    for(int v=0;v<m;v++){en[1<<v][v]=1;st[1<<v][v]=1;}
    for(int s=1;s<=full;s++) for(int v=0;v<m;v++){ u64 c=en[s][v]; if(c){unsigned nx=out[v]&~s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;en[s|1u<<u][u]+=c;}}
      c=st[s][v]; if(c){unsigned nx=inn[v]&~s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;st[s|1u<<u][u]+=c;}} }
    hs[0]=1; for(int s=1;s<=full;s++){u64 t=0; for(int v=0;v<m;v++) t+=en[s][v]; hs[s]=t;}
    i64 H=hs[full]; i64 x=H-1; if(H<3 || (x&(x-1))) continue;
    for(int y=0;y<m;y++){ i64 t=0; for(int X=1;X<=full;X++) if((X>>y&1)&&st[X][y]) t+=(i64)st[X][y]*(i64)hs[full^X];
      printf("%d %llu %d %lld %lld %llu\n",m,code,y,H,t,en[full][y]); } }
  return 0; }
