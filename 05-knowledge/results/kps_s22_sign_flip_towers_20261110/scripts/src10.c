/* src10.c - selected sets of any size with a Fermat-type offset: prints "m code Omask r h M" when r - 1 is a power of two, r >= 3,
 * and the odd part of M equals D (argv[2]).   usage: src10 m D < cls_m.txt */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXV 10
static u64 en[1<<MAXV][MAXV], st[1<<MAXV][MAXV], hs[1<<MAXV];
int main(int argc,char**argv){ int m=atoi(argv[1]); i64 DD=atoll(argv[2]); char line[1<<16]; int full=(1<<m)-1;
  while(fgets(line,sizeof line,stdin)){ u64 code=strtoull(line,NULL,10); unsigned out[MAXV]={0},inn[MAXV]={0}; int b=0;
    for(int i=0;i<m;i++) for(int j=i+1;j<m;j++){ if(code>>b&1){out[i]|=1u<<j;inn[j]|=1u<<i;} else {out[j]|=1u<<i;inn[i]|=1u<<j;} b++; }
    for(int s=0;s<=full;s++) for(int v=0;v<m;v++){en[s][v]=0;st[s][v]=0;}
    for(int v=0;v<m;v++){en[1<<v][v]=1;st[1<<v][v]=1;}
    for(int s=1;s<=full;s++) for(int v=0;v<m;v++){ u64 c=en[s][v]; if(c){unsigned nx=out[v]&~s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;en[s|1u<<u][u]+=c;}}
      c=st[s][v]; if(c){unsigned nx=inn[v]&~s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;st[s|1u<<u][u]+=c;}} }
    hs[0]=1; for(int s=1;s<=full;s++){u64 t=0; for(int v=0;v<m;v++) t+=en[s][v]; hs[s]=t;}
    i64 H=hs[full], U[MAXV], Q[MAXV][MAXV];
    for(int y=0;y<m;y++){ i64 t=0; for(int X=1;X<=full;X++) if((X>>y&1)&&st[X][y]) t+=(i64)st[X][y]*(i64)hs[full^X]; U[y]=t; }
    for(int a=0;a<m;a++) for(int c=0;c<m;c++) Q[a][c]=0;
    for(int X=1;X<full;X++){ int Y=full^X; for(int a=0;a<m;a++) if((X>>a&1)&&en[X][a]) for(int c=0;c<m;c++) if((Y>>c&1)&&st[Y][c]) Q[a][c]+=(i64)en[X][a]*(i64)st[Y][c]; }
    for(int O=1;O<=full;O++){ i64 w=0, M=0; for(int a=0;a<m;a++) if(O>>a&1){ M+=U[a]; for(int c=0;c<m;c++) if(c!=a&&(O>>c&1)) w+=Q[a][c]; }
      i64 r=H-w; if(r<3) continue; i64 x=r-1; if(x&(x-1)) continue; i64 c=M; while(c%2==0) c/=2; if(c==DD) printf("%d %llu %d %lld %lld %lld\n",m,code,O,r,H,M); } }
  return 0; }
