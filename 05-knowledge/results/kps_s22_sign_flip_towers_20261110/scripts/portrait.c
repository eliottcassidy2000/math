/* portrait.c - states of all pairs (B, O): B a tournament class of order m (code = first field of each input line), O any vertex set.
 *   h = H(B), r = h - W_B(O), M = U_B(O), E = End_B(O).
 * usage: portrait m RMAX < cls_m.txt      prints every pair with |r| <= RMAX as  "m code Omask r h M E"  */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXV 10
static u64 en[1<<MAXV][MAXV], st[1<<MAXV][MAXV], hs[1<<MAXV];
int main(int argc,char**argv){ int m=atoi(argv[1]); i64 RMAX=atoll(argv[2]); char line[1<<16]; int full=(1<<m)-1;
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
    for(int O=0;O<=full;O++){ i64 w=0, M=0, E=0;
      for(int a=0;a<m;a++) if(O>>a&1){ M+=U[a]; E+=en[full][a]; for(int c=0;c<m;c++) if(c!=a&&(O>>c&1)) w+=Q[a][c]; }
      i64 r=H-w; if(r<=RMAX&&r>=-RMAX) printf("%d %llu %d %lld %lld %lld %lld\n",m,code,O,r,H,M,E); } }
  return 0; }
