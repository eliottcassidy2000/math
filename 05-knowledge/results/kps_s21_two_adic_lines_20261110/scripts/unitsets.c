/* unitsets.c - unit-offset roots through the cut matrix. For each tournament class B (order m, canonical code) compute
 *   H(B), the responses U_B(y), and W_ef = Q(e,f) + Q(f,e) (Q = ordered two-path covers, first ending at e, second starting at f),
 * and list every vertex set O with 2 <= |O| <= KMAX and  sum_{pairs in O} W_ef = H(B) - 1.
 * By THM-5126 (A') the root v of L = B + v with N+(v) = O then has offset 1 and multiplier M = sum_{y in O} U_B(y).
 * usage: unitsets m KMAX < cls_m.txt        output: "m+1 code Omask H(B) M" */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXV 10
static u64 en[1<<MAXV][MAXV], st[1<<MAXV][MAXV], hs[1<<MAXV];
int main(int argc,char**argv){ int m=atoi(argv[1]), KMAX=atoi(argv[2]); char line[1<<16]; int full=(1<<m)-1;
  while(fgets(line,sizeof line,stdin)){ u64 code=strtoull(line,NULL,10); unsigned out[MAXV]={0},inn[MAXV]={0}; int b=0;
    for(int i=0;i<m;i++) for(int j=i+1;j<m;j++){ if(code>>b&1){out[i]|=1u<<j;inn[j]|=1u<<i;} else {out[j]|=1u<<i;inn[i]|=1u<<j;} b++; }
    for(int s=0;s<=full;s++) for(int v=0;v<m;v++){en[s][v]=0;st[s][v]=0;}
    for(int v=0;v<m;v++){en[1<<v][v]=1;st[1<<v][v]=1;}
    for(int s=1;s<=full;s++) for(int v=0;v<m;v++){ u64 c=en[s][v]; if(c){unsigned nx=out[v]&~s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;en[s|1u<<u][u]+=c;}}
      c=st[s][v]; if(c){unsigned nx=inn[v]&~s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;st[s|1u<<u][u]+=c;}} }
    hs[0]=1; for(int s=1;s<=full;s++){u64 t=0; for(int v=0;v<m;v++) t+=en[s][v]; hs[s]=t;}
    i64 H=hs[full], U[MAXV], W[MAXV][MAXV];
    for(int y=0;y<m;y++){ i64 t=0; for(int X=1;X<=full;X++) if((X>>y&1)&&st[X][y]) t+=(i64)st[X][y]*(i64)hs[full^X]; U[y]=t; }
    for(int a=0;a<m;a++) for(int c=0;c<m;c++) W[a][c]=0;
    for(int X=1;X<full;X++){ int Y=full^X; for(int a=0;a<m;a++) if((X>>a&1)&&en[X][a]) for(int c=0;c<m;c++) if((Y>>c&1)&&st[Y][c]){ i64 q=(i64)en[X][a]*(i64)st[Y][c]; W[a][c]+=q; W[c][a]+=q; } }
    for(int O=1;O<=full;O++){ int k=__builtin_popcount(O); if(k<2||k>KMAX) continue; i64 w=0, M=0;
      for(int a=0;a<m;a++) if(O>>a&1){ M+=U[a]; for(int c=a+1;c<m;c++) if(O>>c&1) w+=W[a][c]; }
      if(w==H-1) printf("%d %llu %d %lld %lld\n",m+1,code,O,H,M); } }
  return 0; }
