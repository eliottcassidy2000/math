/* carriers10.c - ports of unit pairs. Input lines "L codeB Omask H M" (unitsets.c output; base B of order m = L - 1).
 * For every vertex u of B outside O with in-degree exactly 1 in B: replace u by TT_2 (u1 -> u2) and recompute the state of (B2, O).
 * If the offset is still 1, print "m codeB Omask u M1 M2".  (A module at such a port is visited at most twice, so the slope family is
 * M_a = (M1 - 2 de) + de 2^a with de = (M2 - M1)/2, and the unit offset holds for every a.)   usage: carriers10 < us_m.txt */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXV 11
static u64 en[1<<MAXV][MAXV], st[1<<MAXV][MAXV], hs[1<<MAXV];
static int A[MAXV][MAXV];
static void state(int n, unsigned O, i64*r, i64*M){ int full=(1<<n)-1; unsigned out[MAXV],inn[MAXV];
  for(int i=0;i<n;i++){out[i]=0;inn[i]=0;for(int j=0;j<n;j++){if(A[i][j])out[i]|=1u<<j; if(A[j][i])inn[i]|=1u<<j;}}
  for(int s=0;s<=full;s++) for(int v=0;v<n;v++){en[s][v]=0;st[s][v]=0;}
  for(int v=0;v<n;v++){en[1<<v][v]=1;st[1<<v][v]=1;}
  for(int s=1;s<=full;s++) for(int v=0;v<n;v++){ u64 c=en[s][v]; if(c){unsigned nx=out[v]&~s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;en[s|1u<<u][u]+=c;}}
    c=st[s][v]; if(c){unsigned nx=inn[v]&~s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;st[s|1u<<u][u]+=c;}} }
  hs[0]=1; for(int s=1;s<=full;s++){u64 t=0; for(int v=0;v<n;v++) t+=en[s][v]; hs[s]=t;}
  i64 W=0,U=0;
  for(int o=0;o<n;o++) if(O>>o&1) for(int X=1;X<=full;X++) if((X>>o&1)&&st[X][o]) U+=(i64)st[X][o]*(i64)hs[full^X];
  for(int X=1;X<full;X++){ int Y=full^X; i64 a=0,b=0; for(int e=0;e<n;e++) if(O>>e&1){ if(X>>e&1) a+=en[X][e]; else b+=st[Y][e]; } W+=a*b; }
  *r=(i64)hs[full]-W; *M=U; }
int main(void){ char line[1<<12]; int K[MAXV][MAXV];
  while(fgets(line,sizeof line,stdin)){ int L; u64 code; unsigned O; i64 h0,M0; if(sscanf(line,"%d %llu %u %lld %lld",&L,&code,&O,&h0,&M0)!=5) continue;
    int m=L-1; int b=0; for(int i=0;i<m;i++) for(int j=i+1;j<m;j++){ K[i][j]=(code>>b&1); K[j][i]=!K[i][j]; b++; } for(int i=0;i<m;i++) K[i][i]=0;
    for(int u=0;u<m;u++){ if(O>>u&1) continue; int indeg=0; for(int x=0;x<m;x++) indeg+=K[x][u]; if(indeg!=1) continue;
      /* B2: vertices 0..m-1 as in B, plus vertex m = second twin of u (u -> m) */
      int n=m+1; for(int i=0;i<m;i++) for(int j=0;j<m;j++) A[i][j]=K[i][j];
      for(int i=0;i<m;i++){ if(i==u){A[u][m]=1;A[m][u]=0;} else {A[i][m]=K[i][u];A[m][i]=K[u][i];} } A[m][m]=0;
      i64 r,M; state(n,O,&r,&M); if(r==1) printf("%d %llu %u %d %lld %lld\n",m,code,O,u,M0,M); } }
  return 0; }
