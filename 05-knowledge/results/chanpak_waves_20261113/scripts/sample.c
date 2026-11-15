// Wave 25 (P73): (h, f) = (P_T(1), P_T(2)) of uniform random labeled tournaments of order n. Usage: sample n M seed
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long i64; typedef unsigned long long u64;
static u64 s; static u64 rnd(void){ s^=s<<13; s^=s>>7; s^=s<<17; return s; }
int main(int argc,char**argv){ int n=atoi(argv[1]); long M=atol(argv[2]); s=strtoull(argv[3],0,10)|1; int N=1<<n; int adj[16][16];
  i64 *dp=malloc(sizeof(i64)*N*n), *H=malloc(sizeof(i64)*N);
  for(long t=0;t<M;t++){
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ int b=rnd()>>63; adj[i][j]=b; adj[j][i]=!b; }
    memset(dp,0,sizeof(i64)*N*n); for(int v=0;v<n;v++) dp[(1<<v)*n+v]=1;
    for(int S=1;S<N;S++) for(int v=0;v<n;v++){ i64 x=dp[S*n+v]; if(!x) continue; for(int u=0;u<n;u++) if(!((S>>u)&1)&&adj[v][u]) dp[(S|1<<u)*n+u]+=x; }
    H[0]=1; for(int S=1;S<N;S++){ i64 q=0; for(int v=0;v<n;v++) q+=dp[S*n+v]; H[S]=q; }
    i64 f=0; for(int S=0;S<N;S++) f+=H[S]*H[(N-1)^S];
    printf("%lld %lld\n",H[N-1],f); }
  return 0; }
