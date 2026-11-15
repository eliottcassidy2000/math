// Wave 25: path-cover profile census. Reads gentourng output for order n (upper triangle, '1' = i -> j),
// prints  h f F3  per class, with h = H(T) = P_T(1), f = F_2 = P_T(2), F3 = P_T(3) (ordered q-slot covers,
// empty slots allowed), computed by subset DP: H(S) for all S, F2(R) = sum_{S<=R} H(S)H(R\S), F3 = sum_S H(S)F2(V\S).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
int main(int argc,char**argv){
  int n=atoi(argv[1]); char line[4096]; int adj[16][16];
  int N=1<<n; i64 *dp=malloc(sizeof(i64)*N*n), *H=malloc(sizeof(i64)*N), *F2=malloc(sizeof(i64)*N);
  while(fgets(line,sizeof line,stdin)){
    if(line[0]=='>') continue;
    int k=0; memset(adj,0,sizeof adj);
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(line[k]=='1') adj[i][j]=1; else adj[j][i]=1; k++; }
    memset(dp,0,sizeof(i64)*N*n);
    for(int v=0;v<n;v++) dp[(1<<v)*n+v]=1;
    for(int S=1;S<N;S++) for(int v=0;v<n;v++){ i64 x=dp[S*n+v]; if(!x) continue;
      for(int u=0;u<n;u++) if(!((S>>u)&1) && adj[v][u]) dp[(S|1<<u)*n+u]+=x; }
    H[0]=1; for(int S=1;S<N;S++){ i64 s=0; for(int v=0;v<n;v++) s+=dp[S*n+v]; H[S]=s; }
    for(int R=0;R<N;R++){ i64 s=0; for(int S=R;;S=(S-1)&R){ s+=H[S]*H[R^S]; if(!S) break; } F2[R]=s; }
    int full=N-1; i64 F3=0; for(int S=0;S<N;S++) F3+=H[S]*F2[full^S];
    printf("%lld %lld %lld\n",H[full],F2[full],F3);
  }
  return 0; }
