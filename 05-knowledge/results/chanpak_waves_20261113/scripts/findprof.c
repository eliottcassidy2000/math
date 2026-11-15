// print the first gentourng line of order n whose profile (h, f) equals the target
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
int main(int argc,char**argv){ int n=atoi(argv[1]); i64 th=atoll(argv[2]), tf=atoll(argv[3]); char line[4096]; int adj[16][16]; int N=1<<n;
  i64 *dp=malloc(sizeof(i64)*N*n), *H=malloc(sizeof(i64)*N);
  while(fgets(line,sizeof line,stdin)){ if(line[0]=='>') continue; int k=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ int b=line[k++]=='1'; adj[i][j]=b; adj[j][i]=!b; }
    memset(dp,0,sizeof(i64)*N*n); for(int v=0;v<n;v++) dp[(1<<v)*n+v]=1;
    for(int S=1;S<N;S++) for(int v=0;v<n;v++){ i64 x=dp[S*n+v]; if(!x) continue; for(int u=0;u<n;u++) if(!((S>>u)&1)&&adj[v][u]) dp[(S|1<<u)*n+u]+=x; }
    H[0]=1; for(int S=1;S<N;S++){ i64 q=0; for(int v=0;v<n;v++) q+=dp[S*n+v]; H[S]=q; }
    if(H[N-1]!=th) continue; i64 f=0; for(int S=0;S<N;S++) f+=H[S]*H[(N-1)^S];
    if(f==tf){ fputs(line,stdout); return 0; } }
  return 1; }
