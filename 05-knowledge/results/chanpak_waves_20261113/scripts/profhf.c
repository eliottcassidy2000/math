// Wave 25: (h, f) = (P_T(1), P_T(2)) for every class of order n from gentourng; prints the DISTINCT pairs with class counts.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
#define HS (1u<<26)
static u64 *key; static unsigned *cnt;
static void add(u64 k){ u64 i=(k*0x9E3779B97F4A7C15ULL)>>38; while(key[i] && key[i]!=k) i=(i+1)&(HS-1); key[i]=k; cnt[i]++; }
int main(int argc,char**argv){
  int n=atoi(argv[1]); char line[4096]; int adj[16][16]; int N=1<<n;
  i64 *dp=malloc(sizeof(i64)*N*n), *H=malloc(sizeof(i64)*N);
  key=calloc(HS,sizeof(u64)); cnt=calloc(HS,sizeof(unsigned));
  while(fgets(line,sizeof line,stdin)){
    if(line[0]=='>') continue;
    int k=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(line[k]=='1'){adj[i][j]=1;adj[j][i]=0;} else {adj[j][i]=1;adj[i][j]=0;} k++; }
    memset(dp,0,sizeof(i64)*N*n);
    for(int v=0;v<n;v++) dp[(1<<v)*n+v]=1;
    for(int S=1;S<N;S++) for(int v=0;v<n;v++){ i64 x=dp[S*n+v]; if(!x) continue;
      for(int u=0;u<n;u++) if(!((S>>u)&1) && adj[v][u]) dp[(S|1<<u)*n+u]+=x; }
    H[0]=1; for(int S=1;S<N;S++){ i64 s=0; for(int v=0;v<n;v++) s+=dp[S*n+v]; H[S]=s; }
    i64 f=0; for(int S=0;S<N;S++) f+=H[S]*H[(N-1)^S];
    add(((u64)H[N-1]<<32)|(u64)f);
  }
  for(u64 i=0;i<HS;i++) if(key[i]) printf("%llu %llu %u\n",key[i]>>32,key[i]&0xffffffffULL,cnt[i]);
  return 0; }
