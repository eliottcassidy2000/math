// read tournaments (n, then n rows of out-masks) from stdin; print H and P(2)=sum_S H(S)H(V\S)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
int main(){
  int n;
  while(scanf("%d",&n)==1){
    unsigned out[32];
    for(int i=0;i<n;i++) scanf("%u",&out[i]);
    int N=1<<n;
    // dp[S][v]: number of paths covering S ending at v
    u64 *dp=calloc((size_t)N*n,sizeof(u64));
    u64 *HS=calloc(N,sizeof(u64));
    for(int v=0;v<n;v++) dp[((size_t)1<<v)*n+v]=1;
    for(int S=1;S<N;S++){
      u64 tot=0;
      for(int v=0;v<n;v++){ u64 c=dp[(size_t)S*n+v]; if(!c) continue; tot+=c;
        unsigned o=out[v]&~S;
        while(o){ int w=__builtin_ctz(o); o&=o-1; dp[((size_t)(S|(1<<w)))*n+w]+=c; } }
      HS[S]=tot;
    }
    HS[0]=1;
    u64 P2=0; for(int S=0;S<N;S++) P2+=HS[S]*HS[(N-1)^S];
    printf("%llu %llu\n",HS[N-1],P2);
    free(dp); free(HS);
  }
  return 0;
}
