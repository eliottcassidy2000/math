#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
int main(int argc,char**argv){ int n=atoi(argv[1]); char line[256]; int N=1<<n;
  ll *dp=malloc(sizeof(ll)*N*n); ll *cnt=malloc(sizeof(ll)*N); uint32_t out[32];
  while(fgets(line,sizeof line,stdin)){ if(line[0]=='>') continue; int k=0;
    for(int i=0;i<n;i++) out[i]=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(line[k]=='1') out[i]|=1u<<j; else out[j]|=1u<<i; k++; }
    memset(dp,0,sizeof(ll)*N*n); memset(cnt,0,sizeof(ll)*N);
    for(int v=0;v<n;v++) dp[(size_t)(1<<v)*n+v]=1;
    for(int S=1;S<N;S++) for(int v=0;v<n;v++){ ll c=dp[(size_t)S*n+v]; if(!c) continue; cnt[S]+=c;
       uint32_t nb=out[v]&~S; while(nb){int u=__builtin_ctz(nb); nb&=nb-1; dp[(size_t)(S|1u<<u)*n+u]+=c;} }
    ll H=cnt[N-1], P2=0; for(int S=1;S<N-1;S++) P2+=cnt[S]*cnt[(N-1)^S];
    printf("%lld %lld\n",H,P2); }
  return 0; }
