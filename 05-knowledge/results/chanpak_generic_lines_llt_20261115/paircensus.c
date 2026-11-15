// (H, P2) for every tournament read from gentourng (upper triangle, '1' = i->j); prints "H P2" per class.
// P2 = number of ordered pairs of nonempty vertex-disjoint paths covering V = sum_{0<S<V} H(S) H(V\S).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
int main(int argc,char**argv){
  int n=atoi(argv[1]); char line[4096]; unsigned out[32];
  int N=1<<n; u64 *dp=malloc(sizeof(u64)*N*n), *H=malloc(sizeof(u64)*N);
  while(fgets(line,sizeof line,stdin)){
    if(line[0]=='>') continue;
    int k=0; memset(out,0,sizeof out);
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(line[k]=='1') out[i]|=1u<<j; else out[j]|=1u<<i; k++; }
    memset(dp,0,sizeof(u64)*N*n);
    for(int v=0;v<n;v++) dp[(size_t)(1<<v)*n+v]=1;
    H[0]=1;
    for(int S=1;S<N;S++){ u64 tot=0; u64*row=dp+(size_t)S*n;
      unsigned m=S; while(m){ int v=__builtin_ctz(m); m&=m-1; u64 x=row[v]; if(!x) continue; tot+=x;
        unsigned nb=out[v]&~S; while(nb){ int u=__builtin_ctz(nb); nb&=nb-1; dp[(size_t)(S|(1<<u))*n+u]+=x; } }
      H[S]=tot; }
    u64 P2=0; for(int S=1;S<N-1;S++) P2+=H[S]*H[(N-1)^S];
    printf("%llu %llu\n",H[N-1],P2);
  }
  return 0; }
