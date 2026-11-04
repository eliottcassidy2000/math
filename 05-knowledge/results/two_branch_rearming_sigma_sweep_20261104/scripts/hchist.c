// hchist.c: histogram of the number of Hamiltonian cycles (0..K) over all tournaments of order n (gentourng input),
// plus, for each count k <= K, the number of those tournaments having at least one vertex of out-degree 1.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define K 24
static int n; static u64 adj[16]; static u64 *dp;
int main(int argc,char**argv){ n=atoi(argv[1]); u64 N=1ULL<<n; dp=malloc(sizeof(u64)*N*n);
  char line[256]; u64 total=0, hist[K+2]={0}, hist1[K+2]={0};
  while(fgets(line,sizeof line,stdin)){ int L=strlen(line); while(L&&(line[L-1]=='\n'||line[L-1]=='\r')) line[--L]=0; if(L!=n*(n-1)/2) continue;
    memset(adj,0,sizeof adj); int idx=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(line[idx++]=='1') adj[i]|=1ULL<<j; else adj[j]|=1ULL<<i; }
    total++;
    memset(dp,0,sizeof(u64)*N*n); dp[1*n+0]=1;
    for(u64 S=1;S<N;S+=2){ for(int v=0;v<n;v++){ u64 c=dp[S*n+v]; if(!c) continue; u64 nx=adj[v]&~S; while(nx){int u=__builtin_ctzll(nx); nx&=nx-1; dp[(S|(1ULL<<u))*n+u]+=c;} } }
    u64 cyc=0; for(int v=1;v<n;v++) if(adj[v]&1ULL) cyc+=dp[(N-1)*n+v];
    int k = cyc>K ? K+1 : (int)cyc; hist[k]++;
    int has1=0; for(int v=0;v<n;v++) if(__builtin_popcountll(adj[v])==1) has1=1; if(has1) hist1[k]++;
  }
  printf("order %d: %llu tournaments; #HC histogram k:count(with an out-degree-1 vertex):",n,total);
  for(int k=0;k<=K;k++) printf(" %d:%llu(%llu)",k,hist[k],hist1[k]); printf(" >%d:%llu\n",K,hist[K+1]); return 0; }
