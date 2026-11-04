// twohc.c: for every tournament on stdin (gentourng format, order n) count Hamiltonian cycles via subset DP
// (paths from vertex 0 over subsets containing 0; close with an arc back to 0). Report the numbers of tournaments
// with exactly 1 and exactly 2 Hamiltonian cycles, and for the exactly-2 ones the number of out-degree-1 vertices.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n; static u64 adj[16]; static u64 *dp;
int main(int argc,char**argv){ n=atoi(argv[1]); u64 N=1ULL<<n; dp=malloc(sizeof(u64)*N*n);
  char line[256]; u64 total=0,one=0,two=0; u64 od1[16]={0};
  while(fgets(line,sizeof line,stdin)){ int L=strlen(line); while(L&&(line[L-1]=='\n'||line[L-1]=='\r')) line[--L]=0; if(L!=n*(n-1)/2) continue;
    memset(adj,0,sizeof adj); int idx=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(line[idx++]=='1') adj[i]|=1ULL<<j; else adj[j]|=1ULL<<i; }
    total++;
    // dp[S][v]: # paths starting at 0, visiting exactly S (S contains 0), ending at v
    memset(dp,0,sizeof(u64)*N*n); dp[1*n+0]=1;
    for(u64 S=1;S<N;S+=2){ for(int v=0;v<n;v++){ u64 c=dp[S*n+v]; if(!c) continue; u64 nx=adj[v]&~S; while(nx){int u=__builtin_ctzll(nx); nx&=nx-1; dp[(S|(1ULL<<u))*n+u]+=c;} } }
    u64 cyc=0; for(int v=1;v<n;v++) if(adj[v]&1ULL) cyc+=dp[(N-1)*n+v];
    if(cyc==1) one++; else if(cyc==2){ two++; int k=0; for(int v=0;v<n;v++) if(__builtin_popcountll(adj[v])==1) k++; od1[k]++; }
  }
  printf("order %d: %llu tournaments; exactly one Hamiltonian cycle: %llu; exactly two: %llu; out-degree-1 vertex counts among the two-cycle ones:",n,total,one,two);
  for(int k=0;k<n;k++) if(od1[k]) printf(" %d:%llu",k,od1[k]); printf("\n"); return 0; }
