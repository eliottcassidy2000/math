// One-apex tournament: apex 0, spine 1..n (i -> j for i<j), apex beats spine vertex i iff bit i of word is 1 (O).
// For each single spine-arc flip (i,j), compute the change in H. Report, per word, which even increments in [-G, G]
// occur, and aggregate over random words.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
static int n, V; static unsigned adj[32]; // adj[v] bitmask of out-neighbours
static i64 *dp;
static i64 ham(void){
  int full=(1<<V)-1; memset(dp,0,sizeof(i64)*((size_t)1<<V)*V);
  for(int v=0;v<V;v++) dp[((size_t)1<<v)*V+v]=1;
  for(int m=1;m<=full;m++) for(int v=0;v<V;v++){ i64 c=dp[(size_t)m*V+v]; if(!c) continue; unsigned out=adj[v]&~m; while(out){ int w=__builtin_ctz(out); out&=out-1; dp[((size_t)(m|(1<<w)))*V+w]+=c; } }
  i64 s=0; for(int v=0;v<V;v++) s+=dp[(size_t)full*V+v]; return s;
}
static u64 rng=0x1234567ULL; static u64 xr(void){ rng^=rng<<13; rng^=rng>>7; rng^=rng<<17; return rng; }
int main(int argc,char**argv){
  n=atoi(argv[1]); int trials=atoi(argv[2]); int G=atoi(argv[3]); V=n+1;
  dp=malloc(sizeof(i64)*((size_t)1<<V)*V);
  int *hist=calloc(2*G+1,sizeof(int)); int covered_all=0;
  for(int t=0;t<trials;t++){
    u64 word=xr()&((1ULL<<n)-1);
    for(int v=0;v<V;v++) adj[v]=0;
    for(int i=1;i<=n;i++) for(int j=i+1;j<=n;j++) adj[i]|=1u<<j;
    for(int i=1;i<=n;i++){ if((word>>(i-1))&1) adj[0]|=1u<<i; else adj[i]|=1u; }
    i64 H0=ham();
    char seen[1024]={0};
    for(int i=1;i<=n;i++) for(int j=i+2;j<=n;j++){ // non-adjacent spine flips (adjacent flips just permute)
      adj[i]&=~(1u<<j); adj[j]|=1u<<i; i64 H1=ham(); adj[j]&=~(1u<<i); adj[i]|=1u<<j;
      i64 d=H1-H0; if(d>=-G && d<=G){ seen[d+G]=1; } }
    int all=1; for(int d=-G; d<=G; d+=2) if(d!=0 && !seen[d+G]) all=0;
    if(all) covered_all++;
    for(int d=-G; d<=G; d++) if(seen[d+G]) hist[d+G]++;
    if(t<3){ printf("word %llx H=%lld small increments:",(unsigned long long)word,H0); for(int d=-G; d<=G; d++) if(seen[d+G]) printf(" %d",d); printf("\n"); }
  }
  printf("n=%d, %d random words: fraction of words whose single spine flips realize every nonzero even increment in [-%d,%d]: %.3f\n",n,trials,G,G,(double)covered_all/trials);
  printf("per increment, fraction of words realizing it:"); for(int d=-G; d<=G; d+=2) printf(" %d:%.2f",d,(double)hist[d+G]/trials); printf("\n");
  return 0; }
