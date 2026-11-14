// For random one-apex words (spine n), over all single non-adjacent spine-arc flips: the smallest nonzero |change| of H,
// as a fraction of H, and the number of flips with |change| <= 20. Aggregated.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef long long i64; typedef unsigned long long u64;
static int n, V; static unsigned adj[32]; static i64 *dp;
static i64 ham(void){ int full=(1<<V)-1; memset(dp,0,sizeof(i64)*((size_t)1<<V)*V); for(int v=0;v<V;v++) dp[((size_t)1<<v)*V+v]=1;
  for(int m=1;m<=full;m++) for(int v=0;v<V;v++){ i64 c=dp[(size_t)m*V+v]; if(!c) continue; unsigned out=adj[v]&~m; while(out){ int w=__builtin_ctz(out); out&=out-1; dp[((size_t)(m|(1<<w)))*V+w]+=c; } }
  i64 s=0; for(int v=0;v<V;v++) s+=dp[(size_t)full*V+v]; return s; }
static u64 rng=0x9876543ULL; static u64 xr(void){ rng^=rng<<13; rng^=rng>>7; rng^=rng<<17; return rng; }
int main(int argc,char**argv){
  n=atoi(argv[1]); int trials=atoi(argv[2]); V=n+1; dp=malloc(sizeof(i64)*((size_t)1<<V)*V);
  double sumfrac=0; int words_small=0; long small_total=0, flips_total=0; double minfrac_overall=1e9;
  for(int t=0;t<trials;t++){ u64 word=xr()&((1ULL<<n)-1);
    for(int v=0;v<V;v++) adj[v]=0; for(int i=1;i<=n;i++) for(int j=i+1;j<=n;j++) adj[i]|=1u<<j;
    for(int i=1;i<=n;i++){ if((word>>(i-1))&1) adj[0]|=1u<<i; else adj[i]|=1u; }
    i64 H0=ham(); i64 mind=-1; int nsmall=0;
    for(int i=1;i<=n;i++) for(int j=i+2;j<=n;j++){ adj[i]&=~(1u<<j); adj[j]|=1u<<i; i64 d=ham()-H0; adj[j]&=~(1u<<i); adj[i]|=1u<<j;
      flips_total++; i64 a=d<0?-d:d; if(a && (mind<0||a<mind)) mind=a; if(a && a<=20) nsmall++; }
    if(nsmall) words_small++; small_total+=nsmall; double f=(double)mind/H0; sumfrac+=f; if(f<minfrac_overall) minfrac_overall=f; }
  printf("spine n=%d, %d random one-apex words: mean (smallest nonzero |dH|)/H = %.3f, min over words %.4f; words with some flip |dH|<=20: %.2f; flips with |dH|<=20: %.4f of all flips\n",
    n,trials,sumfrac/trials,minfrac_overall,(double)words_small/trials,(double)small_total/flips_total);
  return 0; }
