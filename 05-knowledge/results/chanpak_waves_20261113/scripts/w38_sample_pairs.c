/* Wave 38: (H, P2) of uniformly random labelled tournaments.
   P2 = number of ordered pairs of nonempty vertex-disjoint paths covering the tournament
      = sum over nonempty proper S of H(T[S]) H(T[V \ S]).
   Usage: w38_sample_pairs n samples seed  ->  lines "H P2"
   Build: cc -O3 -o w38_sample_pairs w38_sample_pairs.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
static uint64_t s[2];
static inline uint64_t rnd(void){ uint64_t a=s[0], b=s[1]; s[0]=b; a^=a<<23; a^=a>>17; a^=b^(b>>26); s[1]=a; return a+b; }
int main(int argc, char **argv){
  int n=atoi(argv[1]); long samples=atol(argv[2]); uint64_t seed=strtoull(argv[3],0,10);
  s[0]=seed*0x9E3779B97F4A7C15ULL+1; s[1]=seed^0xD1B54A32D192ED03ULL; for(int i=0;i<20;i++) rnd();
  int N=1<<n; ll *dp=malloc(sizeof(ll)*(size_t)N*n); ll *cnt=malloc(sizeof(ll)*N); uint32_t out[32];
  for(long t=0;t<samples;t++){
    for(int i=0;i<n;i++) out[i]=0;
    uint64_t r=rnd(); int used=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){
      if(used==64){ r=rnd(); used=0; }
      if(r&1) out[i]|=1u<<j; else out[j]|=1u<<i; r>>=1; used++; }
    memset(dp,0,sizeof(ll)*(size_t)N*n); memset(cnt,0,sizeof(ll)*N);
    for(int v=0;v<n;v++) dp[(size_t)(1<<v)*n+v]=1;
    for(int S=1;S<N;S++) for(int v=0;v<n;v++){ ll c=dp[(size_t)S*n+v]; if(!c) continue; cnt[S]+=c;
      uint32_t nb=out[v]&~S; while(nb){ int u=__builtin_ctz(nb); nb&=nb-1; dp[(size_t)(S|1u<<u)*n+u]+=c; } }
    ll H=cnt[N-1], P2=0; for(int S=1;S<N-1;S++) P2+=cnt[S]*cnt[(N-1)^S];
    printf("%lld %lld\n",H,P2);
  }
  return 0;
}
