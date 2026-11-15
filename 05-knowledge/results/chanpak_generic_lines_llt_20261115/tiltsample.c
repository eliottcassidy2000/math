// Tilted random tournaments: i->j with probability q for i<j. Prints "H P2 seedindex" for each sample.
// usage: tiltsample n q count seed [hlo hhi]   (optional: only print samples with hlo <= H <= hhi)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static u64 s[2];
static u64 rnd(void){ u64 s1=s[0]; const u64 s0=s[1]; s[0]=s0; s1^=s1<<23; s[1]=s1^s0^(s1>>17)^(s0>>26); return s[1]+s0; }
int main(int argc,char**argv){
  int n=atoi(argv[1]); double q=atof(argv[2]); long cnt=atol(argv[3]); u64 seed=strtoull(argv[4],0,10);
  u64 hlo=argc>6?strtoull(argv[5],0,10):0, hhi=argc>6?strtoull(argv[6],0,10):~0ULL;
  s[0]=seed*0x9E3779B97F4A7C15ULL+1; s[1]=seed^0xD1B54A32D192ED03ULL; for(int i=0;i<20;i++) rnd();
  int N=1<<n; u64 *dp=malloc(sizeof(u64)*N*n), *H=malloc(sizeof(u64)*N); unsigned out[32];
  u64 thr=(u64)(q*18446744073709551615.0);
  for(long it=0; it<cnt; it++){
    memset(out,0,sizeof out);
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(rnd()<thr) out[i]|=1u<<j; else out[j]|=1u<<i; }
    memset(dp,0,sizeof(u64)*N*n);
    for(int v=0;v<n;v++) dp[(size_t)(1<<v)*n+v]=1;
    H[0]=1;
    for(int S=1;S<N;S++){ u64 tot=0; u64*row=dp+(size_t)S*n; unsigned m=S;
      while(m){ int v=__builtin_ctz(m); m&=m-1; u64 x=row[v]; if(!x) continue; tot+=x;
        unsigned nb=out[v]&~S; while(nb){ int u=__builtin_ctz(nb); nb&=nb-1; dp[(size_t)(S|(1<<u))*n+u]+=x; } }
      H[S]=tot; }
    if(H[N-1]<hlo||H[N-1]>hhi) continue;
    u64 P2=0; for(int S=1;S<N-1;S++) P2+=H[S]*H[(N-1)^S];
    printf("%llu %llu %ld\n",H[N-1],P2,it);
  }
  return 0; }
