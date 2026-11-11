// hampaths_exact.c - exact number of Hamiltonian paths of tournaments on N <= 22 vertices (subset dynamic programme
// over (set, last vertex) with 64-bit counts; the count must be below 2^64, true for N <= 22 since H <= N!/2^(N-1) * e
// roughly and 22!/2^21 < 2^64).
// Input (stdin): any number of blocks "N" followed by N decimal out-neighbour masks.
// Output: one line per block: "N H" (H exact).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
int main(void){ int N; u64 *dp=NULL; size_t cap=0;
  while(scanf("%d",&N)==1){ u64 out[32]; for(int v=0;v<N;v++) if(scanf("%llu",&out[v])!=1) return 1;
    if(N>22){ fprintf(stderr,"N too large\n"); return 1; }
    size_t NM=(size_t)1<<N, need=NM*(size_t)N; if(need>cap){ free(dp); dp=malloc(need*sizeof(u64)); cap=need; if(!dp){fprintf(stderr,"memory\n");return 1;} }
    memset(dp,0,need*sizeof(u64));
    for(int v=0;v<N;v++) dp[((size_t)1<<v)*N+v]=1;
    for(size_t mask=1;mask<NM;mask++){ u64 *row=dp+mask*N;
      for(int v=0;v<N;v++){ u64 c=row[v]; if(!c) continue; u64 nb=out[v]&~(u64)mask&(NM-1);
        while(nb){ int w=__builtin_ctzll(nb); nb&=nb-1; dp[(mask|((size_t)1<<w))*N+w]+=c; } } }
    u64 H=0; for(int v=0;v<N;v++) H+=dp[(NM-1)*N+v];
    printf("%d %llu\n",N,H); fflush(stdout); }
  free(dp); return 0; }
