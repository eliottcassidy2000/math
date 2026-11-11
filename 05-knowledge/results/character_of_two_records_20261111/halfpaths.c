// halfpaths.c - invariant odd cycles of a tournament with an anti-involution psi having exactly one fixed vertex v0.
//
// An odd cycle C with psi(C) = C passes through v0 and has the form
//     v0 -> a_1 -> ... -> a_k -> psi(a_k) -> ... -> psi(a_1) -> v0,
// where (a_1, ..., a_k) is a directed path that uses at most one vertex from each pair {a, psi(a)}, v0 -> a_1 and
// a_k -> psi(a_k).  This program counts these "half-paths" by k with a dynamic programme over
// (set of pairs used, last vertex), counts modulo 256.
//
// Input (stdin): m (number of pairs), then 2m lines "out_mask start", vertex 2i and 2i+1 forming pair i; out_mask is a
// decimal 64-bit mask of out-neighbours among the 2m vertices; start = 1 iff v0 -> vertex.
// Output: for k = 1..m: F_k mod 256, split by the type of the half-path:
//   type (i)  a_k -> v0  (then the triple {v0, a_k, psi(a_k)} is transitive),
//   type (ii) v0 -> a_k  (then it is cyclic);
// and the parity of F = sum_k F_k.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef unsigned char u8;
int main(void){ int m; if(scanf("%d",&m)!=1) return 1; int n=2*m; u64 out[64]; int st[64];
  for(int v=0;v<n;v++){ if(scanf("%llu %d",&out[v],&st[v])!=2) return 1; }
  size_t NM=(size_t)1<<m; u8 *dp=calloc(NM*(size_t)n,1); if(!dp){fprintf(stderr,"memory\n");return 1;}
  for(int v=0;v<n;v++) if(st[v]) dp[((size_t)1<<(v>>1))*n+v]=1;
  unsigned Fi[64], Fii[64]; memset(Fi,0,sizeof Fi); memset(Fii,0,sizeof Fii);
  for(size_t mask=1;mask<NM;mask++){ int k=__builtin_popcountll(mask); u8 *row=dp+mask*n;
    for(int v=0;v<n;v++){ u8 c=row[v]; if(!c) continue;
      int closes=(out[v]>>(v^1))&1;          // v -> psi(v)
      if(closes){ if(st[v]) Fii[k]+=c; else Fi[k]+=c; }
      u64 nb=out[v];
      for(int j=0;j<m;j++){ if(mask>>j&1) continue; size_t nm=mask|((size_t)1<<j);
        if(nb>>(2*j)&1) dp[nm*n+2*j]+=c;
        if(nb>>(2*j+1)&1) dp[nm*n+2*j+1]+=c; } } }
  unsigned tot=0;
  for(int k=1;k<=m;k++){ printf("k %d length %d F_k mod 256 = %u (type i %u, type ii %u)\n", k, 2*k+1, (Fi[k]+Fii[k])&255, Fi[k]&255, Fii[k]&255); tot+=Fi[k]+Fii[k]; }
  printf("F mod 2 = %u ; F_1 mod 256 = %u ; sum_{k>=2} F_k mod 2 = %u\n", tot&1, (Fi[1]+Fii[1])&255, (tot-(Fi[1]+Fii[1]))&1);
  free(dp); return 0; }
