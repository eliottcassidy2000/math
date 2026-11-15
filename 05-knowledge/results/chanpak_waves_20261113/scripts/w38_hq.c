/* Wave 38: exact H and Q = Pi_2 for tournaments given by out-neighbour masks.
   Input lines:  n r_0 r_1 ... r_{n-1}   (r_i = bitmask of the vertices that i beats)
   Output lines: H Q   with Q = sum over all S (empty and full included) of H(T[S]) H(T[V \ S]), H(empty) = 1,
   i.e. Q = P2 + 2H counts ordered pairs of possibly empty vertex-disjoint paths covering T.
   The triangle reader of chanpak_pair_line_census is F = H + Q = H(C3[T, *, *]).
   Build: cc -O3 -o w38_hq w38_hq.c      (n <= 21) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
int main(void){
  int n; static uint32_t out[32]; ll *dp=0, *cnt=0; int cap=0;
  while(scanf("%d",&n)==1){
    for(int i=0;i<n;i++){ unsigned x; if(scanf("%u",&x)!=1) return 1; out[i]=x; }
    int N=1<<n;
    if(n>cap){ free(dp); free(cnt); dp=malloc(sizeof(ll)*(size_t)N*n); cnt=malloc(sizeof(ll)*N); cap=n; }
    memset(dp,0,sizeof(ll)*(size_t)N*n); memset(cnt,0,sizeof(ll)*N);
    for(int v=0;v<n;v++) dp[(size_t)(1<<v)*n+v]=1;
    cnt[0]=1;
    for(int S=1;S<N;S++) for(int v=0;v<n;v++){ ll c=dp[(size_t)S*n+v]; if(!c) continue; cnt[S]+=c;
      uint32_t nb=out[v]&~S; while(nb){ int u=__builtin_ctz(nb); nb&=nb-1; dp[(size_t)(S|1u<<u)*n+u]+=c; } }
    ll H=cnt[N-1], Q=0; for(int S=0;S<N;S++) Q+=cnt[S]*cnt[(N-1)^S];
    printf("%lld %lld\n",H,Q); fflush(stdout);
  }
  return 0;
}
