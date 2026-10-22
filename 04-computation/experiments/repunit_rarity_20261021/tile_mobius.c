/* tile_mobius.c -- H as a multilinear polynomial on the tiling cube.  Backbone 0->1->...->n-1; tiles = pairs (i,j),
 * j >= i+2; x_t = 1 means the tile is reversed (back arc j->i).  H(x) for all 2^m tilings, then the Moebius transform
 * c_S = sum_{U subseteq S} (-1)^{|S|-|U|} H(U).  Reports, per |S|, the distribution of v2(c_S) (c_S != 0), and the
 * minimum of v2(c_S) - (number of non-tile-endpoint vertices strictly inside the union of the spans of S).
 * mac-mini-2026-10-21-fanopsp.  Usage: ./tile_mobius n */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n,m; static int ti[64],tj[64];
static int64_t hp(uint64_t mask){ /* number of Hamiltonian paths of the tiling tournament */
  uint32_t out[16]; for(int i=0;i<n;i++){ out[i]=0; for(int j=i+1;j<n;j++) out[i]|=1u<<j; }
  for(int t=0;t<m;t++) if(mask>>t&1){ int i=ti[t],j=tj[t]; out[i]&=~(1u<<j); out[j]|=1u<<i; }
  static int64_t f[1<<12][12]; int full=(1<<n)-1; memset(f,0,sizeof(int64_t)*(full+1)*12);
  for(int v=0;v<n;v++) f[1<<v][v]=1;
  for(int s=1;s<=full;s++) for(int v=0;v<n;v++){ int64_t c=f[s][v]; if(!c) continue; uint32_t cand=out[v]&~s;
      while(cand){int w=__builtin_ctz(cand); cand&=cand-1; f[s|(1<<w)][w]+=c;} }
  int64_t H=0; for(int v=0;v<n;v++) H+=f[full][v]; return H; }
static int v2(int64_t x){ if(x==0) return 99; if(x<0) x=-x; int c=0; while(!(x&1)){x>>=1;c++;} return c; }
int main(int argc,char**argv){ n=atoi(argv[1]); m=0; for(int i=0;i<n;i++) for(int j=i+2;j<n;j++){ ti[m]=i; tj[m]=j; m++; }
  uint64_t N=1ULL<<m; int64_t *c=malloc(N*8);
  for(uint64_t s=0;s<N;s++) c[s]=hp(s);
  for(int t=0;t<m;t++) for(uint64_t s=0;s<N;s++) if(s>>t&1) c[s]-=c[s^(1ULL<<t)];   /* Moebius */
  printf("n=%d tiles=%d  c_empty=%lld\n",n,m,(long long)c[0]);
  /* stats */
  long cnt[32][40]; memset(cnt,0,sizeof cnt); int minexcess=99; long nonzero=0; long minex_count=0;
  for(uint64_t s=1;s<N;s++){ if(c[s]==0) continue; nonzero++; int k=__builtin_popcountll(s); int v=v2(c[s]); if(v<40) cnt[k][v]++;
    /* free interior vertices: inside [min i, max j] of the tiles, not an endpoint of any tile in S */
    int lo=n,hi=-1; uint32_t ends=0; for(int t=0;t<m;t++) if(s>>t&1){ if(ti[t]<lo) lo=ti[t]; if(tj[t]>hi) hi=tj[t]; ends|=1u<<ti[t]; ends|=1u<<tj[t]; }
    int fr=0; for(int x=lo+1;x<hi;x++) if(!(ends>>x&1)) fr++;
    int ex=v-fr; if(ex<minexcess){ minexcess=ex; minex_count=0; } if(ex==minexcess) minex_count++; }
  printf("nonzero Moebius coefficients: %ld of %llu\n",nonzero,(unsigned long long)(N-1));
  for(int k=1;k<=m;k++){ long tot=0; for(int v=0;v<40;v++) tot+=cnt[k][v]; if(!tot) continue; printf(" |S|=%2d: ",k); for(int v=0;v<16;v++) if(cnt[k][v]) printf("v2=%d:%ld ",v,cnt[k][v]); printf("\n"); }
  printf("min over S of v2(c_S) - #free interior vertices = %d (attained %ld times)\n",minexcess,minex_count);
  return 0; }
