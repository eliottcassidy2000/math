// hampaths_layered.c - number of Hamiltonian paths of a tournament on N <= 26 vertices that start in a given set,
// modulo 256, by a subset dynamic programme stored layer by layer (two layers of C(N,k) * N bytes and a rank table
// of 2^N 32-bit integers).
// Input (stdin): N, then N lines "out_mask start" (decimal out-neighbour mask; start = 1 if a path may start there).
// Output: the number of Hamiltonian paths starting in the start set, mod 256, and the same count by end vertex.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef unsigned char u8; typedef unsigned int u32;
int main(void){ int N; if(scanf("%d",&N)!=1) return 1; u64 out[32]; int st[32];
  for(int v=0;v<N;v++) if(scanf("%llu %d",&out[v],&st[v])!=2) return 1;
  size_t NM=(size_t)1<<N; u32 *idx=malloc(NM*sizeof(u32)); if(!idx){fprintf(stderr,"memory (rank table)\n");return 1;}
  static u64 cnt[40]; memset(cnt,0,sizeof cnt);
  for(size_t mask=0;mask<NM;mask++){ int k=__builtin_popcountll(mask); idx[mask]=(u32)cnt[k]++; }
  u8 *cur=calloc(cnt[1]*(size_t)N,1), *nxt=NULL;
  for(int v=0;v<N;v++) if(st[v]) cur[(size_t)idx[(size_t)1<<v]*N+v]=1;
  for(int k=1;k<N;k++){ nxt=calloc(cnt[k+1]*(size_t)N,1); if(!nxt){fprintf(stderr,"memory (layer %d)\n",k+1);return 1;}
    // enumerate masks of popcount k in increasing order (Gosper)
    u64 mask=((u64)1<<k)-1; size_t r=0;
    while(mask<NM){ u8 *row=cur+r*N;
      for(int v=0;v<N;v++){ u8 c=row[v]; if(!c) continue; u64 nb=out[v]&~mask&(NM-1);
        while(nb){ int w=__builtin_ctzll(nb); nb&=nb-1; nxt[(size_t)idx[mask|((u64)1<<w)]*N+w]+=c; } }
      r++; u64 c0=mask&-mask, rr=mask+c0; mask=(((rr^mask)>>2)/c0)|rr; }
    free(cur); cur=nxt; }
  unsigned tot=0; printf("paths by end vertex (mod 256):");
  for(int v=0;v<N;v++){ printf(" %u", cur[v]); tot+=cur[v]; } printf("\n");
  printf("number of Hamiltonian paths starting in the start set, mod 256: %u\n", tot&255);
  return 0; }
