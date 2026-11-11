// cycle_lengths.c - numbers of directed cycles of every length of a tournament (or digraph) on N <= 27 vertices,
// modulo 256.  A cycle is counted at its smallest vertex s: paths s -> ... -> v inside {s+1, ..., N-1} are counted by
// a subset dynamic programme stored layer by layer (two layers of C(m,k) * m bytes, m = N - 1 - s, and one rank
// table of 2^(N-1) 32-bit integers), and closed by an arc v -> s.
// Input (stdin): N, then N decimal out-neighbour masks.   Output: "len count_mod_256" for len = 3..N.
// With the argument "root0": only the cycles through vertex 0.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef unsigned char u8; typedef unsigned int u32;
int main(int argc,char**argv){ int N; if(scanf("%d",&N)!=1) return 1; u64 out[32]; int root0=(argc>1&&!strcmp(argv[1],"root0"));
  for(int v=0;v<N;v++) if(scanf("%llu",&out[v])!=1) return 1;
  if(N>27||N<3){ fprintf(stderr,"N out of range\n"); return 1; }
  int M=N-1; size_t NM=(size_t)1<<M; u32 *idx=malloc(NM*sizeof(u32)); if(!idx){fprintf(stderr,"memory (rank table)\n");return 1;}
  static u64 cnt[40]; memset(cnt,0,sizeof cnt);
  for(size_t mask=0;mask<NM;mask++){ int k=__builtin_popcountll(mask); idx[mask]=(u32)cnt[k]++; }
  static u64 binom[40][40]; for(int a=0;a<40;a++){ binom[a][0]=1; for(int b=1;b<=a;b++) binom[a][b]=binom[a-1][b-1]+(b<=a-1?binom[a-1][b]:0); }
  u32 c[40]; memset(c,0,sizeof c);
  for(int s=0;s<(root0?1:N-2);s++){ int m=N-1-s; u64 full=((u64)1<<m)-1; u64 o[32]; u64 start=0, back=0;
    for(int i=0;i<m;i++){ int v=s+1+i; o[i]=(out[v]>>(s+1))&full; if((out[s]>>v)&1) start|=(u64)1<<i; if((out[v]>>s)&1) back|=(u64)1<<i; }
    u8 *cur=calloc(binom[m][1]*(size_t)m,1), *nxt=NULL; if(!cur){fprintf(stderr,"memory\n");return 1;}
    for(int i=0;i<m;i++) if((start>>i)&1) cur[(size_t)idx[(size_t)1<<i]*m+i]=1;
    for(int k=1;k<=m;k++){ nxt=(k<m)?calloc(binom[m][k+1]*(size_t)m,1):NULL; if(k<m&&!nxt){fprintf(stderr,"memory (layer)\n");return 1;}
      u64 mask=((u64)1<<k)-1; size_t r=0;
      while(mask<=full){ u8 *row=cur+r*m;
        for(int v=0;v<m;v++){ u8 x=row[v]; if(!x) continue;
          if(k>=2 && ((back>>v)&1)) c[k+1]+=x;
          if(k<m){ u64 nb=o[v]&~mask; while(nb){ int w=__builtin_ctzll(nb); nb&=nb-1; nxt[(size_t)idx[mask|((u64)1<<w)]*m+w]+=x; } } }
        r++; if(k==m) break; u64 c0=mask&-mask, rr=mask+c0; mask=(((rr^mask)>>2)/c0)|rr; }
      free(cur); cur=nxt; if(!cur) break; }
    if(cur) free(cur); }
  for(int len=3;len<=N;len++) printf("%d %u\n",len,c[len]&255);
  return 0; }
