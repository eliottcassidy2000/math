/* allags.c: q_inf(T) = P(y merges with NO translate y+r, r>=1, by Terras time T) = #groups/2^T, where words w of length T
   are grouped by (a, c(w) mod 3^a) [c(w) = 2^T T^T(y) - 3^a y]; each group has exactly one 'upward-uncertified' word
   (the one with minimal c). Also counts groups with >= 2 members etc.  Brute force over all 2^T words. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128;
static uint64_t P3[64];
static int cmp(const void *a, const void *b){ uint64_t x=*(const uint64_t*)a, y=*(const uint64_t*)b; return x<y?-1:x>y; }
int main(int argc,char**argv){
  int TM=atoi(argv[1]);
  P3[0]=1; for(int i=1;i<41;i++) P3[i]=P3[i-1]*3;
  for(int T=1;T<=TM;T++){
    /* for each weight a, collect c mod 3^a ; c up to 3^a 2^T fits in u128; store residue */
    uint64_t N=1ULL<<T; long groups=0; 
    /* bucket by weight */
    uint64_t *cnt=calloc(T+2,sizeof(uint64_t));
    for(uint64_t r=0;r<N;r++){ int a=__builtin_popcountll(r); cnt[a]++; }
    uint64_t **buf=malloc(sizeof(uint64_t*)*(T+1)); uint64_t *fill=calloc(T+1,sizeof(uint64_t));
    for(int a=0;a<=T;a++) buf[a]=malloc(sizeof(uint64_t)*(cnt[a]?cnt[a]:1));
    /* enumerate words directly: bit j of r = parity at step j (word), c(w) = sum_j 3^{a-1-idx} 2^{t_j} */
    for(uint64_t w=0; w<N; w++){
      int a=__builtin_popcountll(w); u128 c=0; int idx=0;
      for(int j=0;j<T;j++) if(w>>j&1){ c = c*3 + (((u128)1)<<j); idx++; }
      /* c computed as Horner: c = sum 3^{a-1-idx} 2^{t_idx} */
      uint64_t m = (a<=40)? (uint64_t)(c % P3[a]) : 0;
      buf[a][fill[a]++]=m;
    }
    long multi=0;
    for(int a=0;a<=T;a++){
      qsort(buf[a],fill[a],sizeof(uint64_t),cmp);
      for(uint64_t i=0;i<fill[a];i++){ if(i==0||buf[a][i]!=buf[a][i-1]) groups++; }
      free(buf[a]);
    }
    double bound=0; double binom=1;
    for(int a=0;a<=T;a++){ double c3=1; for(int i=0;i<a;i++) c3*=3; bound += (binom<c3?binom:c3); binom = binom*(T-a)/(a+1); }
    printf("%2d groups %12ld  q_inf=%.6e  sqrtT*q=%.4f  pigeonhole bound %.6e\n",T,groups,(double)groups/N,
           sqrt((double)T)*(double)groups/N, bound/N);
    free(buf); free(fill); free(cnt);
  }
}
