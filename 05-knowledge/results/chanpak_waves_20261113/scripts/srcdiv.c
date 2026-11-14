// Wave 24, P58: in the invertible language (free one-apex source, then invertible blocks), record for each word
// (ending in J or Z, with at least one J) whether its source state (H,E,L) at the first J is p-divisible, p in
// {3,5,7,11,13}. Per dyadic window of the count: pi_p = share of such words; measured ratio
// (mean r over N prime to p)/(mean r) against the prediction (1 - pi_p) p^3/(p^3 - 1).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static int K; static int P[5]={3,5,7,11,13};
static u64 tot[64], divc[64][5], cop[64][5], nocj[64];
static void rec(i64 H,int mask,int started){ int k=63-__builtin_clzll((u64)H); if(!started){ nocj[k]++; return; }
  tot[k]++; for(int t=0;t<5;t++){ if(mask>>t&1) divc[k][t]++; if(H%P[t]) cop[k][t]++; } }
static void dfs(i64 H,i64 E,i64 L,i64 h,int last,int started,int mask){
  if(last!=1) rec(H,mask,started);
  i64 D=L-E; int okO=1,okZ=1,okJ=1;
  if(started){ if(last==3) okZ=0; if(last==1) okJ=0; }
  if(okO && 2*L<=X-H) dfs(H,0,2*L,h,1,started,mask);
  if(okZ && H+D<=X && !(D==0 && L==h && E==h)) dfs(H+D,h,2*L-h,h,2,started,mask);
  if(okJ && H+D<=X){ int m=mask; if(!started){ m=0; for(int t=0;t<5;t++) if(H%P[t]==0 && E%P[t]==0 && L%P[t]==0) m|=1<<t; }
    dfs(H+D,H-E,H+L,H,3,1,m); }
}
int main(int argc,char**argv){ K=argc>1?atoi(argv[1]):22; X=(1LL<<K)-1; dfs(1,1,1,1,0,0,0);
  printf("window  words(with J)  noJ   |  per p: pi_p  measured-share-of-words-with-p|H-coprime  predicted ratio (1-pi)p^3/(p^3-1)\n");
  for(int k=12;k<K;k++){ printf("2^%-2d %12llu %6llu |",k,tot[k],nocj[k]);
    for(int t=0;t<5;t++){ double pi=(double)divc[k][t]/tot[k]; double p=P[t]; double pred=(1-pi)*p*p*p/(p*p*p-1);
      double meas=((double)cop[k][t]/tot[k])/(1-1/p); printf("  p=%d pi=%.4f pred=%.4f wordratio=%.4f",P[t],pi,pred,meas); }
    printf("\n"); }
  return 0; }
