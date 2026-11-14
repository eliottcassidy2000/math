// Misses of root-change (O/Z/J) words; letters a=O, b=Z, c=J; with at most CMAX root changes, reported per dyadic window (persistence test).
#include <stdio.h>
#include <stdlib.h>
typedef long long i64;
static i64 X; static unsigned char *cov; static int CMAX;
static void dfs(i64 H,i64 E,i64 L,i64 h,int last,int nc){
  cov[H>>1]=1; i64 D=L-E; if(D>X-H) return;
  if(2*L<=X-H) dfs(H,0,2*L,h,1,nc);
  if(!(D==0 && L==h && E==h)) dfs(H+D,h,2*L-h,h,2,nc);
  if(nc<CMAX && last!=1) dfs(H+D,H-E,H+L,H,3,nc+1);
}
int main(int argc,char**argv){
  X=atoll(argv[1]); CMAX=atoi(argv[2]); cov=calloc(X/2+2,1); dfs(1,1,1,1,0,0);
  printf("CMAX=%d X=%lld: misses per dyadic window [2^k,2^(k+1)):\n",CMAX,X);
  for(int k=1;(1LL<<k)<X;k++){ i64 lo=1LL<<k, hi=(1LL<<(k+1)); if(hi>X) hi=X; i64 m=0,last=0;
    for(i64 n=lo|1;n<hi;n+=2) if(!cov[n>>1]){m++;last=n;}
    printf("  k=%2d misses %8lld  (odd numbers %lld, fraction %.2e) last %lld\n",k,m,(hi-lo)/2,(double)m/((hi-lo)/2),last); }
  return 0; }
