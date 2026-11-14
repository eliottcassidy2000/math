// Exhaustive values of root-change (O/Z/J) words; letters a=O, b=Z, c=J; with at most CMAX root changes (c only after b/c or at start; JZ=OJ normal form
// makes 'c after a' redundant), from the single-vertex start (1,1,1,1). Prints misses above a threshold.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
static i64 X; static unsigned char *cov; static int CMAX;
static void dfs(i64 H,i64 E,i64 L,i64 h,int last,int nc){
  cov[H>>1]=1;
  i64 D=L-E; if(D>X-H) return;
  if(2*L<=X-H) dfs(H,0,2*L,h,1,nc);
  if(!(D==0 && L==h && E==h)) dfs(H+D,h,2*L-h,h,2,nc);
  if(nc<CMAX && last!=1) dfs(H+D,H-E,H+L,H,3,nc+1);
}
int main(int argc,char**argv){
  X=atoll(argv[1]); CMAX=atoi(argv[2]); cov=calloc(X/2+2,1);
  dfs(1,1,1,1,0,0);
  i64 miss=0,last=0; int shown=0; printf("X=%lld CMAX=%d\nmisses:",X,CMAX);
  for(i64 n=1;n<=X;n+=2) if(!cov[n>>1]){ miss++; last=n; if(shown<60){printf(" %lld",n);shown++;} }
  printf("\ntotal %lld, last %lld\n",miss,last);
  // misses in the top half
  i64 top=0; for(i64 n=X/2|1;n<=X;n+=2) if(!cov[n>>1]) top++; printf("misses in (X/2, X]: %lld\n",top);
  return 0; }
