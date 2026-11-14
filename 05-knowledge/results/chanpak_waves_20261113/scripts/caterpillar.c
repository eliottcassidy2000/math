// Caterpillar tournaments (moves a=X0stay, b=X1stay, c=X1move); DFS over words, no dedupe.
// Notation: these are the repository's O/Z/J root-change words (THM-4904, THM-5148, HYP-9494): a = O, b = Z,
// c = J.  The reader blocks c a^i b^j = J O^i Z^j are not transitive; their counts are (2^i-1)(2^j-1)+2.
// State (H,E,L,h). a:(H,0,2L,h) b:(H+L-E,h,2L-h,h) c:(H+L-E,H-E,H+L,H). Start (1,1,1,1).
#include <stdio.h>
#include <stdlib.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static unsigned char *cov; static u64 nodes;
static void dfs(i64 H,i64 E,i64 L,i64 h,int depth){
  nodes++; cov[H>>1]=1;
  if(depth>200) return;
  i64 D=L-E; if(D>X-H) return;
  if(2*L<=X-H) dfs(H,0,2*L,h,depth+1);                       // a
  if(!(D==0 && L==h && E==h)) dfs(H+D,h,2*L-h,h,depth+1);     // b (skip the fixed point)
  dfs(H+D,H-E,H+L,H,depth+1);                                 // c
}
int main(int argc,char**argv){
  X=atoll(argv[1]); cov=calloc(X/2+2,1);
  dfs(1,1,1,1,0);
  u64 miss=0,last=0; int shown=0; printf("X=%lld nodes=%llu\nmisses:",X,nodes);
  for(i64 n=1;n<=X;n+=2) if(!cov[n>>1]){ miss++; last=n; if(shown<40){printf(" %lld",n);shown++;} }
  printf("\ntotal misses %llu, last miss %llu\n",miss,last);
  return 0; }
