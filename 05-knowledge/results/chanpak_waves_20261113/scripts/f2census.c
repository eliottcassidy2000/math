// Exhaustive census of root-change (O/Z/J) words with at most CMAX root changes, bitset over odd values <= X.
// Letters a=O, b=Z, c=J; c only at the start or after b/c (JZ = OJ normal form); start (1,1,1,1).
// Prints every miss above LOW and the count below.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static u64 *bits; static int CMAX;
#define SETB(n) (bits[(n)>>7] |= 1ULL<<(((n)>>1)&63))
#define GETB(n) ((bits[(n)>>7]>>(((n)>>1)&63))&1ULL)
static void dfs(i64 H,i64 E,i64 L,i64 h,int last,int nc){
  SETB(H); i64 D=L-E; if(D>X-H) return;
  if(2*L<=X-H) dfs(H,0,2*L,h,1,nc);
  if(!(D==0 && L==h && E==h)) dfs(H+D,h,2*L-h,h,2,nc);
  if(nc<CMAX && last!=1) dfs(H+D,H-E,H+L,H,3,nc+1);
}
int main(int argc,char**argv){
  X=atoll(argv[1]); CMAX=atoi(argv[2]); i64 LOW=atoll(argv[3]);
  bits=calloc((X>>7)+2,sizeof(u64)); if(!bits){fprintf(stderr,"alloc\n");return 1;}
  dfs(1,1,1,1,0,0);
  i64 miss=0,last=0; printf("X=%lld CMAX=%d; misses above %lld:",X,CMAX,LOW);
  for(i64 n=1;n<=X;n+=2) if(!GETB(n)){ miss++; last=n; if(n>LOW) printf(" %lld",n); }
  printf("\ntotal misses %lld, largest %lld\n",miss,last); return 0; }
