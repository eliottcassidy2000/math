// As runcap.c, with the exact JZ = OJ domination: a Z as the first letter after J is skipped when the run before
// that J had length < RM (then "O J" reaches the same root state with a shorter current run, a superset of futures).
#include <stdio.h>
#include <stdlib.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static u64 *bits; static int RM; static u64 nodes;
#define SETB(n) (bits[(n)>>7] |= 1ULL<<(((n)>>1)&63))
#define GETB(n) ((bits[(n)>>7]>>(((n)>>1)&63))&1ULL)
static void dfs(i64 H,i64 E,i64 L,i64 h,int run,int prevrun){
  nodes++; SETB(H); i64 D=L-E; if(D>X-H) return;
  if(run<RM){
    if(2*L<=X-H) dfs(H,0,2*L,h,run+1,prevrun);
    if(!(D==0 && L==h && E==h) && !(run==0 && prevrun>=0)) dfs(H+D,h,2*L-h,h,run+1,prevrun);
  }
  dfs(H+D,H-E,H+L,H,0,run);
}
int main(int argc,char**argv){
  X=atoll(argv[1]); RM=atoi(argv[2]); i64 LOW=atoll(argv[3]);
  bits=calloc((X>>7)+2,sizeof(u64)); dfs(1,1,1,1,0,-1);
  i64 miss=0,last=0; printf("X=%lld runs<=%d nodes=%llu; misses above %lld:",X,RM,nodes,LOW);
  for(i64 n=1;n<=X;n+=2) if(!GETB(n)){ miss++; last=n; if(n>LOW) printf(" %lld",n); }
  printf("\ntotal misses %lld, largest %lld\n",miss,last); return 0; }
