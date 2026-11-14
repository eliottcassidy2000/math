// bounded-run census with per-dyadic-window miss counts (persistence test); same pruning as runcap2.c
#include <stdio.h>
#include <stdlib.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static u64 *bits; static int RM;
#define SETB(n) (bits[(n)>>7] |= 1ULL<<(((n)>>1)&63))
#define GETB(n) ((bits[(n)>>7]>>(((n)>>1)&63))&1ULL)
static void dfs(i64 H,i64 E,i64 L,i64 h,int run,int prevrun){
  SETB(H); i64 D=L-E; if(D>X-H) return;
  if(run<RM){ if(2*L<=X-H) dfs(H,0,2*L,h,run+1,prevrun);
    if(!(D==0 && L==h && E==h) && !(run==0 && prevrun>=0 && prevrun<RM)) dfs(H+D,h,2*L-h,h,run+1,prevrun); }
  dfs(H+D,H-E,H+L,H,0,run);
}
int main(int argc,char**argv){ X=atoll(argv[1]); RM=atoi(argv[2]); bits=calloc((X>>7)+2,sizeof(u64)); dfs(1,1,1,1,0,-1);
  printf("runs<=%d, X=%lld: misses per window [2^k,2^(k+1)):",RM,X);
  for(int k=1;(1LL<<k)<X;k++){ i64 lo=1LL<<k,hi=1LL<<(k+1); if(hi>X) hi=X; i64 m=0,last=0; for(i64 n=lo|1;n<hi;n+=2) if(!GETB(n)){m++;last=n;} printf(" %d:%lld",k,m); }
  printf("\n"); return 0; }
