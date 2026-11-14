// Count bounded-run root-change words with value <= X, excluding words ending in O (they repeat their prefix's value),
// with the exact JZ = OJ domination pruning: an estimate of the orbit-growth exponent delta_L.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static int RM; static u64 cnt;
static void dfs(i64 H,i64 E,i64 L,i64 h,int run,int prevrun,int lastO){
  if(!lastO) cnt++; i64 D=L-E; if(D>X-H) return;
  if(run<RM){ if(2*L<=X-H) dfs(H,0,2*L,h,run+1,prevrun,1);
    if(!(D==0 && L==h && E==h) && !(run==0 && prevrun>=0 && prevrun<RM)) dfs(H+D,h,2*L-h,h,run+1,prevrun,0); }
  dfs(H+D,H-E,H+L,H,0,run,0);
}
int main(int argc,char**argv){ RM=atoi(argv[1]); double prev=0;
  for(int k=14;k<=22;k+=2){ X=1LL<<k; cnt=0; dfs(1,1,1,1,0,-1,0);
    printf("runs<=%d X=2^%d: words %llu  (per odd number %.3f)%s",RM,k,cnt,(double)cnt/(X/2), prev>0? "":"\n");
    if(prev>0) printf("  local exponent %.3f\n", log(cnt/prev)/log(4.0)); prev=cnt; }
  return 0; }
