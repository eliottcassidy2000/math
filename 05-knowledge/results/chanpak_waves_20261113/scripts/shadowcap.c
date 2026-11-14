// Bounded-run root-change words (runs of O/Z <= RM, with the exact JZ=OJ domination pruning): count distinct
// (H, lam) pairs with H in [X/2, X] and H <= lam <= 4H, and distinct H values there. Growth exponent between scales.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X, LMAX; static unsigned char *bm; static u64 cnt; static int RM; static unsigned char *hv;
static void mark(i64 H,i64 L){ if(H<=X && H>=X/2 && L>=H && L<=4*H){ u64 idx=(u64)((H-X/2)>>1)*(u64)(3*X+2)+(u64)(L-H); if(!(bm[idx>>3]&(1<<(idx&7)))){ bm[idx>>3]|=1<<(idx&7); cnt++; } hv[(H-X/2)>>1]=1; } }
static void dfs(i64 H,i64 E,i64 L,i64 h,int run,int prevrun){
  mark(H,L); i64 D=L-E;
  if(run<RM){ if(2*L<=LMAX) dfs(H,0,2*L,h,run+1,prevrun);
    if(D<=X-H && !(D==0 && L==h && E==h) && !(run==0 && prevrun>=0 && prevrun<RM)) dfs(H+D,h,2*L-h,h,run+1,prevrun); }
  if(D<=X-H) dfs(H+D,H-E,H+L,H,0,run);
}
int main(int argc,char**argv){
  RM=atoi(argv[1]); double prev=0, prevh=0;
  for(int k=atoi(argv[2]);k<=atoi(argv[3]);k+=2){ X=1LL<<k; LMAX=4*X; cnt=0;
    u64 bits=(u64)(X/4+2)*(u64)(3*X+2); bm=calloc(bits/8+1,1); hv=calloc(X/4+2,1);
    dfs(1,1,1,1,0,-1);
    u64 nh=0; for(i64 i=0;i<X/4+1;i++) nh+=hv[i];
    printf("runs<=%d X=2^%d: distinct (H,lam) in cone %llu, distinct H %llu of %lld odd",RM,k,cnt,nh,X/4);
    if(prev>0) printf("  | pair exponent %.3f",log(cnt/prev)/log(4.0));
    printf("\n"); prev=cnt; prevh=nh; free(bm); free(hv); }
  return 0; }
