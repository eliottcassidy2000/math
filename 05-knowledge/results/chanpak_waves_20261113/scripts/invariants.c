// Wave 25 (no Heisenberg orbit): check on every reachable state of the full O/Z/J grammar with H <= 2^K that
// 0 <= E <= H, 1 <= h <= H, L >= h, L >= E, and H(J^3 s) = 3H + 3L - E >= 3H. Also: H never decreases.
#include <stdio.h>
#include <stdlib.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static u64 nodes, bad;
static void chk(i64 H,i64 E,i64 L,i64 h){
  if(!(E>=0 && E<=H && h>=1 && h<=H && L>=h && L>=E)) bad++;
  // J^3 closed form
  i64 H1=H+L-E,E1=H-E,L1=H+L,h1=H; i64 H2=H1+L1-E1,E2=H1-E1,L2=H1+L1; i64 H3=H2+L2-E2;
  if(H3!=3*H+3*L-E || H3<3*H) bad++;
}
static void dfs(i64 H,i64 E,i64 L,i64 h,int d){
  nodes++; chk(H,E,L,h); if(d>60) return;
  i64 D=L-E; if(D<0){bad++; return;}
  if(2*L<=X-H) dfs(H,0,2*L,h,d+1);
  if(H+D<=X && !(D==0&&L==h&&E==h)) dfs(H+D,h,2*L-h,h,d+1);
  if(H+D<=X) dfs(H+D,H-E,H+L,H,d+1);
}
int main(int argc,char**argv){ int K=atoi(argv[1]); X=1LL<<K; dfs(1,1,1,1,0);
  printf("K=%d: reachable states visited %llu, invariant or J^3-bound violations %llu\n",K,nodes,bad); return 0; }
