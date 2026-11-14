// Full root-change grammar (letters a=O, b=Z, c=J; c only at start or after b/c): which (H, lam) pairs are reached,
// with H <= X and lam <= R*X? Reports the fraction of lattice points (H odd) in that box that are reached, per scale.
#include <stdio.h>
#include <stdlib.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X, LMAX; static unsigned char *bm; static u64 cnt;
static void mark(i64 H,i64 L){ if(H<=X && L<=LMAX && L>=1){ u64 idx=(u64)(H>>1)*(u64)(LMAX+1)+(u64)L; if(!(bm[idx>>3]&(1<<(idx&7)))){ bm[idx>>3]|=1<<(idx&7); cnt++; } } }
static void dfs(i64 H,i64 E,i64 L,i64 h,int last){
  mark(H,L); i64 D=L-E; if(D>X-H) { if(2*L<=LMAX) dfs(H,0,2*L,h,1); return; }
  if(2*L<=LMAX || 2*L<=X-H) dfs(H,0,2*L,h,1);
  if(!(D==0 && L==h && E==h)) dfs(H+D,h,2*L-h,h,2);
  if(last!=1) dfs(H+D,H-E,H+L,H,3);
}
int main(int argc,char**argv){
  X=atoll(argv[1]); i64 R=atoll(argv[2]); LMAX=R*X;
  u64 bits=(u64)(X/2+1)*(u64)(LMAX+1); bm=calloc(bits/8+1,1);
  dfs(1,1,1,1,0);
  // fraction of reached lattice points among (H odd <= X, 1 <= L <= LMAX) and in the sub-box H in [X/2, X], L in [H, R*H]
  u64 tot=0, hit=0;
  for(i64 H=(X/2)|1; H<=X; H+=2) for(i64 L=H; L<=R*H && L<=LMAX; L++){ tot++; u64 idx=(u64)(H>>1)*(u64)(LMAX+1)+(u64)L; if(bm[idx>>3]&(1<<(idx&7))) hit++; }
  printf("X=%lld: reached (H,lam) pairs %llu; in the cone H in [X/2,X], H <= lam <= %lld H: %llu of %llu lattice points (%.5f)\n",X,cnt,R,hit,tot,(double)hit/tot);
  return 0; }
