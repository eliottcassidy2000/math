// Coverage of odd numbers by one-apex-derived families.
// g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v&1) [v>=1]; one-apex T_v: h=g(v), port x: e=[v==0 || v even], B=1+v.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
int main(int argc,char**argv){
  int LOGX = argc>1?atoi(argv[1]):20;
  ll X = 1LL<<LOGX; ll VMAX = X;
  ll *g = malloc(sizeof(ll)*(VMAX+1));
  g[0]=1;
  for(ll v=1; v<=VMAX; v++){ if(v&1) g[v]=g[v>>1]; else { ll u=v>>1; g[v]=g[u]+u+(u&1);} }
  char *c1=calloc(X+1,1),*c2=calloc(X+1,1),*c3=calloc(X+1,1),*ctow=calloc(X+1,1);
  for(ll v=0; v<=VMAX; v++){
    ll h=g[v], e=(v==0||!(v&1))?1:0, B=1+v;
    if(h<=X) c1[h]=1;
    // near-join / C2 at x: h + 2^j B, j>=1
    for(int j=1; j<62; j++){ ll N=h+(B<<j); if(N>X||N<0) break; c2[N]=1; }
    // tower at x: h1=h-e+B, then tribonacci h_k = h_{k-1}+h_{k-2}+h_{k-3} with h_0=h, h_1, h_2=h+2B
    ll a0=h, a1=h-e+B, a2=h+2*B;
    if(a1<=X) ctow[a1]=1; if(a2<=X) ctow[a2]=1;
    for(int k=3;k<200;k++){ ll a3=a2+a1+a0; if(a3>X) break; ctow[a3]=1; a0=a1;a1=a2;a2=a3; }
  }
  // report
  const char* names[4]={"one-apex h=g(v)","h+2^j B (j>=1) at x","tower at x (k>=1)","union"};
  for(int which=0; which<4; which++){
    ll miss=0, last=-1, first=-1; ll cnt=0;
    for(ll N=1; N<=X; N+=2){
      int in = which==0?c1[N]:which==1?c2[N]:which==2?ctow[N]:(c1[N]|c2[N]|ctow[N]);
      if(in) cnt++; else { miss++; last=N; if(first<0) first=N; }
    }
    printf("%-24s odd<=2^%d: hit %lld of %lld, missed %lld (largest missed %lld)\n",names[which],LOGX,cnt,X/2,miss,last);
  }
  // list missed by union below 2000
  printf("union misses <= 3000:");
  for(ll N=1; N<=3000; N+=2) if(!(c1[N]|c2[N]|ctow[N])) printf(" %lld",N);
  printf("\n");
  // density of misses in dyadic blocks for union
  for(int b=4; b<LOGX; b++){ ll lo=1LL<<b, hi=1LL<<(b+1), m=0, tot=0; for(ll N=lo+1; N<hi; N+=2){ tot++; if(!(c1[N]|c2[N]|ctow[N])) m++; } printf("[2^%d,2^%d): miss %lld/%lld = %.4f\n",b,b+1,m,tot,(double)m/tot); }
  return 0;
}
