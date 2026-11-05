// Coverage of odd N <= X by H(u,v) = (1+u)(1+v) + g(u)g(v) - e(u)e(v), u,v >= 0.
#include <stdio.h>
#include <stdlib.h>
typedef long long ll;
int main(int argc,char**argv){
  int LOGX=argc>1?atoi(argv[1]):22; ll X=1LL<<LOGX;
  ll *g=malloc(sizeof(ll)*(X+2)); g[0]=1;
  for(ll v=1; v<=X; v++){ if(v&1) g[v]=g[v>>1]; else { ll u=v>>1; g[v]=g[u]+u+(u&1);} }
  char *hit=calloc(X+1,1); int *wu=malloc(sizeof(int)*(X+1));
  for(ll u=0; u<=X; u++){
    ll gu=g[u], eu=(u==0||!(u&1));
    for(ll v=u; (1+u)*(1+v)<=X+1; v++){   // symmetric in u,v
      ll ev=(v==0||!(v&1));
      ll H=(1+u)*(1+v)+gu*g[v]-eu*ev;
      if(H<=X && H>0 && !hit[H]){ hit[H]=1; wu[H]=(int)u; }
    }
  }
  ll miss=0,last=-1;
  printf("misses <= 5000:");
  for(ll N=1; N<=X; N+=2) if(!hit[N]){ miss++; last=N; if(N<=5000) printf(" %lld",N);}
  printf("\nodd <= 2^%d missed: %lld of %lld; largest missed %lld\n",LOGX,miss,X/2,last);
  for(int b=4;b<LOGX;b++){ ll lo=1LL<<b,hi=1LL<<(b+1),m=0,t=0; for(ll N=lo+1;N<hi;N+=2){t++; if(!hit[N]) m++;} printf("[2^%d,2^%d): miss %lld/%lld = %.5f\n",b,b+1,m,t,(double)m/t);} 
  return 0;
}
