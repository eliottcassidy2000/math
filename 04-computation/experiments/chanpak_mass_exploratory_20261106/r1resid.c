#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
typedef long long ll;
static uint32_t *gtab; static ll X; static uint16_t *r1;
int main(int argc,char**argv){int LOGX=atoi(argv[1]);X=1LL<<LOGX;
 gtab=malloc(sizeof(uint32_t)*(size_t)(X+2));gtab[0]=1;for(ll v=1;v<=X+1;v++){if(v&1)gtab[v]=gtab[v>>1];else{ll u=v>>1;gtab[v]=(uint32_t)((ll)gtab[u]+u+(u&1));}}
 r1=calloc((size_t)(X/2+1),2);ll umax=0;while((umax+2)*(umax+2)<=X+1)umax++;
 #pragma omp parallel for schedule(dynamic,1)
 for(ll u=0;u<=umax;u++){ll gu=gtab[u],eu=(u==0||!(u&1));for(ll v=u;(1+u)*(1+v)<=X;v++){ll ev=(v==0||!(v&1));ll N=(1+u)*(1+v)+gu*(ll)gtab[v]-eu*ev;if(N>=1&&N<=X&&(N&1))__atomic_fetch_add(&r1[N>>1],1,__ATOMIC_RELAXED);}}
 int K=LOGX-1;ll lo=1LL<<K,hi=1LL<<(K+1);
 int mods[]={3,9,27,81};
 for(int mi=0;mi<4;mi++){int m=mods[mi];double s[100]={0};ll c[100]={0},z[100]={0};
  for(ll N=lo+1;N<hi;N+=2){int r=N%m;s[r]+=r1[N>>1];c[r]++;if(!r1[N>>1])z[r]++;}
  double mn=1e9,mx=0;int rmn=-1;for(int r=0;r<m;r++)if(c[r]){double a=s[r]/c[r];if(a<mn){mn=a;rmn=r;}if(a>mx)mx=a;}
  printf("mod %d: min class mean %.3f (r=%d, F1-miss frac %.5f), max class mean %.3f\n",m,mn,rmn,(double)z[rmn]/c[rmn],mx);
  if(m<=9){printf("   class means:");for(int r=0;r<m;r++)if(c[r])printf(" %d:%.3f/%.5f",r,s[r]/c[r],(double)z[r]/c[r]);printf("\n");}
 }
 return 0;}
