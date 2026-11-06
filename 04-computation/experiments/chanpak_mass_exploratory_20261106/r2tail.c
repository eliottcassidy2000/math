#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
typedef long long ll;
static uint32_t *gtab; static ll X; static uint16_t *r2;
int main(int argc,char**argv){int LOGX=atoi(argv[1]);X=1LL<<LOGX;FILE*f=fopen(argv[2],"r");int nb=0;static ll yh[4096],ye[4096],yR[4096];
 while(fscanf(f,"%lld %lld %lld",&yh[nb],&ye[nb],&yR[nb])==3)nb++;fclose(f);
 gtab=malloc(sizeof(uint32_t)*(size_t)(X+2));gtab[0]=1;for(ll v=1;v<=X+1;v++){if(v&1)gtab[v]=gtab[v>>1];else{ll u=v>>1;gtab[v]=(uint32_t)((ll)gtab[u]+u+(u&1));}}
 r2=calloc((size_t)(X/2+1),2);
 #pragma omp parallel for schedule(dynamic,1)
 for(int i=0;i<nb;i++)for(ll v=0;(1+v)*yR[i]<=X;v++){ll ev=(v==0||!(v&1));ll N=(ll)gtab[v]*yh[i]-ev*ye[i]+(1+v)*yR[i];if(N>=1&&N<=X&&(N&1))__atomic_fetch_add(&r2[N>>1],1,__ATOMIC_RELAXED);}
 printf("# k: frac(r2<=t) for t=6,8,10,12,14,16,20 | mean\n");
 for(int k=16;k<LOGX;k++){ll lo=1LL<<k,hi=1LL<<(k+1);ll c[64]={0};double s=0;ll n=0;for(ll N=lo+1;N<hi;N+=2){int t=r2[N>>1];if(t<64)c[t]++;s+=t;n++;}
  ll cum=0;double fr[64];for(int t=0;t<64;t++){cum+=c[t];fr[t]=(double)cum/n;}
  printf("k=%d: %.2e %.2e %.2e %.2e %.2e %.2e %.2e | %.4f\n",k,fr[6],fr[8],fr[10],fr[12],fr[14],fr[16],fr[20],s/n);}
 return 0;}
