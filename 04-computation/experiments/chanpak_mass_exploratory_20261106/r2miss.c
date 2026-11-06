#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
typedef long long ll;
static uint32_t *gtab; static ll X; static uint16_t *r2;
int main(int argc,char**argv){int LOGX=atoi(argv[1]);X=1LL<<LOGX;FILE*f=fopen(argv[2],"r");int nb=0;ll yh[4096],ye[4096],yR[4096];
 while(fscanf(f,"%lld %lld %lld",&yh[nb],&ye[nb],&yR[nb])==3)nb++;fclose(f);
 gtab=malloc(sizeof(uint32_t)*(size_t)(X+2));gtab[0]=1;for(ll v=1;v<=X+1;v++){if(v&1)gtab[v]=gtab[v>>1];else{ll u=v>>1;gtab[v]=(uint32_t)((ll)gtab[u]+u+(u&1));}}
 r2=calloc((size_t)(X/2+1),2);
 #pragma omp parallel for schedule(dynamic,1)
 for(int i=0;i<nb;i++)for(ll v=0;(1+v)*yR[i]<=X;v++){ll ev=(v==0||!(v&1));ll N=(ll)gtab[v]*yh[i]-ev*ye[i]+(1+v)*yR[i];if(N>=1&&N<=X&&(N&1))__atomic_fetch_add(&r2[N>>1],1,__ATOMIC_RELAXED);}
 printf("# F2 alone (bank %d): misses per dyadic block and min r2\n",nb);
 ll tot=0; printf("misses <= 4096:"); for(ll N=1;N<=4096;N+=2) if(!r2[N>>1]){printf(" %lld",N);} printf("\n");
 for(int k=4;k<LOGX;k++){ll lo=1LL<<k,hi=1LL<<(k+1),m=0,mn=1<<30;for(ll N=lo+1;N<hi;N+=2){if(N==7||N==21)continue;if(!r2[N>>1])m++;if(r2[N>>1]<mn)mn=r2[N>>1];}tot+=m;printf("k=%2d F2-misses %lld min r2 %lld\n",k,m,mn);}
 return 0;}
