// Absolute-target solver: find two-hub words (x,y,e) with H = N by simulated annealing on log|H-N|.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128;
static int L; static int X[64],Y[64],E;
static u128 H2(void){
  u128 c[2][3][3][2]; memset(c,0,sizeof c); c[0][0][0][0]=1; c[1][0][0][0]=1;
  for(int i=0;i<L;i++){
    u128 n[2][3][3][2]; memset(n,0,sizeof n);
    for(int o=0;o<2;o++){ int h1=o?Y[i]:X[i], h2=o?X[i]:Y[i];
      for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){ u128 v=c[o][a][b][d]; if(!v) continue;
        n[o][h1?2:1][b][d]+=v;
        if(b||h1) n[o][a][h2?2:1][d]+=v;
        if(d||h2) n[o][a][b][1]+=v; } }
    memcpy(c,n,sizeof c); }
  u128 H=0;
  for(int o=0;o<2;o++){ int h12=o?!E:E;
    for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){ if(a==2||b==2||(b==0&&!h12)) continue; H+=c[o][a][b][d]; } }
  return H;
}
static double lerr(u128 H,u128 N){ u128 d=H>N?H-N:N-H; return log((double)d+1.0); }
static unsigned long long rng=88172645463325252ULL;
static unsigned long long xr(void){ rng^=rng<<13; rng^=rng>>7; rng^=rng<<17; return rng; }
static double ur(void){ return (xr()>>11)*(1.0/9007199254740992.0); }
static long long solve(u128 N,long long maxev){
  long long ev=0;
  while(ev<maxev){
    for(int i=0;i<L;i++){X[i]=xr()&1;Y[i]=xr()&1;} E=xr()&1;
    double cur=lerr(H2(),N); ev++; double T=2.0;
    for(int it=0;it<200000 && ev<maxev;it++){
      int k=xr()%(2*L+1); int *p = k<L? &X[k] : (k<2*L? &Y[k-L] : &E);
      *p^=1; u128 h=H2(); ev++; double nw=lerr(h,N);
      if(nw==0.0) return ev;
      if(nw<=cur || ur()<exp((cur-nw)/T)) cur=nw; else *p^=1;
      T*=0.9997; if(T<0.02) T=0.02;
    }
  }
  return -1;
}
int main(int argc,char**argv){
  int trials=atoi(argv[1]); double scale=atof(argv[2]); long long maxev=atoll(argv[3]);
  int ok=0; long long tot=0; rng^=(unsigned long long)(scale*1000);
  for(int t=0;t<trials;t++){
    unsigned long long n=(unsigned long long)(scale*(0.5+ur())); n|=1; u128 N=n;
    L=(int)ceil(log((double)n)/log(3.0))+3;
    long long ev=solve(N,maxev);
    if(ev>0){ok++; tot+=ev;}
    else printf("  FAILED target %llu (L=%d)\n",n,L);
  }
  printf("scale %.0e: solved %d/%d, mean evaluations %.0f\n",scale,ok,trials,ok?(double)tot/ok:0.0);
  return 0;
}
