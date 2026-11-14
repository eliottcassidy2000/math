#include <stdio.h>
#include <string.h>
typedef unsigned long long u64;
static u64 Hof(int*x,int*y,int len,int e){
  u64 c[2][3][3][2]; memset(c,0,sizeof c); c[0][0][0][0]=1; c[1][0][0][0]=1;
  for(int i=0;i<len;i++){ u64 n[2][3][3][2]; memset(n,0,sizeof n);
    for(int o=0;o<2;o++){ int h1=o?y[i]:x[i], h2=o?x[i]:y[i];
      for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){ u64 v=c[o][a][b][d]; if(!v) continue;
        n[o][h1?2:1][b][d]+=v; if(b||h1) n[o][a][h2?2:1][d]+=v; if(d||h2) n[o][a][b][1]+=v; } }
    memcpy(c,n,sizeof c); }
  u64 H=0; for(int o=0;o<2;o++){ int h12=o?!e:e; for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++) if(!(a==2||b==2||(b==0&&!h12))) H+=c[o][a][b][d]; }
  return H; }
int main(void){
  int X[16],Y[16],xr[16],yr[16]; long long nw=0, hit357=0; u64 vals[64]; int nv=0;
  for(int L=8;L<=10;L++) for(u64 code=0;code<(1ULL<<(2*L));code++){
    for(int i=0;i<L;i++){X[i]=(code>>(2*i))&1; Y[i]=(code>>(2*i+1))&1;}
    for(int e=0;e<2;e++){ if(Hof(X,Y,L,e)!=1297) continue; nw++;
      int got=0;
      for(int r=1;r<L;r++){ for(int i=0;i<L;i++){xr[i]=X[(i+r)%L]; yr[i]=Y[(i+r)%L];}
        for(int f=0;f<2;f++){ u64 h=Hof(xr,yr,L,f); if(h==357) got=1; int s=0; for(int q=0;q<nv;q++) if(vals[q]==h) s=1; if(!s&&nv<64) vals[nv++]=h; } }
      if(got){ hit357++; if(hit357==1){ printf("example: L=%d e=%d x=",L,e); for(int i=0;i<L;i++) printf("%d",X[i]); printf(" y="); for(int i=0;i<L;i++) printf("%d",Y[i]); printf("\n"); } } } }
  printf("words with H=1297 (length 8-10): %lld; with a cyclic spine rotation giving 357: %lld\n",nw,hit357);
  return 0; }
