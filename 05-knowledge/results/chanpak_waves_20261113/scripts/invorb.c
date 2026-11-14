// Wave 24: census of the invertible root-change language (P51, P53, P54).
// Letters on (H,E,L,h), start (1,1,1,1):  O:(H,0,2L,h)  Z:(H+L-E,h,2L-h,h)  J:(H+L-E,H-E,H+L,H).
// mode 1 (invertible): before the first J any letter; after it, Z may not follow J and J may not follow O,
//   so every block J w has w empty or w = O..Z (the invertible blocks, det 1 or 2^|w|); no reset (THM-5120).
// mode 0: the full O/Z/J grammar (for small comparisons only).
// r(N) = number of words ending in J or Z (or the empty word) with count N, saturating at 65535.
// Usage: invorb K mode   (census of odd N < 2^K)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static uint16_t *r; static u64 nodes; static int mode; static int cap=1000;
static void dfs(i64 H,i64 E,i64 L,i64 h,int last,int started,int run){
  nodes++;
  if(last!=1){ if(r[H>>1]<65535) r[H>>1]++; }
  i64 D=L-E;
  int okO=1, okZ=1, okJ=1;
  if(mode>=1 && started){ if(last==3||last==0) okZ=0; if(last==1) okJ=0; if(run>=cap){ okO=0; okZ=0; } }
  if(okO && 2*L<=X-H) dfs(H,0,2*L,h,1,started,started?run+1:0);
  if(okZ && H+D<=X && !(D==0 && L==h && E==h)) dfs(H+D,h,2*L-h,h,2,started,started?run+1:0);
  if(okJ && H+D<=X) dfs(H+D,H-E,H+L,H,3,1,0);
}
int main(int argc,char**argv){
  int K=argc>1?atoi(argv[1]):20; mode=argc>2?atoi(argv[2]):1;
  X=(1LL<<K)-1; r=calloc(X/2+2,sizeof(uint16_t));
  if(getenv("CAP")) cap=atoi(getenv("CAP"));
  dfs(1,1,1,1,0,mode==2?1:0,0);
  // sieve
  uint8_t *comp=calloc(X+1,1); for(i64 i=2;i*i<=X;i++) if(!comp[i]) for(i64 j=i*i;j<=X;j+=i) comp[j]=1;
  printf("mode=%d K=%d cap=%d nodes=%llu\n",mode,K,cap,nodes);
  u64 miss=0; i64 last=0; int shown=0; printf("misses:");
  for(i64 n=1;n<=X;n+=2) if(!r[n>>1]){ miss++; last=n; if(shown<60){printf(" %lld",n);shown++;} }
  printf("\ntotal misses %llu, largest %lld\n",miss,last);
  int ps[]={3,5,7,11,13,17,19,23,29,31};
  printf("window      odd     miss    mean   primeMean  min(at)        ratio p∤N for p=3,5,7,11,13,17,19,23,29,31\n");
  for(int k=10;k<K;k++){
    i64 a=1LL<<k, b=(1LL<<(k+1))-1; u64 n=0,m=0; double s=0,sp=0; u64 np_=0; unsigned mn=1u<<30; i64 at=0;
    double sc[10]={0}; u64 nc[10]={0};
    for(i64 N=a+1;N<=b;N+=2){ unsigned v=r[N>>1]; n++; s+=v; if(!v) m++; if(v<mn){mn=v;at=N;}
      if(!comp[N]){ sp+=v; np_++; }
      for(int t=0;t<10;t++) if(N%ps[t]){ sc[t]+=v; nc[t]++; } }
    double mean=s/n; printf("2^%-2d %9llu %8llu %8.2f %9.2f %6u(%lld) ",k,n,m,mean,np_?sp/np_:0,mn,at);
    for(int t=0;t<10;t++) printf(" %.3f",(sc[t]/nc[t])/mean); printf("\n");
  }
  if(argc>3){ FILE*f=fopen(argv[3],"wb"); fwrite(r,sizeof(uint16_t),X/2+1,f); fclose(f); }
  return 0; }
