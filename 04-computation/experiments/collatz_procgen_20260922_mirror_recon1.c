// collatz_procgen_20260922_mirror_recon1.c
// Q1 mirror of loops_recon1.c + loops_hubs.c: closed walks through a negative hub h in the E-graph
// (moves v -> 3v+1 always, v -> v/2 for even v), layered by halvings, minimising the number q of
// multiplications (then the height, max |multiplication point|, if heights are on).
//   * h = -1, MODE 0: loops through -1 (pruning uses the loop bound a >= ceil((K+1) log_3 2)).
//   * any h, MODE 1: cycles through h (pruning uses the cycle bound 3^q |h| >= 2^p |h| + B, i.e.
//     q >= p log_3 2; a cycle through h with (p,q) has ratio 3^q/2^p = prod 3|v|/(3|v|-1) > 1).
// Pruning keeps a state w at layer L only if A_L(w) < base_L + log_3|w| where
//   MODE 0: base_L = (L+1) log_3 2 + 1 + SLACK;  MODE 1: base_L = L log_3 2 - log_3|h| + 1 + SLACK;
// this is exact for walks with q <= ceil(bound) + SLACK (SLACK < 0: no pruning).
// Usage:
//   mirror_recon1 H P HCAP SLACK MODE scan       : for p = 1..P print min q at h (hub table)
//   mirror_recon1 H P HCAP SLACK MODE recon [CK] : reconstruct one q-minimal closed walk with exactly
//        P halvings (checkpoints every CK layers, default ~sqrt(P)); memory ~ (P/CK + CK) * 2 * 1.5 HCAP bytes
//   (H is |h|, a positive integer.)  Words are printed as forward runs "M3 H1 M1 ...".
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#define NONE 0xFFFF
static double TH; static int64_t HCAP, N, HH; static int SLACK, MODE;
static int64_t *minW; static int MAXA=70000;
static void set_thresholds(int L){
  double base = (MODE==0)? (L+1)*TH+1.0+SLACK : L*TH - log((double)HH)/log(3.0) + 1.0 + SLACK;
  for(int a=0;a<MAXA;a++){
    if(SLACK<0){ minW[a]=0; continue; }
    double T=a-base; if(T<0){ minW[a]=0; continue; }
    double w=pow(3.0,T)*(1.0-1e-9);
    minW[a]=(w>(double)N)? N+1 : (int64_t)floor(w);
  }
}
// one layer step: from cur (layer L-1) to nxt (layer L); returns live count
static long long step(const uint16_t *cur, uint16_t *nxt, int L){
  set_thresholds(L);
  memset(nxt,255,2*N); long long live=0;
  for(int64_t i=1;i<N;i++){
    uint16_t a=cur[i]; if(a==NONE) continue; live++;
    int64_t x=-i; int aa=a;
    if(x&1){ if(-x>HCAP) continue; x=3*x+1; aa++; }
    for(;;){
      int64_t w=-(x/2);
      if(aa<MAXA && w>minW[aa] && aa<nxt[w]) nxt[w]=(uint16_t)aa;
      int64_t m2=-(3*x+1); if(-x>HCAP || m2>HCAP) break;
      x=9*x+4; aa+=2;
    }
  }
  return live;
}
static int ceil_ptheta(int p){ return (int)ceil(p*TH); }
int main(int argc,char**argv){
  if(argc<7){fprintf(stderr,"usage: %s H P HCAP SLACK MODE scan|recon [CK]\n",argv[0]);return 1;}
  HH=atoll(argv[1]); int P=atoi(argv[2]); HCAP=atoll(argv[3]); SLACK=atoi(argv[4]); MODE=atoi(argv[5]);
  int recon=(strcmp(argv[6],"recon")==0); int CK=argc>7?atoi(argv[7]):(int)ceil(sqrt((double)P));
  TH=log(2.0)/log(3.0); N=(3*HCAP)/2+2; if(HH>=N){fprintf(stderr,"hub beyond cap\n");return 1;}
  minW=malloc(sizeof(int64_t)*MAXA);
  if(!recon){
    uint16_t *c=malloc(2*N),*n=malloc(2*N); memset(c,255,2*N); c[HH]=0;
    for(int L=1;L<=P;L++){
      long long live=step(c,n,L);
      if(n[HH]!=NONE){
        int q=n[HH]; double r=exp2(q*log2(3.0)-L); int cb=(MODE==0)?(int)ceil((L+1)*TH):ceil_ptheta(L);
        printf("h=-%lld p=%4d q=%4d ratio=%.6f bound=%d %s live=%lld\n",(long long)HH,L,q,r,cb,q==cb?"=bound":">bound",live);
      }
      fflush(stdout); uint16_t*t=c;c=n;n=t;
    }
    return 0;
  }
  // recon: checkpoints at layers 0, CK, 2CK, ...
  int nck=P/CK+1;
  uint16_t **ck=malloc(sizeof(uint16_t*)*(nck+1));
  uint16_t *c=malloc(2*N),*n=malloc(2*N); memset(c,255,2*N); c[HH]=0;
  ck[0]=malloc(2*N); memcpy(ck[0],c,2*N);
  for(int L=1;L<=P;L++){
    step(c,n,L); uint16_t*t=c;c=n;n=t;
    if(L%CK==0 && L/CK<=nck){ ck[L/CK]=malloc(2*N); if(!ck[L/CK]){fprintf(stderr,"oom ck\n");return 1;} memcpy(ck[L/CK],c,2*N); }
  }
  if(c[HH]==NONE){ printf("RECON h=-%lld P=%d none within cap/slack\n",(long long)HH,P); return 0; }
  int Aend=c[HH];
  // backtrack: segment by segment
  int *tl=malloc(sizeof(int)*(P+1));
  uint16_t **seg=malloc(sizeof(uint16_t*)*(CK+1)); for(int i=0;i<=CK;i++){ seg[i]=malloc(2*N); if(!seg[i]){fprintf(stderr,"oom seg\n");return 1;} }
  int64_t w=-HH; int A=Aend; int L=P;
  while(L>0){
    int s0=((L-1)/CK)*CK;              // segment start layer (checkpoint)
    memcpy(seg[0],ck[s0/CK],2*N);
    for(int j=s0+1;j<=L;j++) step(seg[j-s0-1],seg[j-s0],j);
    for(;L>s0;L--){
      int64_t x=2*w; int found=0;
      for(int t=0;t<=80;t++){
        int64_t u=x; int good=1; int64_t mx=0;
        for(int s=0;s<t;s++){ int64_t y=u-1; if(y%3!=0){good=0;break;} u=y/3; if(-u>mx) mx=-u; }
        if(!good) break;
        if(t>0 && mx>HCAP) continue;
        if(((u&1)!=0)!=((t&1)!=0)) continue;
        if(-u>=N || -u<1) continue;
        uint16_t au=seg[L-1-s0][-u];
        if(au!=NONE && au+t==A){ tl[L]=t; A=au; w=u; found=1; break; }
      }
      if(!found){ printf("RECON FAILED at layer %d\n",L); return 1; }
    }
  }
  if(w!=-HH || A!=0){ printf("RECON FAILED (end %lld, A %d)\n",(long long)w,A); return 1; }
  // forward word
  int64_t v=-HH, hmax=0; int aa=0; char last=0; int run=0; int ok=1;
  printf("RECON h=-%lld P=%d q=%d word:",(long long)HH,P,Aend);
  for(int L2=1;L2<=P;L2++){
    for(int s=0;s<tl[L2];s++){ if(-v>hmax) hmax=-v; v=3*v+1; aa++; if(last!='M'){ if(run) printf(" %c%d",last,run); last='M'; run=0;} run++; }
    if(v&1){ok=0;break;} v/=2; if(last!='H'){ if(run) printf(" %c%d",last,run); last='H'; run=0;} run++;
  }
  if(run) printf(" %c%d",last,run);
  printf("\nRECON_CHECK h=-%lld P=%d q=%d H=%lld end=%lld %s\n",(long long)HH,P,aa,(long long)hmax,(long long)v,(ok&&v==-HH&&aa==Aend)?"OK":"BAD");
  return 0;
}
