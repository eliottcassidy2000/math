// collatz_procgen_20260922_mirror_mitm.c
// Memory-light reconstruction (Hirschberg-style divide and conquer) of one q-minimal closed walk through a
// negative hub h with exactly P halvings, for walks too tall for mirror_recon1 (e.g. the record loop B_1054,
// height 1.9e7).  Moves: v -> 3v+1 (always), v -> v/2 (even v); all multiplication points |v| <= HCAP.
//   forward layer values F_L(w) = least #multiplications from (L1,w1) to (L,w)   [absolute: + a1]
//   backward values     C_L(w) = least #multiplications from (L,w) to (L2,w2)
// Split at mid = (L1+L2)/2: pick w with F_mid(w) + C_mid(w) = a2, recurse on both halves; segments of length
// <= 3 are solved by a forward DP storing their layers in the same four arrays.  Memory 4 * 2 * 1.5 HCAP bytes.
// Pruning (exact):
//   MODE 0 (loops through -1): forward keeps A < (L+1)theta + log3|w| + 1 + SLACK (as mirror_loops_dp);
//          any path from -1 to (L,w) has A > L theta + log3|w|, so backward keeps C < OPT - L theta - log3|w|.
//   MODE 1 (cycles through h): forward keeps A < L theta + log3|w| - log3|h| + 1 + SLACK;
//          backward keeps C < OPT - (L theta + log3|w| - log3|h|).
// Usage: mirror_mitm H P HCAP SLACK MODE   (H = |h|; prints the word and a self-check)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#define NONE 0xFFFF
static double TH, L3; static int64_t HCAP, N, HH; static int SLACK, MODE, OPT, BASE=3;
static int64_t *minW; static const int MAXA=70000;
static uint16_t *FA,*FB,*BA,*BB; static uint16_t **SEG;
static int *TL;     // TL[L] = multiplications before the L-th halving (L = 1..P)
static void thresholds(int L){
  double base=(MODE==0)?(L+1)*TH+1.0+SLACK : L*TH-log((double)HH)/L3+1.0+SLACK;
  for(int a=0;a<MAXA;a++){ double T=a-base; if(SLACK<0||T<0){minW[a]=0;continue;}
    double w=pow(3.0,T)*(1.0-1e-9); minW[a]=(w>(double)N)?N+1:(int64_t)floor(w); }
}
static void fstep(const uint16_t *cur,uint16_t *nxt,int L){      // layer L-1 -> L
  thresholds(L); memset(nxt,255,2*N);
  for(int64_t i=1;i<N;i++){ uint16_t a=cur[i]; if(a==NONE) continue;
    int64_t x=-i; int aa=a;
    if(x&1){ if(-x>HCAP) continue; x=3*x+1; aa++; }
    for(;;){ int64_t w=-(x/2); if(aa<MAXA && w>minW[aa] && aa<nxt[w]) nxt[w]=(uint16_t)aa;
      int64_t m2=-(3*x+1); if(-x>HCAP||m2>HCAP) break; x=9*x+4; aa+=2; } }
}
static int64_t *maxI;   // backward pruning: keep C_(L-1)(u) = c only if |u| < maxI[c]
static void bstep(const uint16_t *cur,uint16_t *prv,int L){      // C_L -> C_(L-1)
  double off=(MODE==0)?0.0:log((double)HH)/L3;
  // any path to (L-1,u) has more than (L-1)theta + log3|u| - off multiplications, so on an optimal
  // path c < OPT - (L-1)theta - log3|u| + off, i.e. |u| < 3^(OPT - (L-1)theta + off - c) (permissive)
  for(int c=0;c<=OPT;c++){ double T=OPT-(L-1)*TH+off-c; double v=pow(3.0,T)*(1.0+1e-9)+1.0;
    maxI[c]=(v>(double)N)?N:(int64_t)ceil(v); }
  memset(prv,255,2*N);
  for(int64_t i=1;i<N;i++){
    int64_t x=-i; int t=0; int best=NONE;
    if(x&1){ if(-x>HCAP) continue; x=3*x+1; t=1; }
    for(;;){ int64_t w=-(x/2); if(w<N && cur[w]!=NONE && t+cur[w]<best) best=t+cur[w];
      int64_t m2=-(3*x+1); if(-x>HCAP||m2>HCAP) break; x=9*x+4; t+=2; }
    if(best==NONE || best>OPT) continue;
    if(i<maxI[best]) prv[i]=(uint16_t)best;
  }
}
static void solve(int L1,int64_t w1,int a1,int L2,int64_t w2,int a2){
  if(L2-L1<=BASE){
    memset(SEG[0],255,2*N); SEG[0][-w1]=(uint16_t)a1;
    for(int L=L1+1;L<=L2;L++) fstep(SEG[L-1-L1],SEG[L-L1],L);
    if(SEG[L2-L1][-w2]!=a2){ fprintf(stderr,"base mismatch L=%d..%d got %d want %d\n",L1,L2,SEG[L2-L1][-w2],a2); exit(1); }
    int64_t w=w2; int A=a2;
    for(int L=L2;L>L1;L--){ int64_t x=2*w; int found=0;
      for(int t=0;t<=80;t++){ int64_t u=x; int good=1; int64_t mx=0;
        for(int s=0;s<t;s++){ int64_t y=u-1; if(y%3!=0){good=0;break;} u=y/3; if(-u>mx) mx=-u; }
        if(!good) break; if(t>0&&mx>HCAP) continue; if(((u&1)!=0)!=((t&1)!=0)) continue; if(-u>=N||-u<1) continue;
        uint16_t au=SEG[L-1-L1][-u]; if(au!=NONE && au+t==A){ TL[L]=t; A=au; w=u; found=1; break; } }
      if(!found){ fprintf(stderr,"backtrack failed at L=%d\n",L); exit(1); } }
    if(w!=w1||A!=a1){ fprintf(stderr,"segment end mismatch\n"); exit(1); }
    return;
  }
  int mid=(L1+L2)/2;
  memset(FA,255,2*N); FA[-w1]=(uint16_t)a1; uint16_t *c=FA,*n=FB;
  for(int L=L1+1;L<=mid;L++){ fstep(c,n,L); uint16_t*t=c;c=n;n=t; }
  memset(BA,255,2*N); BA[-w2]=0; uint16_t *bc=BA,*bn=BB;
  for(int L=L2;L>mid;L--){ bstep(bc,bn,L); uint16_t*t=bc;bc=bn;bn=t; }
  int64_t wm=0; int am=-1;
  for(int64_t i=1;i<N;i++) if(c[i]!=NONE && bc[i]!=NONE && c[i]+bc[i]==a2){ wm=-i; am=c[i]; break; }
  if(am<0){ fprintf(stderr,"no midpoint at L=%d (segment %d..%d)\n",mid,L1,L2); exit(1); }
  fprintf(stderr,"  split %d..%d at %d: state %lld, a=%d\n",L1,L2,mid,(long long)wm,am);
  solve(L1,w1,a1,mid,wm,am); solve(mid,wm,am,L2,w2,a2);
}
int main(int argc,char**argv){
  if(argc<6){fprintf(stderr,"usage: %s H P HCAP SLACK MODE [BASE]\n",argv[0]);return 1;}
  HH=atoll(argv[1]); int P=atoi(argv[2]); HCAP=atoll(argv[3]); SLACK=atoi(argv[4]); MODE=atoi(argv[5]);
  TH=log(2.0)/log(3.0); L3=log(3.0); N=(3*HCAP)/2+2; minW=malloc(sizeof(int64_t)*MAXA);
  FA=malloc(2*N); FB=malloc(2*N); BA=malloc(2*N); BB=malloc(2*N);
  if(!FA||!FB||!BA||!BB){fprintf(stderr,"oom\n");return 1;}
  SEG=malloc(sizeof(uint16_t*)*4); SEG[0]=FA; SEG[1]=FB; SEG[2]=BA; SEG[3]=BB;   // reused by the base case
  TL=calloc(P+2,sizeof(int)); maxI=malloc(sizeof(int64_t)*(MAXA+1));
  // optimum: forward pass over all P layers
  memset(FA,255,2*N); FA[HH]=0; uint16_t *c=FA,*n=FB;
  for(int L=1;L<=P;L++){ fstep(c,n,L); uint16_t*t=c;c=n;n=t; }
  if(c[HH]==NONE){ printf("MITM h=-%lld P=%d: none within cap/slack\n",(long long)HH,P); return 0; }
  OPT=c[HH]; fprintf(stderr,"optimum q=%d\n",OPT);
  solve(0,-HH,0,P,-HH,OPT);
  // print and self-check
  int64_t v=-HH,hmax=0; int aa=0,ok=1; char last=0; int run=0;
  printf("RECON h=-%lld P=%d q=%d word:",(long long)HH,P,OPT);
  for(int L=1;L<=P;L++){
    for(int s=0;s<TL[L];s++){ if(-v>hmax) hmax=-v; v=3*v+1; aa++; if(last!='M'){ if(run) printf(" %c%d",last,run); last='M'; run=0;} run++; }
    if(v&1){ok=0;break;} v/=2; if(last!='H'){ if(run) printf(" %c%d",last,run); last='H'; run=0;} run++;
  }
  if(run) printf(" %c%d",last,run);
  printf("\nRECON_CHECK h=-%lld P=%d q=%d H=%lld end=%lld %s\n",(long long)HH,P,aa,(long long)hmax,(long long)v,(ok&&v==-HH&&aa==OPT)?"OK":"BAD");
  return 0;
}
