// collatz_procgen_20260922_loops_seesaw.c
// Backward-E (Q2) descent game with exact certificate bookkeeping, level r (classes mod 3^(r+1)).
// For a 3-adic unit class c, a certificate of length <= r is a legal reverse path whose legality is decided
// by c mod 3^(r+1); its multiplier after a prefix of s moves with K halvings is 2^K/3^s.
// best_r(c) = min over legal paths of length <= r and over their prefixes (s >= 1) of 2^K/3^s,
// stored exactly as the pair (s,K).  Exceptional: best_r(c) >= 1.
// Reports: the exceptional list; the worst certified factor max{best_r(c) : best_r(c) < 1};
// the same worst factor restricted to classes at 3-adic distance >= 3^-j from every exceptional class
// (i.e. c mod 3^j differs from e mod 3^j for all exceptional e), for j = 1..r+1;
// the fraction of classes with best factor >= theta for several theta; and the certified factors of the
// classes on the hostile threads 1 and 1/2 (c = 1 + 3^k u, c = (1 + 3^k w)/2) as functions of k.
// Usage: loops_seesaw R [dumpfile]   (R <= 14; memory ~ 3^(R+1) * 4 bytes * 2); dumpfile lists all classes
//        with best factor >= 1/2 as 'class s K d cert|EXC'.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
typedef struct { int16_t s, K; } SK;      // s = 0 means "no certificate / empty"
static long double val(SK a){ return a.s? a.K*logl(2.0L)-a.s*logl(3.0L) : INFINITY; }
int main(int argc,char**argv){
  int R=atoi(argv[1]);
  uint64_t P[40]; P[0]=1; for(int i=1;i<40;i++) P[i]=P[i-1]*3;
  SK *prev=NULL,*cur=NULL;
  for(int r=1;r<=R;r++){
    uint64_t M=P[r+1];
    cur=malloc(sizeof(SK)*M);
    int kmax=(int)(r*log(3.0)/log(2.0))+2;
    for(uint64_t c=0;c<M;c++){
      SK best={0,0}; long double bv=INFINITY;
      if(c%3==0){ cur[c]=best; continue; }
      int k0=(c%3==1)?0:1; uint64_t pw=(k0==0)?1:2;
      for(int k=k0;k<=kmax;k+=2){
        uint64_t t=(pw*c)%M; uint64_t y=((t+M-1)%M)/3;   // y mod 3^r
        if(y%3!=0){
          // one move: (1,k); then optionally continue with the best certificate of y at level r-1 (if < 1)
          SK one={1,(int16_t)k}; long double v1=val(one);
          if(v1<bv){ bv=v1; best=one; }
          if(r>1){ SK nx=prev[y%P[r]]; if(nx.s){ long double v2=v1+val(nx); if(val(nx)<0 && v2<bv){ bv=v2; best.s=1+nx.s; best.K=k+nx.K; } } }
        }
        pw=(pw*4)%M;
      }
      cur[c]=best;
    }
    if(prev) free(prev); prev=cur;
  }
  uint64_t M=P[R+1];
  // exceptional classes
  uint64_t nexc=0; uint64_t *exc=malloc(sizeof(uint64_t)*M);
  for(uint64_t c=0;c<M;c++) if(c%3 && !(val(cur[c])<0)) exc[nexc++]=c;
  printf("level r=%d (mod 3^%d): exceptional classes %llu\n",R,R+1,(unsigned long long)nexc);
  printf("  exceptional:"); for(uint64_t i=0;i<nexc;i++) printf(" %llu",(unsigned long long)exc[i]); printf("\n");
  // distance: d(c) = max j<=R+1 with c == e mod 3^j for some exceptional e
  // compute via marking residues of exceptional classes mod 3^j
  uint8_t *dist=calloc(M,1);
  for(int j=1;j<=R+1;j++){
    uint8_t *mark=calloc(P[j],1);
    for(uint64_t i=0;i<nexc;i++) mark[exc[i]%P[j]]=1;
    for(uint64_t c=0;c<M;c++) if(mark[c%P[j]]) dist[c]=j;
    free(mark);
  }
  // worst factor overall and by distance
  long double worst=-INFINITY; SK wsk={0,0}; uint64_t wc=0;
  long double wd[40]; SK wdk[40]; uint64_t wdc[40]; uint64_t cnt[40];
  for(int j=0;j<=R+1;j++){ wd[j]=-INFINITY; cnt[j]=0; }
  double th[]={0.5,0.6,0.7,0.8,0.9,0.95,0.98}; uint64_t nth[7]={0}; uint64_t nunit=0,nok=0;
  for(uint64_t c=0;c<M;c++){
    if(c%3==0) continue; nunit++;
    long double v=val(cur[c]); if(!(v<0)) continue; nok++;
    long double f=expl(v);
    for(int t=0;t<7;t++) if(f>=th[t]) nth[t]++;
    if(v>worst){ worst=v; wsk=cur[c]; wc=c; }
    int d=dist[c]; cnt[d]++;
    if(v>wd[d]){ wd[d]=v; wdk[d]=cur[c]; wdc[d]=c; }
  }
  printf("  units %llu, certified %llu\n",(unsigned long long)nunit,(unsigned long long)nok);
  printf("  worst certified factor overall: 2^%d/3^%d = %.6Lf at class %llu (dist %d)\n",wsk.K,wsk.s,expl(worst),(unsigned long long)wc,dist[wc]);
  for(int t=0;t<7;t++) printf("  fraction of certified classes with factor >= %.2f: %.3e (%llu)\n",th[t],(double)nth[t]/nok,(unsigned long long)nth[t]);
  printf("  worst certified factor among classes whose 3-adic distance to the exceptional set is exactly 3^-d\n");
  printf("  (d = largest j with c = e mod 3^j):\n");
  for(int j=0;j<=R+1;j++) if(cnt[j]) printf("   d=%2d  classes %9llu  worst 2^%d/3^%d = %.6Lf (class %llu)\n",j,(unsigned long long)cnt[j],wdk[j].K,wdk[j].s,expl(wd[j]),(unsigned long long)wdc[j]);
  // cumulative: max over d <= j
  printf("  cumulative worst over classes with d <= j:\n");
  long double cm=-INFINITY; SK cmk={0,0};
  for(int j=0;j<=R+1;j++){ if(cnt[j] && wd[j]>cm){ cm=wd[j]; cmk=wdk[j]; } if(cnt[j]) printf("   d<=%2d: 2^%d/3^%d = %.6Lf\n",j,cmk.K,cmk.s,expl(cm)); }
  // best descent rate per consumed digit: for each class c and each length s <= R let Kmin_s(c) be the least
  // number of halvings over legal paths of exactly s moves (determined by c mod 3^(s+1)); the class's best
  // rate is max_s (s ln3 - Kmin_s ln2)/s (nats per move = per consumed 3-adic digit).
  {
    uint16_t **T=malloc(sizeof(uint16_t*)*(R+1));
    T[0]=malloc(2*3); T[0][0]=T[0][1]=T[0][2]=0;
    for(int s=1;s<=R;s++){
      uint64_t Ms=P[s+1]; T[s]=malloc(2*Ms);
      for(uint64_t c=0;c<Ms;c++){
        uint16_t best=0xFFFF;
        if(c%3){ int k0=(c%3==1)?0:1; uint64_t pw=(k0==0)?1:2;
          for(int k=k0;k<=2*s+8;k+=2){ uint64_t t=(pw*c)%Ms; uint64_t y=((t+Ms-1)%Ms)/3;
            if(y%3!=0){ uint16_t prevK=T[s-1][y%P[s]]; if(prevK!=0xFFFF && k+prevK<best) best=k+prevK; }
            pw=(pw*4)%Ms; } }
        T[s][c]=best;
      }
    }
    long double mr[40]; uint64_t mc[40]; int ms[40];
    for(int j=0;j<=R+1;j++) mr[j]=INFINITY;
    for(uint64_t c=0;c<M;c++){ if(c%3==0) continue; if(!(val(cur[c])<0)) continue;
      long double br=-INFINITY; int bs=0;
      for(int s=1;s<=R;s++){ uint16_t K=T[s][c%P[s+1]]; if(K==0xFFFF) continue;
        long double r=(s*logl(3.0L)-K*logl(2.0L))/s; if(r>br){ br=r; bs=s; } }
      int d=dist[c]; if(br<mr[d]){ mr[d]=br; mc[d]=c; ms[d]=bs; } }
    printf("  best descent rate per consumed digit, max_s (s ln3 - Kmin_s ln2)/s, worst class among d <= j:\n");
    long double cm=INFINITY; uint64_t cc=0; int cs=0;
    for(int j=0;j<=R+1;j++){ if(mr[j]<cm){ cm=mr[j]; cc=mc[j]; cs=ms[j]; }
      if(cnt[j]) printf("   d<=%2d: worst best-rate %.5Lf nats/move (class %llu, at length %d)\n",j,cm,(unsigned long long)cc,cs); }
  }
  // hostile threads: class of 1 + 3^k u (u = 1, 2 mod 3) and of (1 + 3^k w)/2 (w = 1, 2 mod 3)
  printf("  thread 1: classes 1 + 3^k u mod 3^%d, certified factor (worst over u mod 3^(R+1-k)):\n",R+1);
  for(int k=1;k<=R;k++){
    long double w=-INFINITY; SK ws={0,0}; int allok=1;
    for(uint64_t u=1;u<P[R+1-k];u++){ if(u%3==0) continue; uint64_t c=(1+P[k]*u)%M; long double v=val(cur[c]); if(!(v<0)) allok=0; if(v>w){w=v;ws=cur[c];} }
    printf("   k=%2d  %s worst 2^%d/3^%d = %.6Lf\n",k,allok?"all certified":"SOME EXCEPTIONAL",ws.K,ws.s,expl(w));
  }
  printf("  thread 1/2: classes (1 + 3^k w)/2 mod 3^%d, certified factor (worst over w):\n",R+1);
  uint64_t inv2=(M+1)/2;
  for(int k=1;k<=R;k++){
    long double w=-INFINITY; SK ws={0,0}; uint64_t nbad=0, ntot=0;
    for(uint64_t u=1;u<P[R+1-k];u++){ if(u%3==0) continue; ntot++; uint64_t c=((1+P[k]*u)%M)*inv2%M; long double v=val(cur[c]); if(!(v<0)) nbad++; else if(v>w){w=v;ws=cur[c];} }
    printf("   k=%2d  exceptional %llu/%llu, worst certified 2^%d/3^%d = %.6Lf\n",k,(unsigned long long)nbad,(unsigned long long)ntot,ws.K,ws.s,expl(w));
  }
  // dump classes with best factor >= 1/2 (certified or not): class, s, K, distance d
  if(argc>2){
    FILE*f=fopen(argv[2],"w");
    for(uint64_t c=0;c<M;c++){ if(c%3==0) continue; long double v=val(cur[c]);
      if(!(v<logl(0.5L))) fprintf(f,"%llu %d %d %d %s\n",(unsigned long long)c,cur[c].s,cur[c].K,dist[c],(v<0)?"cert":"EXC"); }
    fclose(f);
  }
  return 0;
}
