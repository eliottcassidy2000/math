// Althofer 3n+-1 game (= Conway's Beans-Don't-Talk on odd positions): no-draw verification in
// BLOCKS, for ranges too large for the single-batch exact-height run (game_sparse.c).
//
// PHASE 1 (dense): as in game_dense.c / game_sparse.c (incremental bottleneck retrograde up to CAP,
//   2 bits per odd non-multiple of 3).
// PHASE 2 (sparse, block mode): odd starts in [B1,B2) are processed in consecutive blocks.  For each
//   block the unresolved starts are targets; positions > CAP are inserted on demand into a hash map H
//   (candidates in increasing order via a min-heap), exactly as in game_sparse.c, with the explicit
//   "both children N" rule (no counters).  Explored unresolved table positions are marked with the
//   (otherwise unused in phase 2) table state 1.  When H exceeds HLIMIT, and at every block start,
//   H, the heap and the marks are discarded ("flush"); every label already written into the table is
//   a sound fact and is kept.  All labels come from the sound rules
//      N <- some child P,      P <- both children N,
//   so every resolved start has finite remoteness.  Heights are NOT exact in this mode: the reported
//   maximum popped value Tmax is an upper bound for the height of every resolved start.
//   If the heap empties with a target unresolved, the reachable unresolved set is closed (certified
//   finite DRAW trap) -- reported.  If two consecutive flushes resolve no target, the run aborts
//   (reported; no claim beyond the frontier).
//
// usage: game_blocks CAP B1 B2 [HLIMIT] [TMAXB]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

typedef unsigned __int128 u128;
static uint64_t *S; static uint64_t CAP, MAXIDX, WORDS;
static inline uint64_t idx_of(uint64_t n){ return n/3; }
static inline uint64_t n_of(uint64_t i){ return 3*i + 1 + (i&1); }
static inline int getS(uint64_t i){ return (S[i>>5] >> ((i&31)<<1)) & 3; }
static inline void setS(uint64_t i, int v){ uint64_t sh=(i&31)<<1; S[i>>5] = (S[i>>5] & ~(3ULL<<sh)) | ((uint64_t)v<<sh); }
static inline uint64_t oddpart(uint64_t x){ return x >> __builtin_ctzll(x); }
static inline uint64_t dchild(uint64_t n){ return oddpart((n&3)==1 ? 3*n+1 : 3*n-1); }
static inline uint64_t uchild(uint64_t n){ return ((n&3)==1 ? 3*n-1 : 3*n+1) >> 1; }

// hash map for positions > CAP: key -> state byte
static uint64_t *HK=NULL; static uint8_t *HV=NULL; static uint64_t Hcap=0,Hn=0,Hmask=0,Hmaxkey=0;
static inline uint64_t hsh(uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; x*=0xc4ceb9fe1a85ec53ULL; x^=x>>33; return x; }
static inline int64_t hfind(uint64_t x){ if(!Hcap) return -1; uint64_t h=hsh(x)&Hmask; while(HK[h]){ if(HK[h]==x) return (int64_t)h; h=(h+1)&Hmask; } return -1; }
static void hrehash(uint64_t ncap){
  uint64_t *ok=HK; uint8_t *ov=HV; uint64_t oc=Hcap;
  HK=calloc(ncap,8); HV=calloc(ncap,1); if(!HK||!HV){fprintf(stderr,"oom hash\n");exit(1);}
  Hcap=ncap; Hmask=ncap-1;
  for(uint64_t i=0;i<oc;i++) if(ok[i]){ uint64_t h=hsh(ok[i])&Hmask; while(HK[h]) h=(h+1)&Hmask; HK[h]=ok[i]; HV[h]=ov[i]; }
  free(ok); free(ov);
}
static inline void hput(uint64_t x, uint8_t v){
  if((Hn+1)*10 > Hcap*7) hrehash(Hcap? 2*Hcap : (1ULL<<20));
  uint64_t h=hsh(x)&Hmask; while(HK[h]){ if(HK[h]==x){HV[h]=v;return;} h=(h+1)&Hmask; }
  HK[h]=x; HV[h]=v; Hn++; if(x>Hmaxkey) Hmaxkey=x;
}
static void hclear(void){ if(Hcap){ memset(HK,0,Hcap*8); memset(HV,0,Hcap); } Hn=0; Hmaxkey=0; }

static inline int st(uint64_t y){ if(y<=CAP) return getS(idx_of(y)); int64_t h=hfind(y); return h<0?4:HV[h]; }

static uint64_t *Q; static size_t qn=0,qcap=0;
static inline void qpush(uint64_t x){ if(qn==qcap){ qcap=qcap?2*qcap:(1<<20); Q=realloc(Q,qcap*8); if(!Q){fprintf(stderr,"oom Q\n");exit(1);} } Q[qn++]=x; }
static uint64_t *HP; static size_t hpn=0,hpcap=0,hpmax=0;
static void hpush(uint64_t x){ if(hpn==hpcap){ hpcap=hpcap?2*hpcap:(1<<20); HP=realloc(HP,hpcap*8); if(!HP){fprintf(stderr,"oom heap\n");exit(1);} }
  size_t i=hpn++; while(i){ size_t p=(i-1)/2; if(HP[p]<=x) break; HP[i]=HP[p]; i=p; } HP[i]=x; if(hpn>hpmax) hpmax=hpn; }
static uint64_t hpop(void){ uint64_t top=HP[0], x=HP[--hpn]; size_t i=0; for(;;){ size_t l=2*i+1; if(l>=hpn) break; size_t c=(l+1<hpn && HP[l+1]<HP[l])? l+1:l; if(HP[c]>=x) break; HP[i]=HP[c]; i=c; } if(hpn) HP[i]=x; return top; }

static uint64_t now;
static inline uint64_t maxpresent(void){ return Hmaxkey>now? Hmaxkey: now; }
static inline void setres(uint64_t y,int v){ if(y<=CAP) setS(idx_of(y),v); else { int64_t h=hfind(y); HV[h]=(uint8_t)v; } qpush(y); }

static void propagate1(void){
  while(qn){
    uint64_t x=Q[--qn]; int v=getS(idx_of(x));
    u128 y=(u128)x;
    for(int k=1;k<100;k++){
      y<<=1; if(y > (u128)3*now+1) break;
      uint64_t yy=(uint64_t)y; uint64_t p=(yy%3==1)?(yy-1)/3:(yy+1)/3;
      if(p>now) break; if(p%3==0 || p==x) continue;
      uint64_t pi=idx_of(p); int sp=getS(pi); if(sp>=2) continue;
      if(v==3) setres(p,2); else { if(sp==0) setS(pi,1); else setres(p,3); }
    }
  }
}
static void propagate2(void){
  while(qn){
    uint64_t x=Q[--qn]; int v=st(x);
    uint64_t lim=maxpresent(); u128 y=(u128)x;
    for(int k=1;k<120;k++){
      y<<=1; if(y > (u128)3*lim+1) break;
      uint64_t yy=(uint64_t)y; uint64_t p=(yy%3==1)?(yy-1)/3:(yy+1)/3;
      if(p>lim) break; if(p%3==0 || p==x) continue;
      int sp=st(p); if(sp>=2) continue;
      if(v==3){ setres(p,2); continue; }
      uint64_t a=dchild(p), b=uchild(p); uint64_t o=(a==x)? b:a; int so=(o==1)?3:st(o);
      if(so==2) setres(p,3);
    }
  }
}
static int status_start(uint64_t n){
  if(n==1) return 2;
  if(n%3){ int s=st(n); return (s==2||s==3)? s:0; }
  uint64_t a=dchild(n), b=uchild(n); int sa=(a==1)?3:st(a), sb=st(b);
  if(sa==3||sb==3) return 2; if(sa==2&&sb==2) return 3; return 0;
}
static uint64_t *XS; static size_t xsn=0,xscap=0;
static void xpush(uint64_t x){ if(xsn==xscap){ xscap=xscap?2*xscap:(1<<16); XS=realloc(XS,xscap*8);} XS[xsn++]=x; }
static uint64_t nmarked=0;
static void explore(uint64_t y0){
  xpush(y0);
  while(xsn){
    uint64_t y=XS[--xsn]; if(y==1) continue;
    if(y<=CAP){ uint64_t i=idx_of(y); int s=getS(i); if(s!=0) continue; setS(i,1); nmarked++; xpush(dchild(y)); xpush(uchild(y)); }
    else { if(hfind(y)<0) hpush(y); }
  }
}
static int relevant(uint64_t x){
  uint64_t lim=maxpresent(); u128 y=(u128)x;
  for(int k=1;k<120;k++){
    y<<=1; if(y > (u128)3*lim+1) break;
    uint64_t yy=(uint64_t)y; uint64_t p=(yy%3==1)?(yy-1)/3:(yy+1)/3;
    if(p>lim) break; if(p%3==0) continue;
    int s=st(p); if(p<=CAP){ if(s==1) return 1; } else if(s<2) return 1;
  }
  return 0;
}
static void insert_sparse(uint64_t x){
  uint64_t a=dchild(x), b=uchild(x); int sa=(a==1)?3:st(a), sb=st(b);
  int v; if(sa==3||sb==3) v=2; else if(sa==2&&sb==2) v=3; else v=0;
  hput(x,0);
  if(v){ setres(x,v); propagate2(); }
  else { if(sa<2||sa==4) explore(a); if(sb<2||sb==4) explore(b); }
}
static void unmark_all(void){ // table state 01 -> 00
  for(uint64_t w=0; w<WORDS; w++){ uint64_t x=S[w]; uint64_t lo=x&0x5555555555555555ULL, hi=(x>>1)&0x5555555555555555ULL; uint64_t ones=lo&~hi; if(ones) S[w]=x&~ones; }
  nmarked=0;
}

static double pow10i(double x){ double p=1; while(p*10<=x) p*=10; return p; }
int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage: %s CAP B1 B2 [HLIMIT] [TMAXB]\n",argv[0]);return 1;}
  CAP=strtoull(argv[1],0,10); if(CAP%2==0) CAP--;
  uint64_t B1=strtoull(argv[2],0,10), B2=strtoull(argv[3],0,10);
  uint64_t HLIMIT= argc>4? strtoull(argv[4],0,10) : 11000000ULL;
  uint64_t TMAXB = argc>5? strtoull(argv[5],0,10) : 400000ULL;
  if(B2 > (CAP/3)*2){ fprintf(stderr,"need B2 <= 2CAP/3\n"); return 1; }
  MAXIDX=idx_of(CAP)+1; WORDS=(MAXIDX+31)/32+1;
  S=calloc(WORDS,8); if(!S){fprintf(stderr,"oom S\n");return 1;}
  time_t t0=time(0);
  printf("# game_blocks CAP=%llu B1=%llu B2=%llu HLIMIT=%llu TMAXB=%llu table=%.1fMB\n",(unsigned long long)CAP,(unsigned long long)B1,(unsigned long long)B2,(unsigned long long)HLIMIT,(unsigned long long)TMAXB,WORDS*8/1e6);
  setS(0,3);
  uint64_t f=1;
  for(uint64_t i=1;i<MAXIDX;i++){
    uint64_t m=n_of(i); if(m>CAP) break; now=m;
    int sa=getS(idx_of(dchild(m)));
    if(sa==3) setres(m,2); else if(sa==2) setS(i,1);
    propagate1();
    while(f<=m && status_start(f)) f+=2;
  }
  printf("PHASE1 FINAL CAP=%llu frontier=%llu elapsed=%lds\n",(unsigned long long)CAP,(unsigned long long)f,(long)(time(0)-t0)); fflush(stdout);
  unmark_all();   // phase 2 uses state 1 as the "explored" mark
  if(B1<f) B1=f; if(B1%2==0) B1++;
  uint64_t Tmax=CAP, totT=0, totflush=0, totpops=0; int aborted=0, trap=0;
  uint64_t *T=NULL; size_t tcap=0;
  uint64_t b=B1, BLK=1ULL<<26;
  uint64_t Hlim_alloc=1; while(Hlim_alloc*7 < (HLIMIT+2)*10) Hlim_alloc<<=1; hrehash(Hlim_alloc);
  printf("# hash slots=%llu (%.1f MB)\n",(unsigned long long)Hlim_alloc,Hlim_alloc*9/1e6);
  while(b<B2 && !aborted && !trap){
    // choose block size adaptively: at most TMAXB targets
    uint64_t e; size_t tn;
    for(;;){
      e=b+BLK; if(e>B2) e=B2; tn=0;
      for(uint64_t n=b;n<e;n+=2) if(!status_start(n)) tn++;
      if(tn<=TMAXB || BLK<=4096) break;
      BLK/=2;
    }
    if(tn>tcap){ tcap=tn+1024; T=realloc(T,tcap*8); if(!T){fprintf(stderr,"oom T\n");return 1;} }
    tn=0; for(uint64_t n=b;n<e;n+=2) if(!status_start(n)) T[tn++]=n;
    totT+=tn;
    // fresh start for the block
    hclear(); hpn=0; if(nmarked) unmark_all(); now=CAP;
    for(size_t j=0;j<tn;j++){ uint64_t n=T[j]; if(n%3) explore(n); else { explore(dchild(n)); explore(uchild(n)); } }
    size_t tp=0; uint64_t pops=0; int fl=0; size_t prev_resolved=(size_t)-1;
    uint64_t blockmax=CAP;
    while(tp<tn){
      while(tp<tn && status_start(T[tp])) tp++;
      if(tp>=tn) break;
      if(!hpn){ printf("HEAP EMPTY with target %llu unresolved: certified finite DRAW trap\n",(unsigned long long)T[tp]); trap=1; break; }
      uint64_t x=hpop();
      if(hfind(x)>=0) continue;
      if(!relevant(x)) continue;
      if(x>now) now=x; if(now>blockmax) blockmax=now;
      if(x > (1ULL<<61)){ printf("OVERFLOW GUARD: position %llu too large; abort\n",(unsigned long long)x); aborted=1; break; }
      insert_sparse(x); pops++;
      if(Hn>=HLIMIT){
        // number of resolved targets of this block so far; abort if a whole batch resolved none
        size_t res=tp; for(size_t j=tp;j<tn;j++) if(status_start(T[j])) res++;
        if(prev_resolved!=(size_t)-1 && res==prev_resolved){ printf("NO PROGRESS between flushes in block [%llu,%llu): first unresolved %llu; abort\n",(unsigned long long)b,(unsigned long long)e,(unsigned long long)T[tp]); aborted=1; break; }
        prev_resolved=res; fl++;
        hclear(); hpn=0; unmark_all(); now=CAP;
        for(size_t j=tp;j<tn;j++){ uint64_t n=T[j]; if(status_start(n)) continue; if(n%3) explore(n); else { explore(dchild(n)); explore(uchild(n)); } }
      }
    }
    if(blockmax>Tmax) Tmax=blockmax;
    totflush+=fl; totpops+=pops;
    if(aborted||trap){ printf("STOP in block [%llu,%llu)\n",(unsigned long long)b,(unsigned long long)e); break; }
    printf("BLOCK [%llu,%llu) targets=%zu flushes=%d pops=%llu blockmax=%llu H=%llu heapmax=%zu BLK=%llu elapsed=%lds\n",(unsigned long long)b,(unsigned long long)e,tn,fl,(unsigned long long)pops,(unsigned long long)blockmax,(unsigned long long)Hn,hpmax,(unsigned long long)BLK,(long)(time(0)-t0));
    fflush(stdout);
    b=e;
    if(tn<TMAXB/4 && BLK<(1ULL<<30)) BLK*=2;
  }
  // frontier and counts
  uint64_t F=1; while(F<B2 && status_start(F)) F+=2;
  uint64_t cN=0,cP=0; double nextck=10;
  printf("# P-counts among odd n < X (X = 10^k, 2*10^k, 5*10^k): X  N  P  P/(N+P)\n");
  for(uint64_t n=1;n<F;n+=2){
    if((double)n>nextck){ printf("PCOUNT %.0f %llu %llu %.6f\n",nextck,(unsigned long long)cN,(unsigned long long)cP,(double)cP/(cN+cP)); double l=nextck; nextck = (l==1*pow10i(l))? 2*l : (l==2*pow10i(l))? 5*l/2 : 2*l; }
    int s=status_start(n); if(s==2)cN++; else cP++; }
  // deterministic pseudo-random sample of resolved starts in [10^9, F) for independent certificate checks
  { uint64_t lo=1000000000ULL; if(F>lo+2){ uint64_t x=88172645463325252ULL; printf("# SAMPLE (xorshift64, 64 odd starts in [10^9,%llu)): n value\n",(unsigned long long)F);
      for(int i=0;i<64;i++){ x^=x<<13; x^=x>>7; x^=x<<17; uint64_t n=lo+(x%(F-lo)); n|=1; if(n>=F) n-=2; printf("SAMPLE %llu %c\n",(unsigned long long)n, status_start(n)==2?'N':'P'); } } }
  printf("DONE: every odd start < %llu resolved (finite remoteness, no draw); N=%llu P=%llu (odd n<%llu, n=1 counted N); total phase-2 targets=%llu flushes=%llu pops=%llu; height upper bound Tmax=%llu; aborted=%d trap=%d elapsed=%lds\n",
    (unsigned long long)F,(unsigned long long)cN,(unsigned long long)cP,(unsigned long long)F,(unsigned long long)totT,(unsigned long long)totflush,(unsigned long long)totpops,(unsigned long long)Tmax,aborted,trap,(long)(time(0)-t0));
  return 0;
}
