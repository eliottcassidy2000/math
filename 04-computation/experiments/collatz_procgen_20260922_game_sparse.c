// Althofer 3n+-1 game (= Conway's Beans-Don't-Talk restricted to odd positions).
// Position: odd n.  Moves: n -> oddpart(3n+1), n -> oddpart(3n-1).  Moving to 1 wins.
//
// PHASE 1 (dense): incremental bottleneck retrograde on all odd non-multiples of 3 up to CAP
//   (same algorithm as collatz_procgen_20260922_game_dense.c; 2 bits per position; positions are
//   inserted in increasing order, so the insertion time at which a position resolves is its exact
//   HEIGHT h(n) = least cap admitting a finite win/loss proof tree all of whose nodes are <= cap).
// PHASE 2 (sparse, "climb only where needed"): targets = odd starts in [B1,B2) left unresolved by
//   phase 1.  We continue the SAME least-fixpoint computation on [1,CAP] u H, H = hash map of
//   positions > CAP inserted on demand.  A position x > CAP is a candidate when it is a child of a
//   present unresolved position that is reachable from an unresolved target through unresolved
//   positions; candidates are inserted in increasing order (min-heap), so resolution times are
//   heights (checked against phase-1 heights on overlaps).  Candidates whose present parents are all
//   resolved are skipped (re-pushed if a new unresolved parent appears).  Phase-2 propagation uses an
//   explicit "both children N" check (no counters), so H may be flushed when it exceeds HLIMIT
//   (facts already written into the table stay valid; after a flush reported times are only bounds).
//   Every label comes from the sound rules N <- (some child P), P <- (both children N); hence every
//   resolved start has finite remoteness (is not a draw).  If the heap empties with a target
//   unresolved, the unresolved reachable set is closed: a certified finite DRAW trap (reported).
//
// usage: game_sparse CAP B1 B2 [HLIMIT] [timesfile]
//   timesfile: binary pairs (uint64 n, uint64 time) for every target resolved in phase 2.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

typedef unsigned __int128 u128;
static uint64_t *S; static uint64_t CAP, MAXIDX;
static inline uint64_t idx_of(uint64_t n){ return n/3; }
static inline uint64_t n_of(uint64_t i){ return 3*i + 1 + (i&1); }
static inline int getS(uint64_t i){ return (S[i>>5] >> ((i&31)<<1)) & 3; }
static inline void setS(uint64_t i, int v){ uint64_t sh=(i&31)<<1; S[i>>5] = (S[i>>5] & ~(3ULL<<sh)) | ((uint64_t)v<<sh); }
static inline uint64_t oddpart(uint64_t x){ return x >> __builtin_ctzll(x); }
static inline uint64_t dchild(uint64_t n){ return oddpart((n&3)==1 ? 3*n+1 : 3*n-1); }
static inline uint64_t uchild(uint64_t n){ return ((n&3)==1 ? 3*n-1 : 3*n+1) >> 1; }

// ---------------- generic open-addressing hash (key -> value) ----------------
typedef struct { uint64_t *K; uint64_t *V; uint8_t *B; uint64_t cap,n,mask; int kind; } hmap; // kind 0: byte values, 1: u64 values, 2: set
static inline uint64_t hsh(uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; x*=0xc4ceb9fe1a85ec53ULL; x^=x>>33; return x; }
static inline int64_t hm_find(hmap*m,uint64_t x){ if(!m->cap) return -1; uint64_t h=hsh(x)&m->mask; while(m->K[h]){ if(m->K[h]==x) return (int64_t)h; h=(h+1)&m->mask; } return -1; }
static void hm_free(hmap*m){ free(m->K); free(m->V); free(m->B); m->K=0;m->V=0;m->B=0;m->cap=m->n=m->mask=0; }
static void hm_rehash(hmap*m,uint64_t ncap){
  hmap o=*m; m->cap=ncap; m->mask=ncap-1; m->n=0;
  m->K=calloc(ncap,8); if(m->kind==0) m->B=calloc(ncap,1); if(m->kind==1) m->V=calloc(ncap,8);
  if(!m->K || (m->kind==0&&!m->B) || (m->kind==1&&!m->V)){ fprintf(stderr,"oom hash rehash to %llu\n",(unsigned long long)ncap); exit(1); }
  for(uint64_t i=0;i<o.cap;i++) if(o.K[i]){ uint64_t h=hsh(o.K[i])&m->mask; while(m->K[h]) h=(h+1)&m->mask; m->K[h]=o.K[i]; if(m->kind==0) m->B[h]=o.B[i]; if(m->kind==1) m->V[h]=o.V[i]; m->n++; }
  free(o.K); free(o.V); free(o.B);
}
static inline uint64_t hm_slot(hmap*m,uint64_t x,int*isnew){
  if((m->n+1)*10 > m->cap*7) hm_rehash(m, m->cap? 2*m->cap : (1ULL<<16));
  uint64_t h=hsh(x)&m->mask; while(m->K[h]){ if(m->K[h]==x){ *isnew=0; return h; } h=(h+1)&m->mask; }
  m->K[h]=x; m->n++; *isnew=1; return h;
}
static hmap H={0,0,0,0,0,0,0}, E={0,0,0,0,0,0,2}, T2={0,0,0,0,0,0,1};
static uint64_t Hmaxkey=0;

// state: 0/1 unresolved, 2 N, 3 P, 4 absent
static inline int st(uint64_t y){ if(y<=CAP) return getS(idx_of(y)); int64_t h=hm_find(&H,y); return h<0?4:H.B[h]; }

static uint64_t *Q; static size_t qn=0,qcap=0;
static inline void qpush(uint64_t x){ if(qn==qcap){ qcap=qcap?2*qcap:(1<<20); Q=realloc(Q,qcap*8); if(!Q){fprintf(stderr,"oom Q\n");exit(1);} } Q[qn++]=x; }
static uint64_t *HP; static size_t hpn=0,hpcap=0,hpmax=0;
static void hpush(uint64_t x){ if(hpn==hpcap){ hpcap=hpcap?2*hpcap:(1<<20); HP=realloc(HP,hpcap*8); if(!HP){fprintf(stderr,"oom heap\n");exit(1);} }
  size_t i=hpn++; while(i){ size_t p=(i-1)/2; if(HP[p]<=x) break; HP[i]=HP[p]; i=p; } HP[i]=x; if(hpn>hpmax) hpmax=hpn; }
static uint64_t hpop(void){ uint64_t top=HP[0], x=HP[--hpn]; size_t i=0; for(;;){ size_t l=2*i+1; if(l>=hpn) break; size_t c=(l+1<hpn && HP[l+1]<HP[l])? l+1:l; if(HP[c]>=x) break; HP[i]=HP[c]; i=c; } if(hpn) HP[i]=x; return top; }

static uint64_t now, nres=0, LIM2=0; static int phase=1;
static inline uint64_t maxpresent(void){ return Hmaxkey>now? Hmaxkey: now; }

static inline void setres(uint64_t y,int v){
  if(y<=CAP){ setS(idx_of(y),v); if(phase==2 && y<LIM2){ int nw; uint64_t s=hm_slot(&T2,y,&nw); T2.V[s]=now; } }
  else { int64_t h=hm_find(&H,y); H.B[h]=(uint8_t)v; }
  qpush(y); nres++;
}

static void propagate1(void){   // phase 1: counters
  while(qn){
    uint64_t x=Q[--qn]; int v=getS(idx_of(x));
    u128 y=(u128)x;
    for(int k=1;k<100;k++){
      y<<=1; if(y > (u128)3*now+1) break;
      uint64_t yy=(uint64_t)y; uint64_t p=(yy%3==1)?(yy-1)/3:(yy+1)/3;
      if(p>now) break;
      if(p%3==0 || p==x) continue;
      uint64_t pi=idx_of(p); int sp=getS(pi); if(sp>=2) continue;
      if(v==3) setres(p,2); else { if(sp==0) setS(pi,1); else setres(p,3); }
    }
  }
}
static void propagate2(void){   // phase 2: explicit checks
  while(qn){
    uint64_t x=Q[--qn]; int v=st(x);
    uint64_t lim=maxpresent();
    u128 y=(u128)x;
    for(int k=1;k<120;k++){
      y<<=1; if(y > (u128)3*lim+1) break;
      uint64_t yy=(uint64_t)y; uint64_t p=(yy%3==1)?(yy-1)/3:(yy+1)/3;
      if(p>lim) break;
      if(p%3==0 || p==x) continue;
      int sp=st(p); if(sp>=2) continue;   // resolved or absent
      if(v==3){ setres(p,2); continue; }
      uint64_t a=dchild(p), b=uchild(p); uint64_t o=(a==x)? b : a;
      int so=(o==1)?3:st(o);
      if(so==2) setres(p,3);
    }
  }
}
static int status_start(uint64_t n){
  if(n==1) return 2;
  if(n%3){ int s=st(n); return (s==2||s==3)? s:0; }
  uint64_t a=dchild(n), b=uchild(n);
  int sa=(a==1)?3:st(a), sb=st(b);
  if(sa==3||sb==3) return 2;
  if(sa==2&&sb==2) return 3;
  return 0;
}
static uint64_t *XS; static size_t xsn=0,xscap=0;
static void xpush(uint64_t x){ if(xsn==xscap){ xscap=xscap?2*xscap:(1<<16); XS=realloc(XS,xscap*8);} XS[xsn++]=x; }
static void explore(uint64_t y0){
  xpush(y0);
  while(xsn){
    uint64_t y=XS[--xsn]; if(y==1) continue;
    int s=st(y);
    if(y<=CAP){ if(s>=2) continue; int nw; hm_slot(&E,y,&nw); if(!nw) continue; xpush(dchild(y)); xpush(uchild(y)); }
    else if(s==4) hpush(y);
  }
}
static int relevant(uint64_t x){
  uint64_t lim=maxpresent(); u128 y=(u128)x;
  for(int k=1;k<120;k++){
    y<<=1; if(y > (u128)3*lim+1) break;
    uint64_t yy=(uint64_t)y; uint64_t p=(yy%3==1)?(yy-1)/3:(yy+1)/3;
    if(p>lim) break; if(p%3==0) continue;
    int s=st(p); if(s<2){ if(p>CAP || hm_find(&E,p)>=0) return 1; }
  }
  return 0;
}
static void insert_sparse(uint64_t x){
  uint64_t a=dchild(x), b=uchild(x);
  int sa=(a==1)?3:st(a), sb=st(b);
  int v; if(sa==3||sb==3) v=2; else if(sa==2&&sb==2) v=3; else v=0;
  int nw; uint64_t s=hm_slot(&H,x,&nw); H.B[s]=0; if(x>Hmaxkey) Hmaxkey=x;
  if(v){ setres(x,v); propagate2(); }
  else { if(sa!=2&&sa!=3) explore(a); if(sb!=2&&sb!=3) explore(b); }
}
static uint64_t t2get(uint64_t y){ int64_t h=hm_find(&T2,y); return h<0? CAP : T2.V[h]; }  // phase-1 resolved => <= CAP
static uint64_t target_time(uint64_t n){
  if(n%3) return t2get(n);
  uint64_t a=dchild(n), b=uchild(n); int sa=(a==1)?3:st(a), sb=st(b);
  uint64_t ta=(a==1)?0:t2get(a), tb=t2get(b), best=UINT64_MAX;
  if(sa==3 && ta<best) best=ta; if(sb==3 && tb<best) best=tb;
  if(best==UINT64_MAX && sa==2 && sb==2) best= ta>tb? ta:tb;
  return best;
}

int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage: %s CAP B1 B2 [HLIMIT] [timesfile]\n",argv[0]);return 1;}
  CAP=strtoull(argv[1],0,10); if(CAP%2==0) CAP--;
  uint64_t B1=strtoull(argv[2],0,10), B2=strtoull(argv[3],0,10);
  uint64_t HLIMIT= argc>4? strtoull(argv[4],0,10) : 40000000ULL;
  const char* tfile= argc>5? argv[5]:NULL;
  if(B2 > (CAP/3)*2){ fprintf(stderr,"need B2 <= 2CAP/3\n"); return 1; }
  LIM2 = B2 + B2/2 + 2;
  MAXIDX=idx_of(CAP)+1; size_t words=(MAXIDX+31)/32+1;
  S=calloc(words,8); if(!S){fprintf(stderr,"oom S\n");return 1;}
  time_t t0=time(0);
  printf("# game_sparse CAP=%llu B1=%llu B2=%llu HLIMIT=%llu table=%.1fMB\n",(unsigned long long)CAP,(unsigned long long)B1,(unsigned long long)B2,(unsigned long long)HLIMIT,words*8/1e6);
  setS(0,3);
  uint64_t f=1, nextpow=1024;
  for(uint64_t i=1;i<MAXIDX;i++){
    uint64_t m=n_of(i); if(m>CAP) break; now=m;
    int sa=getS(idx_of(dchild(m)));
    if(sa==3) setres(m,2); else if(sa==2) setS(i,1);
    propagate1();
    while(f<=m && status_start(f)) f+=2;
    if(m>=nextpow){ printf("PHASE1 cap=2^%d frontier=%llu\n",63-__builtin_clzll(nextpow),(unsigned long long)f); nextpow<<=1; fflush(stdout); }
  }
  printf("PHASE1 FINAL CAP=%llu frontier=%llu elapsed=%lds\n",(unsigned long long)CAP,(unsigned long long)f,(long)(time(0)-t0)); fflush(stdout);
  phase=2;
  uint64_t *T=NULL; size_t tn=0,tcap=0;
  if(B1<f) B1=f; if(B1%2==0) B1++;
  for(uint64_t n=B1;n<B2;n+=2) if(!status_start(n)){ if(tn==tcap){tcap=tcap?2*tcap:(1<<16); T=realloc(T,tcap*8);} T[tn++]=n; }
  printf("PHASE2 targets (unresolved odd starts in [%llu,%llu)): %zu\n",(unsigned long long)B1,(unsigned long long)B2,tn); fflush(stdout);
  size_t tp=0; uint64_t pops=0,skips=0,Tmax=now; int nflush=0; time_t tl=time(0);
  for(size_t j=0;j<tn;j++){ uint64_t n=T[j]; if(n%3) explore(n); else { explore(dchild(n)); explore(uchild(n)); } }
  printf("PHASE2 initial explored table positions=%llu candidates=%zu\n",(unsigned long long)E.n,hpn); fflush(stdout);
  int trap=0;
  while(tp<tn){
    while(tp<tn && status_start(T[tp])) tp++;
    if(tp>=tn) break;
    if(!hpn){ printf("HEAP EMPTY with target %llu unresolved: closed unresolved set = certified finite DRAW trap\n",(unsigned long long)T[tp]); trap=1; break; }
    uint64_t x=hpop();
    if(st(x)!=4) continue;
    if(!relevant(x)){ skips++; continue; }
    if(x>now) now=x; if(now>Tmax) Tmax=now;
    insert_sparse(x); pops++;
    if(H.n>HLIMIT){
      nflush++; printf("FLUSH #%d: H=%llu E=%llu heap=%zu now=%llu first-unresolved-target=%llu (%zu/%zu)\n",nflush,(unsigned long long)H.n,(unsigned long long)E.n,hpn,(unsigned long long)now,(unsigned long long)T[tp],tp,tn); fflush(stdout);
      hm_free(&H); hm_free(&E); hpn=0; Hmaxkey=0; now=CAP;
      for(size_t j=tp;j<tn;j++){ uint64_t n=T[j]; if(status_start(n)) continue; if(n%3) explore(n); else { explore(dchild(n)); explore(uchild(n)); } }
    }
    if(time(0)-tl>=60){ tl=time(0); fprintf(stderr,"  pops=%llu skips=%llu now=%.3e H=%llu E=%llu heap=%zu tp=%zu/%zu first-unres=%llu el=%lds\n",(unsigned long long)pops,(unsigned long long)skips,(double)now,(unsigned long long)H.n,(unsigned long long)E.n,hpn,tp,tn,(unsigned long long)(tp<tn?T[tp]:0),(long)(time(0)-t0)); }
  }
  uint64_t Bfinal=(tp<tn)? T[tp] : B2;
  printf("PHASE2 DONE: every odd start in [%llu,%llu) resolved (N or P, finite remoteness)%s; targets %zu/%zu; max time used=%llu; flushes=%d; H=%llu E=%llu heapmax=%zu pops=%llu skips=%llu trap=%d elapsed=%lds\n",
    (unsigned long long)B1,(unsigned long long)Bfinal, (B1<=f+1)? " -- together with phase 1: every odd start below the right end":" (starts between the phase-1 frontier and B1 were not processed in this run)",tp,tn,(unsigned long long)Tmax,nflush,(unsigned long long)H.n,(unsigned long long)E.n,hpmax,(unsigned long long)pops,(unsigned long long)skips,trap,(long)(time(0)-t0));
  // per-target heights
  if(tp==tn){
    uint64_t best=0,bestn=0; double bestr=0; uint64_t bestrn=0;
    FILE*F= tfile? fopen(tfile,"wb"):NULL;
    // histogram of log2(time/n)
    uint64_t hist[64]={0};
    for(size_t j=0;j<tn;j++){ uint64_t n=T[j], t=target_time(n); if(F){ fwrite(&n,8,1,F); fwrite(&t,8,1,F);} if(t>best){best=t;bestn=n;} double r=(double)t/(double)n; if(r>bestr){bestr=r;bestrn=n;} int b=0; while(b<63 && (double)(1ULL<<(b+1))<=r) b++; hist[b]++; }
    if(F) fclose(F);
    printf("PHASE2 per-target: largest height %llu at start n=%llu ; largest height/n = %.1f at n=%llu%s\n",(unsigned long long)best,(unsigned long long)bestn,bestr,(unsigned long long)bestrn, nflush? " (times after a flush are not exact heights)":" (exact: no flush)");
    printf("PHASE2 histogram of floor(log2(height/n)) over phase-2 targets:"); for(int b=0;b<64;b++) if(hist[b]) printf(" %d:%llu",b,(unsigned long long)hist[b]); printf("\n");
  }
  return 0;
}
