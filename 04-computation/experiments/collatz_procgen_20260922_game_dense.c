// Althofer 3n+-1 game (= Conway's Beans-Don't-Talk on odd numbers).
// Position: odd n.  Moves: n -> oddpart(3n+1), n -> oddpart(3n-1).  Moving to 1 wins.
//
// INCREMENTAL BOTTLENECK RETROGRADE ("exact height").
//   Positions (odd, not divisible by 3 -- multiples of 3 are never children) are inserted in
//   increasing order m = 5,7,11,13,...  After inserting all positions <= m, the table holds
//   the least fixpoint of the retrograde rules on the induced subgraph [1,m] (absent children =
//   unknown).  A position n is resolved at time h(n) = the smallest cap m for which n has a
//   finite win/loss proof tree with all nodes <= m ("height" of n).  Each resolution event
//   propagates to the present predecessors p = (x*2^k -+ 1)/3 (k>=1).
//   State: 2 bits per position: 0=U0 (unresolved, no child known N), 1=U1 (one child known N),
//   2=N (mover wins), 3=P (mover loses).  Terminal 1 is stored as P (as a child it is a lost
//   position for the player to move); the START n=1 is N (3*1+1=4 -> 1).
//   Multiples of 3 are evaluated from their two children.
//
// Frontier f(t) = smallest odd start (including multiples of 3) unresolved at cap t.
// Output: frontier checkpoints (t, f), blocker records, final counts, first unresolved start.
// Optional: dump the 2-bit table for n <= DUMPN to a file, and h(n) for idx < HIDX.
//
// usage: game_dense CAP [HIDX] [DUMPN dumpfile] [hfile]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

typedef unsigned __int128 u128;
static uint64_t *S;         // 2-bit states, 32 per word
static uint64_t CAP, MAXIDX;
static inline uint64_t idx_of(uint64_t n){ return n/3; }          // n odd, 3 !| n
static inline uint64_t n_of(uint64_t i){ return 3*i + 1 + (i&1); }
static inline int getS(uint64_t i){ return (S[i>>5] >> ((i&31)<<1)) & 3; }
static inline void setS(uint64_t i, int v){ uint64_t sh=(i&31)<<1; S[i>>5] = (S[i>>5] & ~(3ULL<<sh)) | ((uint64_t)v<<sh); }
static inline uint64_t oddpart(uint64_t x){ return x >> __builtin_ctzll(x); }
static inline uint64_t dchild(uint64_t n){ return oddpart((n&3)==1 ? 3*n+1 : 3*n-1); }
static inline uint64_t uchild(uint64_t n){ return ((n&3)==1 ? 3*n-1 : 3*n+1) >> 1; }

static uint64_t *Q; static size_t qn=0, qcap=0;
static inline void push(uint64_t x){ if(qn==qcap){ qcap=qcap?2*qcap:(1<<20); Q=realloc(Q,qcap*8); if(!Q){fprintf(stderr,"oom\n");exit(1);} } Q[qn++]=x; }

static uint64_t HIDX=0; static uint64_t *Hh=NULL;   // h(n) for idx < HIDX
static uint64_t now;                                // current cap (last inserted position)
static uint64_t nres=0;

static inline void resolve(uint64_t n, int v){ uint64_t i=idx_of(n); setS(i,v); if(i<HIDX) Hh[i]=now; push(n); nres++; }

static void propagate(void){
  while(qn){
    uint64_t x=Q[--qn]; int v=getS(idx_of(x));
    // predecessors p = (x*2^k -+ 1)/3, k>=1, p <= now
    u128 y=(u128)x;
    for(int k=1;k<100;k++){
      y<<=1; if(y > (u128)3*now+1) break;
      uint64_t yy=(uint64_t)y; uint64_t p = (yy%3==1) ? (yy-1)/3 : (yy+1)/3;
      if(p>now) break;
      if(p%3==0) continue;
      if(p==x) continue;          // only x=1
      uint64_t pi=idx_of(p); int sp=getS(pi);
      if(sp>=2) continue;
      if(v==3){ resolve(p,2); }
      else { if(sp==0) setS(pi,1); else resolve(p,3); }
    }
  }
}

// status of an arbitrary odd start n (possibly multiple of 3) at current cap: 0 unresolved, 2 N, 3 P
static inline int status_start(uint64_t n){
  if(n==1) return 2;
  if(n%3) { int s=getS(idx_of(n)); return s>=2? s:0; }
  uint64_t a=dchild(n), b=uchild(n);
  int sa = (a==1)?3:getS(idx_of(a)); if(sa<2) sa=0;
  if(sa==3) return 2;
  int sb = (b<=now)? getS(idx_of(b)) : 0; if(sb<2) sb=0;
  if(sb==3) return 2;
  if(sa==2 && sb==2) return 3;
  return 0;
}

int main(int argc,char**argv){
  if(argc<2){fprintf(stderr,"usage: %s CAP [HIDX] [DUMPN dumpfile] [hfile]\n",argv[0]);return 1;}
  CAP=strtoull(argv[1],0,10); if(CAP%2==0) CAP--;
  HIDX = argc>2? strtoull(argv[2],0,10):0;
  uint64_t DUMPN = argc>4? strtoull(argv[3],0,10):0; const char* dumpf = argc>4? argv[4]:NULL;
  const char* hfile = argc>5? argv[5]:NULL;
  MAXIDX = idx_of(CAP)+1;
  size_t words = (MAXIDX+31)/32 + 1;
  S = calloc(words,8); if(!S){fprintf(stderr,"oom S\n");return 1;}
  if(HIDX){ if(HIDX>MAXIDX) HIDX=MAXIDX; Hh=calloc(HIDX,8); if(!Hh){fprintf(stderr,"oom H\n");return 1;} }
  fprintf(stderr,"CAP=%llu table=%.1f MB  H=%.1f MB\n",(unsigned long long)CAP, words*8/1e6, HIDX*8/1e6);
  time_t t0=time(0);
  setS(0,3); // terminal 1 as a child = P
  if(HIDX) Hh[0]=1;
  uint64_t f=1; // frontier over odd starts: smallest odd start unresolved
  // checkpoints: record when f first exceeds thresholds
  static const double mant[3]={1,2,5}; int ckj=1, ckm=0; double nextck=10; // 10,20,50,100,...
  uint64_t nextpow=1024; // report f at caps 2^j
  double recratio=0; // blocker record h(n)/n
  printf("# frontier checkpoints: cap t at which every odd start < f is resolved (f crosses threshold)\n");
  for(uint64_t i=1;i<MAXIDX;i++){
    uint64_t m=n_of(i); if(m>CAP) break;
    now=m;
    uint64_t a=dchild(m);
    int sa=getS(idx_of(a));
    if(sa==3){ resolve(m,2); }
    else if(sa==2){ setS(i,1); }
    // else U0 (already 0)
    propagate();
    // advance frontier
    uint64_t fold=f;
    while(f<=m){ int st=status_start(f); if(!st) break; f+=2; }
    if(f!=fold){
      // blocker = fold (it was resolved just now at time m); record ratio
      double r=(double)m/(double)fold;
      if(fold>1000 && r>recratio){ recratio=r; printf("BLOCKREC n=%llu h=%llu h/n=%.2f newfrontier=%llu\n",(unsigned long long)fold,(unsigned long long)m,r,(unsigned long long)f); }
      while((double)f>nextck){ printf("CK f>%.0f at cap t=%llu (f=%llu)\n",nextck,(unsigned long long)m,(unsigned long long)f); ckm++; if(ckm==3){ckm=0;ckj++;} double p10=1; for(int q=0;q<ckj;q++) p10*=10; nextck=mant[ckm]*p10; }
    }
    if(m>=nextpow){ printf("CAPPOW cap=2^%d (m=%llu): frontier B=%llu\n",63-__builtin_clzll(nextpow),(unsigned long long)m,(unsigned long long)f); nextpow<<=1; }
    if((i & ((1ULL<<30)-1))==0){ fprintf(stderr,"  i=%llu m=%llu f=%llu resolved=%llu elapsed=%lds\n",(unsigned long long)i,(unsigned long long)m,(unsigned long long)f,(unsigned long long)nres,(long)(time(0)-t0)); }
  }
  printf("FINAL CAP=%llu first unresolved odd start B=%llu ; elapsed %lds\n",(unsigned long long)CAP,(unsigned long long)f,(long)(time(0)-t0));
  // counts
  uint64_t cnt[4]={0,0,0,0};
  for(uint64_t i=1;i<MAXIDX;i++){ if(n_of(i)>CAP) break; cnt[getS(i)]++; }
  printf("non-multiples-of-3 in [5,CAP]: N=%llu P=%llu unresolved=%llu (U1=%llu)\n",(unsigned long long)cnt[2],(unsigned long long)cnt[3],(unsigned long long)(cnt[0]+cnt[1]),(unsigned long long)cnt[1]);
  // all odd 3<=n<=CAP
  uint64_t aN=0,aP=0,aU=0; for(uint64_t n=3;n<=CAP;n+=2){ int s=status_start(n); if(s==2)aN++; else if(s==3)aP++; else aU++; }
  printf("all odd 3<=n<=CAP: N=%llu P=%llu unresolved=%llu\n",(unsigned long long)aN,(unsigned long long)aP,(unsigned long long)aU);
  if(dumpf && DUMPN){
    uint64_t di=idx_of(DUMPN)+1; if(di>MAXIDX) di=MAXIDX; size_t dw=(di+31)/32;
    FILE*F=fopen(dumpf,"wb"); fwrite(&CAP,8,1,F); fwrite(&DUMPN,8,1,F); fwrite(S,8,dw,F); fclose(F);
    fprintf(stderr,"dumped %zu words\n",dw);
  }
  if(hfile && HIDX){ FILE*F=fopen(hfile,"wb"); fwrite(&HIDX,8,1,F); fwrite(Hh,8,HIDX,F); fclose(F); }
  return 0;
}
