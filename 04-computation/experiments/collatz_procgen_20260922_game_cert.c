// Proof-certificate extractor for Althofer's 3n+-1 game (Beans-Don't-Talk on odd positions).
// Runs the same two-phase least-fixpoint computation as game_sparse.c (phase 1 dense up to CAP,
// phase 2 on-demand above CAP, single batch, explicit rules), additionally recording for every
// resolved position the global resolution sequence number.  Then, for each target, it extracts a
// proof DAG:  N-node -> one child that is P (or 1) and was resolved earlier;  P-node -> both children
// (N, resolved earlier).  Nodes are written in increasing sequence number (a topological order), so a
// checker only has to verify the two local rules against EARLIER lines; acyclicity is automatic.
// Checker: collatz_procgen_20260922_game_check.py cert FILE (independent Python arithmetic).
//
// usage: game_cert CAP outfile target1 [target2 ...]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef unsigned __int128 u128;
static uint8_t *S; static uint32_t *SQ; static uint64_t CAP, MAXIDX;   // byte states + seq for table
static inline uint64_t idx_of(uint64_t n){ return n/3; }
static inline uint64_t n_of(uint64_t i){ return 3*i + 1 + (i&1); }
static inline uint64_t oddpart(uint64_t x){ return x >> __builtin_ctzll(x); }
static inline uint64_t dchild(uint64_t n){ return oddpart((n&3)==1 ? 3*n+1 : 3*n-1); }
static inline uint64_t uchild(uint64_t n){ return ((n&3)==1 ? 3*n-1 : 3*n+1) >> 1; }

static uint64_t *HK; static uint8_t *HV; static uint64_t *HS; static uint64_t Hcap=0,Hn=0,Hmask=0,Hmaxkey=0;
static inline uint64_t hsh(uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; x*=0xc4ceb9fe1a85ec53ULL; x^=x>>33; return x; }
static inline int64_t hfind(uint64_t x){ if(!Hcap) return -1; uint64_t h=hsh(x)&Hmask; while(HK[h]){ if(HK[h]==x) return (int64_t)h; h=(h+1)&Hmask; } return -1; }
static void hrehash(uint64_t nc){ uint64_t *ok=HK,*os=HS; uint8_t *ov=HV; uint64_t oc=Hcap; HK=calloc(nc,8); HV=calloc(nc,1); HS=calloc(nc,8); if(!HK||!HV||!HS){fprintf(stderr,"oom\n");exit(1);} Hcap=nc; Hmask=nc-1;
  for(uint64_t i=0;i<oc;i++) if(ok[i]){ uint64_t h=hsh(ok[i])&Hmask; while(HK[h]) h=(h+1)&Hmask; HK[h]=ok[i]; HV[h]=ov[i]; HS[h]=os[i]; } free(ok); free(ov); free(os); }
static inline void hput(uint64_t x){ if((Hn+1)*10>Hcap*7) hrehash(Hcap?2*Hcap:(1ULL<<20)); uint64_t h=hsh(x)&Hmask; while(HK[h]){ if(HK[h]==x) return; h=(h+1)&Hmask; } HK[h]=x; HV[h]=0; HS[h]=0; Hn++; if(x>Hmaxkey) Hmaxkey=x; }
static inline int st(uint64_t y){ if(y<=CAP) return S[idx_of(y)]; int64_t h=hfind(y); return h<0?4:HV[h]; }
static inline uint64_t sq(uint64_t y){ if(y<=CAP) return SQ[idx_of(y)]; return HS[hfind(y)]; }

static uint64_t SEQ=1, now;
static uint64_t *Q; static size_t qn=0,qcap=0;
static inline void qpush(uint64_t x){ if(qn==qcap){ qcap=qcap?2*qcap:(1<<20); Q=realloc(Q,qcap*8);} Q[qn++]=x; }
static inline void setres(uint64_t y,int v){ if(y<=CAP){ S[idx_of(y)]=v; SQ[idx_of(y)]=(uint32_t)SEQ; } else { int64_t h=hfind(y); HV[h]=v; HS[h]=SEQ; } SEQ++; qpush(y); }
static inline uint64_t maxpresent(void){ return Hmaxkey>now? Hmaxkey: now; }
static void propagate(void){   // explicit rules (used in both phases)
  while(qn){
    uint64_t x=Q[--qn]; int v=st(x); uint64_t lim=maxpresent(); u128 y=(u128)x;
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
static uint64_t *HP; static size_t hpn=0,hpcap=0;
static void hpush(uint64_t x){ if(hpn==hpcap){ hpcap=hpcap?2*hpcap:(1<<20); HP=realloc(HP,hpcap*8);} size_t i=hpn++; while(i){ size_t p=(i-1)/2; if(HP[p]<=x) break; HP[i]=HP[p]; i=p; } HP[i]=x; }
static uint64_t hpop(void){ uint64_t top=HP[0], x=HP[--hpn]; size_t i=0; for(;;){ size_t l=2*i+1; if(l>=hpn) break; size_t c=(l+1<hpn && HP[l+1]<HP[l])? l+1:l; if(HP[c]>=x) break; HP[i]=HP[c]; i=c; } if(hpn) HP[i]=x; return top; }
static uint8_t *EX;  // explored marks for table nodes
static uint64_t *XS; static size_t xsn=0,xscap=0;
static void xpush(uint64_t x){ if(xsn==xscap){ xscap=xscap?2*xscap:(1<<16); XS=realloc(XS,xscap*8);} XS[xsn++]=x; }
static void explore(uint64_t y0){ xpush(y0); while(xsn){ uint64_t y=XS[--xsn]; if(y==1) continue; int s=st(y);
  if(y<=CAP){ if(s>=2) continue; uint64_t i=idx_of(y); if(EX[i]) continue; EX[i]=1; xpush(dchild(y)); xpush(uchild(y)); } else if(s==4) hpush(y); } }
static uint64_t *PIN; static int npin=0;   // targets and children of multiple-of-3 targets are always relevant
static int relevant(uint64_t x){ for(int j=0;j<npin;j++) if(PIN[j]==x) return 1; uint64_t lim=maxpresent(); u128 y=(u128)x;
  for(int k=1;k<120;k++){ y<<=1; if(y>(u128)3*lim+1) break; uint64_t yy=(uint64_t)y; uint64_t p=(yy%3==1)?(yy-1)/3:(yy+1)/3; if(p>lim) break; if(p%3==0) continue;
    int s=st(p); if(s<2){ if(p>CAP || EX[idx_of(p)]) return 1; } } return 0; }
static int status_start(uint64_t n){ if(n==1) return 2; if(n%3){ int s=st(n); return (s>=2&&s<4)? s:0; }
  uint64_t a=dchild(n), b=uchild(n); int sa=(a==1)?3:st(a), sb=st(b); if(sa==3||sb==3) return 2; if(sa==2&&sb==2) return 3; return 0; }

// certificate extraction
typedef struct { uint64_t seq, n, c; char v; } rec;
static rec *R; static size_t rn=0, rcap=0;
static uint64_t *VK; static uint64_t Vcap=0,Vn=0,Vmask=0;   // visited set
static int vadd(uint64_t x){ if((Vn+1)*10>Vcap*7){ uint64_t *o=VK, oc=Vcap; Vcap=Vcap?2*Vcap:(1<<16); Vmask=Vcap-1; VK=calloc(Vcap,8); for(uint64_t i=0;i<oc;i++) if(o[i]){ uint64_t h=hsh(o[i])&Vmask; while(VK[h]) h=(h+1)&Vmask; VK[h]=o[i]; } free(o); }
  uint64_t h=hsh(x)&Vmask; while(VK[h]){ if(VK[h]==x) return 0; h=(h+1)&Vmask; } VK[h]=x; Vn++; return 1; }
static void emit(uint64_t n, char v, uint64_t c, uint64_t s){ if(rn==rcap){ rcap=rcap?2*rcap:(1<<16); R=realloc(R,rcap*sizeof(rec)); } R[rn].n=n; R[rn].v=v; R[rn].c=c; R[rn].seq=s; rn++; }
static void extract(uint64_t t){
  xsn=0; xpush(t);
  while(xsn){
    uint64_t x=XS[--xsn]; if(x==1) continue; if(!vadd(x)) continue;
    uint64_t a=dchild(x), b=uchild(x);
    if(x%3==0){ // start that is a multiple of 3: justify from children (seq = max+1 conceptually)
      int s=status_start(x); int sa=(a==1)?3:st(a);
      if(s==2){ uint64_t c = (sa==3)? a : b; emit(x,'N',c,UINT64_MAX-1); if(c!=1) xpush(c); }
      else { emit(x,'P',0,UINT64_MAX-1); xpush(a); xpush(b); }
      continue;
    }
    int v=st(x); uint64_t s=sq(x);
    if(v==2){
      uint64_t c=0;
      if(a==1) c=1;
      else { uint64_t best=UINT64_MAX; int sa=st(a), sb=st(b);
        if(sa==3 && sq(a)<s && sq(a)<best){ best=sq(a); c=a; }
        if(sb==3 && sq(b)<s && sq(b)<best){ best=sq(b); c=b; } }
      if(!c){ fprintf(stderr,"internal: no earlier P child for %llu\n",(unsigned long long)x); exit(1); }
      emit(x,'N',c,s); if(c!=1) xpush(c);
    } else if(v==3){
      if(!(st(a)==2 && st(b)==2 && sq(a)<s && sq(b)<s)){ fprintf(stderr,"internal: P node %llu\n",(unsigned long long)x); exit(1); }
      emit(x,'P',0,s); xpush(a); xpush(b);
    } else { fprintf(stderr,"internal: unresolved %llu in proof\n",(unsigned long long)x); exit(1); }
  }
}
static int cmp(const void*p,const void*q){ const rec*a=p,*b=q; return a->seq<b->seq? -1 : a->seq>b->seq; }

int main(int argc,char**argv){
  if(argc<4){ fprintf(stderr,"usage: %s CAP outfile targets...\n",argv[0]); return 1; }
  CAP=strtoull(argv[1],0,10); if(CAP%2==0) CAP--;
  const char* out=argv[2]; int nt=argc-3; uint64_t *T=malloc(nt*8);
  for(int j=0;j<nt;j++) T[j]=strtoull(argv[3+j],0,10);
  MAXIDX=idx_of(CAP)+1; S=calloc(MAXIDX,1); SQ=calloc(MAXIDX,4); EX=calloc(MAXIDX,1);
  S[0]=3; SQ[0]=0;
  for(uint64_t i=1;i<MAXIDX;i++){ uint64_t m=n_of(i); if(m>CAP) break; now=m;
    uint64_t a=dchild(m), b=uchild(m); int sa=(a==1)?3:st(a); (void)b;
    if(sa==3) setres(m,2);
    propagate(); }
  if(SEQ>=UINT32_MAX){ fprintf(stderr,"seq overflow\n"); return 1; }
  PIN=malloc(3*nt*8); for(int j=0;j<nt;j++){ PIN[npin++]=T[j]; if(T[j]%3==0){ PIN[npin++]=dchild(T[j]); PIN[npin++]=uchild(T[j]); } }
  for(int j=0;j<nt;j++){ uint64_t n=T[j]; if(status_start(n)) continue; if(n%3) explore(n); else { explore(dchild(n)); explore(uchild(n)); } }
  uint64_t pops=0;
  for(;;){ int all=1; for(int j=0;j<nt;j++) if(!status_start(T[j])){ all=0; break; } if(all) break;
    if(!hpn){ printf("HEAP EMPTY: unresolved closed set (DRAW trap)\n"); return 2; }
    uint64_t x=hpop(); if(st(x)!=4) continue; if(!relevant(x)) continue; if(x>now) now=x;
    uint64_t a=dchild(x), b=uchild(x); int sa=(a==1)?3:st(a), sb=st(b);
    hput(x); pops++;
    if(sa==3||sb==3) setres(x,2); else if(sa==2&&sb==2) setres(x,3);
    if(qn) propagate(); else { if(sa<2||sa==4) explore(a); if(sb<2||sb==4) explore(b); }
  }
  for(int j=0;j<nt;j++) extract(T[j]);
  qsort(R,rn,sizeof(rec),cmp);
  FILE*F=fopen(out,"w"); fprintf(F,"# Althofer 3n+-1 game proof certificate (game_cert.c, CAP=%llu). Lines: n V [c], topological order.\n# targets:",(unsigned long long)CAP);
  for(int j=0;j<nt;j++) fprintf(F," %llu",(unsigned long long)T[j]); fprintf(F,"\n");
  uint64_t mx=0; for(size_t r=0;r<rn;r++){ if(R[r].n>mx) mx=R[r].n; if(R[r].v=='N') fprintf(F,"%llu N %llu\n",(unsigned long long)R[r].n,(unsigned long long)R[r].c); else fprintf(F,"%llu P\n",(unsigned long long)R[r].n); }
  fclose(F);
  printf("certificate %s: %zu nodes, largest node %llu, phase-2 insertions %llu, max inserted %llu\n",out,rn,(unsigned long long)mx,(unsigned long long)pops,(unsigned long long)now);
  for(int j=0;j<nt;j++) printf("  target %llu: %s\n",(unsigned long long)T[j], status_start(T[j])==2?"N (mover wins)":"P (mover loses)");
  return 0;
}
