// Structure statistics for Althofer's 3n+-1 game (Beans-Don't-Talk on odd positions).
// 1. Computes the exact value (N = mover wins, P = mover loses) of every odd n < NB using the
//    two-phase least-fixpoint method of game_sparse.c (dense to CAP, on-demand above CAP).
// 2. Prints statistics on the fully determined range odd n < BA (BA = 2*NB/3 so that both children of
//    every analysed n are determined):
//    P-density (cumulative, per octave, per mantissa bin), residues mod 2^k and 3^j, the r <-> -r
//    symmetry, descent exponent, move structure of N-positions (which moves win), the down-only game
//    (parity of the descending chain), sheet statistics, and reflection tests.
// usage: game_struct CAP BA
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
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
static inline int dexp(uint64_t n){ return __builtin_ctzll((n&3)==1 ? 3*n+1 : 3*n-1); }
static inline uint64_t Uplus(uint64_t n){ return oddpart(3*n+1); }
static inline uint64_t Uminus(uint64_t n){ return oddpart(3*n-1); }

static uint64_t *HK; static uint8_t *HV; static uint64_t Hcap=0,Hn=0,Hmask=0,Hmaxkey=0;
static inline uint64_t hsh(uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; x*=0xc4ceb9fe1a85ec53ULL; x^=x>>33; return x; }
static inline int64_t hfind(uint64_t x){ if(!Hcap) return -1; uint64_t h=hsh(x)&Hmask; while(HK[h]){ if(HK[h]==x) return (int64_t)h; h=(h+1)&Hmask; } return -1; }
static void hrehash(uint64_t nc){ uint64_t *ok=HK; uint8_t *ov=HV; uint64_t oc=Hcap; HK=calloc(nc,8); HV=calloc(nc,1); if(!HK||!HV){fprintf(stderr,"oom\n");exit(1);} Hcap=nc; Hmask=nc-1;
  for(uint64_t i=0;i<oc;i++) if(ok[i]){ uint64_t h=hsh(ok[i])&Hmask; while(HK[h]) h=(h+1)&Hmask; HK[h]=ok[i]; HV[h]=ov[i]; } free(ok); free(ov); }
static inline void hput(uint64_t x){ if((Hn+1)*10>Hcap*7) hrehash(Hcap?2*Hcap:(1ULL<<20)); uint64_t h=hsh(x)&Hmask; while(HK[h]){ if(HK[h]==x) return; h=(h+1)&Hmask; } HK[h]=x; HV[h]=0; Hn++; if(x>Hmaxkey) Hmaxkey=x; }
static inline int st(uint64_t y){ if(y<=CAP) return getS(idx_of(y)); int64_t h=hfind(y); return h<0?4:HV[h]; }
static uint64_t now;
static inline uint64_t maxpresent(void){ return Hmaxkey>now? Hmaxkey: now; }
static uint64_t *Q; static size_t qn=0,qcap=0;
static inline void qpush(uint64_t x){ if(qn==qcap){ qcap=qcap?2*qcap:(1<<20); Q=realloc(Q,qcap*8);} Q[qn++]=x; }
static inline void setres(uint64_t y,int v){ if(y<=CAP) setS(idx_of(y),v); else { int64_t h=hfind(y); HV[h]=(uint8_t)v; } qpush(y); }
static void propagate1(void){ while(qn){ uint64_t x=Q[--qn]; int v=getS(idx_of(x)); u128 y=(u128)x;
  for(int k=1;k<100;k++){ y<<=1; if(y>(u128)3*now+1) break; uint64_t yy=(uint64_t)y; uint64_t p=(yy%3==1)?(yy-1)/3:(yy+1)/3; if(p>now) break; if(p%3==0||p==x) continue;
    uint64_t pi=idx_of(p); int sp=getS(pi); if(sp>=2) continue; if(v==3) setres(p,2); else { if(sp==0) setS(pi,1); else setres(p,3); } } } }
static void propagate2(void){ while(qn){ uint64_t x=Q[--qn]; int v=st(x); uint64_t lim=maxpresent(); u128 y=(u128)x;
  for(int k=1;k<120;k++){ y<<=1; if(y>(u128)3*lim+1) break; uint64_t yy=(uint64_t)y; uint64_t p=(yy%3==1)?(yy-1)/3:(yy+1)/3; if(p>lim) break; if(p%3==0||p==x) continue;
    int sp=st(p); if(sp>=2) continue; if(v==3){ setres(p,2); continue; } uint64_t a=dchild(p), b=uchild(p); uint64_t o=(a==x)?b:a; int so=(o==1)?3:st(o); if(so==2) setres(p,3); } } }
static int status_start(uint64_t n){ if(n==1) return 2; if(n%3){ int s=st(n); return (s==2||s==3)?s:0; }
  uint64_t a=dchild(n), b=uchild(n); int sa=(a==1)?3:st(a), sb=st(b); if(sa==3||sb==3) return 2; if(sa==2&&sb==2) return 3; return 0; }
static uint64_t *HP; static size_t hpn=0,hpcap=0;
static void hpush(uint64_t x){ if(hpn==hpcap){ hpcap=hpcap?2*hpcap:(1<<20); HP=realloc(HP,hpcap*8);} size_t i=hpn++; while(i){ size_t p=(i-1)/2; if(HP[p]<=x) break; HP[i]=HP[p]; i=p; } HP[i]=x; }
static uint64_t hpop(void){ uint64_t top=HP[0], x=HP[--hpn]; size_t i=0; for(;;){ size_t l=2*i+1; if(l>=hpn) break; size_t c=(l+1<hpn && HP[l+1]<HP[l])? l+1:l; if(HP[c]>=x) break; HP[i]=HP[c]; i=c; } if(hpn) HP[i]=x; return top; }
static uint8_t *EXB; // explored bitmap for table (1 bit per index)
static inline int exget(uint64_t i){ return (EXB[i>>3]>>(i&7))&1; } static inline void exset(uint64_t i){ EXB[i>>3]|=(uint8_t)(1<<(i&7)); }
static uint64_t *XS; static size_t xsn=0,xscap=0;
static void xpush(uint64_t x){ if(xsn==xscap){ xscap=xscap?2*xscap:(1<<16); XS=realloc(XS,xscap*8);} XS[xsn++]=x; }
static void explore(uint64_t y0){ xpush(y0); while(xsn){ uint64_t y=XS[--xsn]; if(y==1) continue; int s=st(y);
  if(y<=CAP){ if(s>=2) continue; uint64_t i=idx_of(y); if(exget(i)) continue; exset(i); xpush(dchild(y)); xpush(uchild(y)); } else if(s==4) hpush(y); } }
static int relevant(uint64_t x){ uint64_t lim=maxpresent(); u128 y=(u128)x;
  for(int k=1;k<120;k++){ y<<=1; if(y>(u128)3*lim+1) break; uint64_t yy=(uint64_t)y; uint64_t p=(yy%3==1)?(yy-1)/3:(yy+1)/3; if(p>lim) break; if(p%3==0) continue;
    int s=st(p); if(s<2){ if(p>CAP || exget(idx_of(p))) return 1; } } return 0; }

// values of odd n < NB: V[(n-1)/2] = 2 (N) or 3 (P)
static uint8_t *V; static uint64_t NB;
static inline int val(uint64_t n){ return V[(n-1)>>1]; }
static int Ldown(uint64_t n){ int L=0; while(n!=1){ n=dchild(n); L++; } return L; }

int main(int argc,char**argv){
  if(argc<3){ fprintf(stderr,"usage: %s CAP BA\n",argv[0]); return 1; }
  CAP=strtoull(argv[1],0,10); if(CAP%2==0) CAP--;
  uint64_t BA=strtoull(argv[2],0,10); NB = BA + BA/2 + 4;
  if(NB > (CAP/3)*2){ fprintf(stderr,"need 1.5*BA <= 2CAP/3\n"); return 1; }
  MAXIDX=idx_of(CAP)+1; size_t words=(MAXIDX+31)/32+1; S=calloc(words,8); EXB=calloc(MAXIDX/8+2,1);
  time_t t0=time(0);
  setS(0,3); uint64_t f=1;
  for(uint64_t i=1;i<MAXIDX;i++){ uint64_t m=n_of(i); if(m>CAP) break; now=m; int sa=getS(idx_of(dchild(m))); if(sa==3) setres(m,2); else if(sa==2) setS(i,1); propagate1(); while(f<=m && status_start(f)) f+=2; }
  printf("# game_struct CAP=%llu BA=%llu NB=%llu; phase-1 frontier %llu (%lds)\n",(unsigned long long)CAP,(unsigned long long)BA,(unsigned long long)NB,(unsigned long long)f,(long)(time(0)-t0));
  // phase 2: all unresolved odd n < NB
  size_t nt=0; for(uint64_t n=f;n<NB;n+=2) if(!status_start(n)){ nt++; if(n%3) explore(n); else { explore(dchild(n)); explore(uchild(n)); } }
  uint64_t pops=0;
  for(uint64_t n=f;n<NB;n+=2){
    while(!status_start(n)){
      if(!hpn){ printf("HEAP EMPTY at %llu (trap)\n",(unsigned long long)n); return 2; }
      uint64_t x=hpop(); if(st(x)!=4) continue; if(!relevant(x)) continue; if(x>now) now=x;
      uint64_t a=dchild(x), b=uchild(x); int sa=(a==1)?3:st(a), sb=st(b); hput(x); pops++;
      if(sa==3||sb==3) setres(x,2); else if(sa==2&&sb==2) setres(x,3);
      if(qn) propagate2(); else { if(sa<2||sa==4) explore(a); if(sb<2||sb==4) explore(b); }
    }
  }
  printf("# phase 2: %zu targets, %llu insertions, max inserted %llu (%lds)\n",nt,(unsigned long long)pops,(unsigned long long)now,(long)(time(0)-t0));
  V=malloc(NB/2+2); for(uint64_t n=1;n<NB;n+=2) V[(n-1)>>1]=(uint8_t)status_start(n);
  free(S); free(EXB); free(HK); free(HV); free(HP);
  // ---------------- statistics on odd 3 <= n < BA ----------------
  // (S1) densities
  printf("\n## S1 P-density\n");
  { uint64_t cN=1,cP=0; for(uint64_t n=3;n<8;n+=2){ if(val(n)==3) cP++; else cN++; }
    for(int k=3;k<64;k++){ uint64_t lo=1ULL<<k, hi=lo<<1; if(hi>BA) break; uint64_t oN=0,oP=0;
      for(uint64_t n=lo+1;n<hi;n+=2){ if(val(n)==3) oP++; else oN++; }
      cN+=oN; cP+=oP; printf("octave [2^%d,2^%d): P-fraction %.5f   (cumulative odd n<2^%d, n=1 counted N: %.5f)\n",k,k+1,(double)oP/(oN+oP),k+1,(double)cP/(cN+cP)); } }
  // mantissa profile: bins of n/2^k in [1,2), 16 bins, per octave and pooled over octaves 20..top
  { int top=0; while((1ULL<<(top+1))<=BA) top++; top--; // last complete octave [2^top,2^{top+1})
    printf("mantissa profile: P-fraction in 16 bins of n/2^k in [1,2) (rows = octave k)\n");
    double pool[16]={0}, poolt[16]={0};
    for(int k=16;k<=top;k++){ uint64_t lo=1ULL<<k; uint64_t c[16]={0},p[16]={0};
      for(uint64_t n=lo+1;n<(lo<<1);n+=2){ int b=(int)(((n-lo)*16)>>k); c[b]++; if(val(n)==3) p[b]++; }
      printf(" k=%2d:",k); for(int b=0;b<16;b++){ printf(" %.4f",(double)p[b]/c[b]); if(k>=20){ pool[b]+=p[b]; poolt[b]+=c[b]; } } printf("\n"); }
    printf(" pooled k>=20:"); for(int b=0;b<16;b++) printf(" %.4f",pool[b]/poolt[b]); printf("\n"); }
  // (S2) residues mod 2^k : symmetry r <-> -r
  printf("\n## S2 residues (odd 3<=n<BA)\n");
  for(int k=3;k<=12;k++){ uint64_t M=1ULL<<k; uint64_t *c=calloc(M,8),*p=calloc(M,8);
    for(uint64_t n=3;n<BA;n+=2){ c[n&(M-1)]++; if(val(n)==3) p[n&(M-1)]++; }
    double mx=0,mxs=0,lo=1,hi=0; for(uint64_t r=1;r<M;r+=2){ double a=(double)p[r]/c[r], b=(double)p[M-r]/c[M-r]; double se=sqrt(a*(1-a)/c[r]+b*(1-b)/c[M-r]); double d=fabs(a-b); if(d>mx) mx=d; if(d/se>mxs) mxs=d/se; if(a<lo) lo=a; if(a>hi) hi=a; }
    printf("mod 2^%-2d: P-fraction range [%.4f, %.4f]; max |p(r)-p(-r)| = %.5f (max z-score %.2f over %llu pairs)\n",k,lo,hi,mx,mxs,(unsigned long long)(M/4));
    if(k==5){ printf("   mod 32 table r:p(r) :"); for(uint64_t r=1;r<M;r+=2) printf(" %llu:%.4f",(unsigned long long)r,(double)p[r]/c[r]); printf("\n"); }
    free(c); free(p); }
  { for(int j=1;j<=4;j++){ uint64_t M=1; for(int q=0;q<j;q++) M*=3; uint64_t *c=calloc(M,8),*p=calloc(M,8);
      for(uint64_t n=3;n<BA;n+=2){ c[n%M]++; if(val(n)==3) p[n%M]++; }
      double lo=1,hi=0; for(uint64_t r=0;r<M;r++) if(c[r]){ double a=(double)p[r]/c[r]; if(a<lo) lo=a; if(a>hi) hi=a; }
      printf("mod 3^%d: P-fraction range [%.4f, %.4f]",j,lo,hi); if(j<=2){ printf(" :"); for(uint64_t r=0;r<M;r++) if(c[r]) printf(" %llu:%.4f",(unsigned long long)r,(double)p[r]/c[r]); } printf("\n"); free(c); free(p); } }
  // descent exponent
  { uint64_t c[64]={0},p[64]={0}; for(uint64_t n=3;n<BA;n+=2){ int j=dexp(n); c[j]++; if(val(n)==3) p[j]++; }
    printf("descent exponent j=v2(3n+-1) of the descending move: P-fraction:"); for(int j=2;j<64;j++) if(c[j]>1000) printf(" j=%d:%.4f(%.3f%%)",j,(double)p[j]/c[j],100.0*c[j]/(BA/2)); printf("\n"); }
  // (S3) move structure
  printf("\n## S3 winning moves (odd 3<=n<BA)\n");
  { uint64_t N=0,P=0, dOnly=0,uOnly=0,both=0, imm=0; uint64_t firstU[20]; int nf=0;
    uint64_t uOnlyByPlus=0; // forced ascent on the plus sheet
    for(uint64_t n=3;n<BA;n+=2){ int v=val(n); if(v==3){P++; continue;} N++;
      uint64_t a=dchild(n), b=uchild(n); if(a==1){ imm++; dOnly++; continue; }
      int va=val(a), vb=val(b);
      if(va==3 && vb==3) both++; else if(va==3) dOnly++; else if(vb==3){ uOnly++; if(nf<20) firstU[nf++]=n; if(b==Uplus(n)) uOnlyByPlus++; }
      else { printf("INCONSISTENT at %llu\n",(unsigned long long)n); } }
    printf("N=%llu P=%llu; among N: only-descent wins %.4f (incl. immediate wins %llu), only-ascent wins %.4f, both win %.4f\n",(unsigned long long)N,(unsigned long long)P,(double)dOnly/N,(unsigned long long)imm,(double)uOnly/N,(double)both/N);
    printf("forced-ascent positions (only the ascending move wins): %llu = %.4f of all odd n; on plus sheet (n=3 mod 4): %.4f of them\n",(unsigned long long)uOnly,(double)uOnly/(N+P),(double)uOnlyByPlus/uOnly);
    printf("first forced-ascent positions:"); for(int i=0;i<nf;i++) printf(" %llu",(unsigned long long)firstU[i]); printf("\n"); }
  // (S4) down-only game
  printf("\n## S4 down-only game (value = parity of the descending chain length L(n) to 1)\n");
  { uint64_t agree=0,tot=0, cm[2][2]={{0,0},{0,0}}; uint64_t firstdis[20]; int nf=0; int Lmax=0; uint64_t Lsum=0;
    uint64_t byL_c[200]={0}, byL_p[200]={0};
    for(uint64_t n=3;n<BA;n+=2){ int L=Ldown(n); if(L>Lmax) Lmax=L; Lsum+=L; int vd = (L&1)? 2:3; int v=val(n); tot++; if(v==vd) agree++; else if(nf<20) firstdis[nf++]=n; cm[v==3][vd==3]++; if(L<200){ byL_c[L]++; if(v==3) byL_p[L]++; } }
    printf("agreement %.4f; confusion (true,down): NN=%llu NP=%llu PN=%llu PP=%llu; mean L=%.2f max L=%d\n",(double)agree/tot,(unsigned long long)cm[0][0],(unsigned long long)cm[0][1],(unsigned long long)cm[1][0],(unsigned long long)cm[1][1],(double)Lsum/tot,Lmax);
    printf("first disagreements:"); for(int i=0;i<nf;i++) printf(" %llu",(unsigned long long)firstdis[i]); printf("\n");
    printf("P-fraction by L (L even => down-only says P):"); for(int L=1;L<200;L++) if(byL_c[L]>20000) printf(" %d:%.3f",L,(double)byL_p[L]/byL_c[L]); printf("\n"); }
  // one-step corrections: v(n) vs rule "N iff v(d(n))=P" (= descent-only one step), vs "N iff d(n) P or u(n) P" (true)
  // (S7) sheets: the winning move's sheet for N positions, by n mod 4
  printf("\n## S7 sheets (n=1 mod 4: descending move on plus sheet; n=3 mod 4: on minus sheet)\n");
  { uint64_t c[4][4]={{0}}; // [n mod 4][pattern]: pattern bit0 = plus-child P, bit1 = minus-child P
    for(uint64_t n=3;n<BA;n+=2){ uint64_t a=Uplus(n), b=Uminus(n); int pa=(a==1)||(val(a)==3), pb=(b==1)||(val(b)==3); c[n&3][pa|(pb<<1)]++; }
    for(int r=1;r<4;r+=2){ uint64_t t=c[r][0]+c[r][1]+c[r][2]+c[r][3]; printf("n=%d mod 4: P=%.4f  only plus-child P=%.4f  only minus-child P=%.4f  both=%.4f\n",r,(double)c[r][0]/t,(double)c[r][1]/t,(double)c[r][2]/t,(double)c[r][3]/t); } }
  // (S8) reflection tests: v(n) vs v(m), m = 3*2^k - n (m = -n mod 2^k), n in [2^k, 2^{k+1})
  printf("\n## S8 reflection tests\n");
  { int top=0; while((1ULL<<(top+1))<=BA) top++; top--;
    for(int k=top-2;k<=top;k++){ uint64_t lo=1ULL<<k; uint64_t agree=0,tot=0,pP=0; for(uint64_t n=lo+1;n<2*lo;n+=2){ uint64_t m=3*lo-n; int a=val(n), b=val(m); tot++; if(a==b) agree++; if(a==3) pP++; }
      double p=(double)pP/tot; printf("k=%d: agreement of v(n) and v(3*2^k-n) = %.4f (independence baseline %.4f)\n",k,(double)agree/tot,p*p+(1-p)*(1-p)); }
    // and v(n) vs v(-n) statistics conditioned on local data: compare P-fraction given (n mod 2^10) vs (-n mod 2^10) already in S2
  }
  // (S5) ascending rays and the exit index
  printf("\n## S5 ascending rays: k(x) = min{i>=0 : d(u^i(x)) is P or 1}; rule: x is N iff k(x) even\n");
  { uint64_t lim=BA>>12; uint64_t hist[64]={0}, cens=0, viol=0, tot=0; uint64_t tr[2][2]={{0,0},{0,0}}; // side-value transitions along rays
    for(uint64_t x=3;x<lim;x+=2){ uint64_t y=x; int k=-1; for(int i=0;i<60;i++){ uint64_t a=dchild(y); if(a>=NB||y>=NB){ break; } int sa=(a==1)?3:val(a); if(sa==3){ k=i; break; } y=uchild(y); }
      tot++; if(k<0){ cens++; continue; } hist[k<63?k:63]++; int pred=(k%2==0)?2:3; if(pred!=val(x)) viol++; }
    for(uint64_t x=3;x<(BA>>1);x+=2){ uint64_t a=dchild(x), y=uchild(x), b=dchild(y); if(y>=NB) continue; int s0=(a==1)||val(a)==3, s1=(b==1)||val(b)==3; tr[s0][s1]++; }
    printf("x < %llu: exit-index rule violations %llu of %llu (censored %llu)\n",(unsigned long long)lim,(unsigned long long)viol,(unsigned long long)tot,(unsigned long long)cens);
    printf("distribution of k:"); for(int k=0;k<64;k++) if(hist[k]) printf(" %d:%.5f",k,(double)hist[k]/(tot-cens)); printf("\n");
    double q=(double)hist[0]/(tot-cens); printf("P(k=0)=q=%.4f; geometric model would give P(k even)=1/(2-q)=%.4f vs observed %.4f\n",q,1.0/(2-q), ({ double e=0; for(int k=0;k<64;k+=2) e+=hist[k]; e/(tot-cens); }));
    double t00=tr[0][0],t01=tr[0][1],t10=tr[1][0],t11=tr[1][1];
    printf("side values along a ray step (s_i -> s_{i+1}, s = [d-child is P]): P(s1=P|s0=N)=%.4f  P(s1=P|s0=P)=%.4f\n",t01/(t00+t01),t11/(t10+t11)); }
  // (S9) phase predictor on the last complete octave
  printf("\n## S9 phase predictor: value vs theta = frac(log2 n) (majority vote per theta-bin, last octave)\n");
  { int top=0; while((1ULL<<(top+1))<=BA) top++; top--; uint64_t lo=1ULL<<top;
    for(int b=6;b<=14;b+=2){ uint64_t nb=1ULL<<b; uint64_t *c=calloc(nb,8),*p=calloc(nb,8);
      for(uint64_t n=lo+1;n<2*lo;n+=2){ double th=log2((double)n)-top; uint64_t i=(uint64_t)(th*nb); if(i>=nb) i=nb-1; c[i]++; if(val(n)==3) p[i]++; }
      uint64_t ok=0,t=0; for(uint64_t i=0;i<nb;i++){ ok += (2*p[i]>=c[i])? p[i] : c[i]-p[i]; t+=c[i]; }
      // combined with n mod 16
      uint64_t *c2=calloc(nb*16,8),*p2=calloc(nb*16,8);
      for(uint64_t n=lo+1;n<2*lo;n+=2){ double th=log2((double)n)-top; uint64_t i=(uint64_t)(th*nb); if(i>=nb) i=nb-1; uint64_t j=i*16+(n&15); c2[j]++; if(val(n)==3) p2[j]++; }
      uint64_t ok2=0; for(uint64_t j=0;j<nb*16;j++) ok2 += (2*p2[j]>=c2[j])? p2[j] : c2[j]-p2[j];
      printf("octave 2^%d, %llu theta-bins: accuracy %.4f ; with n mod 16 as well: %.4f\n",top,(unsigned long long)nb,(double)ok/t,(double)ok2/t);
      free(c);free(p);free(c2);free(p2); } }
  // (S10) phase twins: the two children have the same phase up to O(1/n)
  { uint64_t same=0,t=0; for(uint64_t n=3;n<BA;n+=2){ uint64_t a=dchild(n), b=uchild(n); if(a==1) continue; t++; if(val(a)==val(b)) same++; }
    printf("\n## S10 the two children d(n),u(n) (same phase theta+log2 3 mod 1, up to O(1/n)) have equal values in %.4f of cases\n",(double)same/t); }
  // (S11) pointwise tests: translation partner n+2^m (same residue mod 2^m, same phase up to 2^m/n)
  //       vs mirror partner n' = 2^m*ceil(2n/2^m) - n (n' = -n mod 2^m, n <= n' < n+2^m)
  { int top=0; while((1ULL<<(top+1))<=BA) top++; top--; uint64_t lo=1ULL<<top;
    printf("\n## S11 pointwise partners in octave 2^%d: mirror partner n' = n+delta, delta = (-2n mod 2^m) (so n' = -n mod 2^m);\n##     control partner n+delta'' with delta'' the mirror distance of an unrelated point (same distance law, no residue relation)\n",top);
    for(int m=4;m<=16;m+=2){ uint64_t M=1ULL<<m; uint64_t ac=0,am=0,t=0; for(uint64_t n=lo+1;n+2*M<2*lo;n+=2){ uint64_t dm=((2*n+M-1)/M)*M-2*n; uint64_t n2=n+(M<<3)+2*((n*2654435761ULL)%M); uint64_t dc=((2*n2+M-1)/M)*M-2*n2; if(n+dm>=NB||n+dc>=NB) continue; t++; if(val(n)==val(n+dm)) am++; if(val(n)==val(n+dc)) ac++; }
      printf(" m=%2d: mirror %.4f  control %.4f\n",m,(double)am/t,(double)ac/t); } }
  printf("\n# done (%lds)\n",(long)(time(0)-t0));
  return 0;
}
