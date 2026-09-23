// collatz_procgen_20260922_endgame_psi.c
// Canonical-escape (Psi) orbits in the backward E-game (Q2), exact 128-bit integer arithmetic.
//
// Psi(m): if m = 1 (mod 3):  k = v_3(m-1),  Psi(m) = 2^{K*(k-1)}   (m-1)/3^k      (1-escape at precision k)
//         if m = 2 (mod 3):  k = v_3(2m-1), Psi(m) = 2^{K*(k-1)} (2m-1)/3^k      (1/2-transfer at precision k)
// K*(0)=0, K*(1)=2, K*(s)=floor((s+1) log2 3) (s>=2; loops through 1 exist for s<=6000, FINITE-EXACT).
// Every Psi step is an explicit legal reverse path of k moves (loop + final move 0); its 3-adic multiplier is
// rho = 2^{K*(k-1)+[half]}/3^k and Psi(m) < rho*m.  Only half-steps with k>=3 have rho > 1 ("L1 links").
// Budget (PROVED): rho <= exp(c* k) for every step, c* = ln(128/81)/4, so m_t/m_0 <= exp(c* D_L1(t)).
//
// modes:
//   scan N [stride_start stride]  : every m = 14 mod 27, 14 <= m <= N: orbit until m_t < m_0
//   sample FILE M NPER SEED        : hostile points "p e" from FILE; random m=(p+3^j W)/2^e <= M near each
//   dfs NT M                       : exhaustive DFS over Psi-alive classes mod 3^NT (m_0 = 14 mod 27 not
//                                    required), least representatives <= M followed exactly
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128;
static int KS[512];
static int FL[1024];   // FL[D] = floor(D log2 3): 2^A < 3^D  <=>  A <= FL[D]  (D >= 1)
static u128 P3[80];
static const double CSTAR = 0.11439540225190346; // ln(128/81)/4
static const double LN2 = 0.6931471805599453, LN3 = 1.0986122886681098;

static void init(void){
  P3[0]=1; for(int i=1;i<80;i++) P3[i]=P3[i-1]*3;
  KS[0]=0; KS[1]=2;
  for(int s=2;s<512;s++){ // floor((s+1) log2 3) exactly: largest K with 2^K < 3^(s+1)
    // use long double estimate then correct with exact check for s+1 <= 79
    long double v=(s+1)*1.5849625007211561814537389439478L; int K=(int)floorl(v);
    if(s+1<80){ while(K>0 && ((u128)1<<(K<127?K:126)) > P3[s+1]) K--; while(K+1<127 && ((u128)1<<(K+1)) < P3[s+1]) K++; }
    KS[s]=K;
  }
  for(int D=1;D<1024;D++){ long double v=D*1.5849625007211561814537389439478L; int K=(int)floorl(v);
    if(D<80){ while(((u128)1<<K) > P3[D]) K--; while(K+1<127 && ((u128)1<<(K+1)) < P3[D]) K++; }
    FL[D]=K; }
}
static void pr128(FILE*f,u128 x){ char b[64]; int n=0; if(!x){fputc('0',f);return;} while(x){b[n++]='0'+(int)(x%10);x/=10;} while(n)fputc(b[--n],f); }

typedef struct { int half, k; } Step;
// one Psi step; returns 0 on overflow (value would exceed 2^126) or at the base point 1.
static inline int psi(u128 m, u128 *y, Step *st){
  u128 t; int half;
  if(m%3==1){ t=m-1; half=0; } else { t=2*m-1; half=1; }
  if(t==0) return 0;
  int k=0; while(t%3==0){ t/=3; k++; }
  int a=KS[k-1];
  if(a>=126 || (t>>(126-a))!=0) return 0;
  *y = t<<a; st->half=half; st->k=k; return 1;
}

typedef struct {
  int steps, L1, DL1, Dtot, runmax_links, runmax_digits; double runmax_logprice;
  double maxlog;      // max_t ln(m_t/m_0)
  double maxlog_budget_excess; // max_t [ln(m_t/m_0) - c* DL1(t)]  (must be <= 0)
  int descended; int maxlog_t;
  int A, D; int mult_at_descent;   // exact multiplier 2^A/3^D of the descending prefix; 1 if 2^A < 3^D
  int first_mult_desc;             // first t (steps) with 2^A < 3^D (0 if none within the orbit)
} Orbit;

static int orbit(u128 m0, int cap, Orbit *o){
  memset(o,0,sizeof *o); o->maxlog_budget_excess=-1e300;
  u128 m=m0; int run_links=0, run_digits=0; double run_lp=0; double lm0=0;
  // ln(m/m0) tracked exactly-ish via long double of ratio of u128
  for(int t=0;t<cap;t++){
    u128 y; Step st;
    if(!psi(m,&y,&st)) return 0;
    o->steps++; o->Dtot+=st.k; o->A += (st.half? KS[st.k-1]+1 : KS[st.k-1]); o->D += st.k;
    if(!o->first_mult_desc && o->A <= FL[o->D]) o->first_mult_desc=o->steps;
    double lr = (st.half? (KS[st.k-1]+1):KS[st.k-1])*LN2 - st.k*LN3;
    if(st.half && st.k>=3){ o->L1++; o->DL1+=st.k; run_links++; run_digits+=st.k; run_lp+=lr;
      if(run_links>o->runmax_links || (run_links==o->runmax_links && run_lp>o->runmax_logprice)){ o->runmax_links=run_links; o->runmax_digits=run_digits; o->runmax_logprice=run_lp; } }
    else { run_links=0; run_digits=0; run_lp=0; }
    m=y;
    long double ratio=(long double)m/(long double)m0; double lg=(double)logl(ratio);
    if(lg>o->maxlog){ o->maxlog=lg; o->maxlog_t=o->steps; }
    double ex=lg - CSTAR*o->DL1; if(ex>o->maxlog_budget_excess) o->maxlog_budget_excess=ex;
    if(m<m0){ o->descended=1; o->mult_at_descent = (o->first_mult_desc!=0); return 1; }
  }
  return 1;
}

static uint64_t rng_s[2];
static inline uint64_t rotl(const uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static uint64_t rng(void){ uint64_t s0=rng_s[0],s1=rng_s[1],r=s0+s1; s1^=s0; rng_s[0]=rotl(s0,55)^s1^(s1<<14); rng_s[1]=rotl(s1,36); return r; }
static u128 rng128(void){ return ((u128)rng()<<64)|rng(); }

// ------------------------------------------------------------------------------------------------ scan
static long long nonmult_scan=0;
static void do_scan(u128 N, u128 start, u128 stride){
  long long cnt=0, nodesc=0, overflow=0; long long hist_steps[64]={0}, hist_L1[64]={0};
  Orbit best_steps={0}, best_L1={0}, best_run={0}, best_max={0}, best_DL1={0}; u128 a_steps=0,a_L1=0,a_run=0,a_max=0,a_DL1=0;
  double worst_budget=-1e300; u128 a_budget=0; double sum_steps=0;
  for(u128 m=start; m<=N; m+=stride){
    Orbit o; int ok=orbit(m,100000,&o); cnt++;
    if(!ok){ overflow++; continue; }
    if(!o.descended){ nodesc++; fprintf(stderr,"no descent within cap: "); pr128(stderr,m); fprintf(stderr,"\n"); continue; }
    sum_steps+=o.steps; if(!o.mult_at_descent) nonmult_scan++;
    hist_steps[o.steps<63?o.steps:63]++; hist_L1[o.L1<63?o.L1:63]++;
    if(o.steps>best_steps.steps){best_steps=o;a_steps=m;}
    if(o.L1>best_L1.L1){best_L1=o;a_L1=m;}
    if(o.DL1>best_DL1.DL1){best_DL1=o;a_DL1=m;}
    if(o.runmax_links>best_run.runmax_links || (o.runmax_links==best_run.runmax_links && o.runmax_logprice>best_run.runmax_logprice)){best_run=o;a_run=m;}
    if(o.maxlog>best_max.maxlog){best_max=o;a_max=m;}
    if(o.maxlog_budget_excess>worst_budget){worst_budget=o.maxlog_budget_excess;a_budget=m;}
  }
  printf("scan m = "); pr128(stdout,start); printf(" (step "); pr128(stdout,stride); printf(") up to "); pr128(stdout,N);
  printf(": %lld orbits, not descended within cap %lld, overflow %lld, mean Psi-steps to descent %.4f\n",cnt,nodesc,overflow,sum_steps/(cnt-nodesc-overflow));
  printf("  histogram of Psi-steps to descent:"); for(int i=0;i<64;i++) if(hist_steps[i]) printf(" %d:%lld",i,hist_steps[i]); printf("\n");
  printf("  histogram of L1 links (half, k>=3) before descent:"); for(int i=0;i<64;i++) if(hist_L1[i]) printf(" %d:%lld",i,hist_L1[i]); printf("\n");
  #define REP(lbl,o,a) do{ printf("  %-28s m=",lbl); pr128(stdout,a); printf("  steps=%d L1=%d D_L1=%d D_tot=%d longest L1-run=%d links/%d digits (price %.4f) max excursion %.4f (ln %.4f) vs budget exp(c* D_L1)=%.4f\n", \
      (o).steps,(o).L1,(o).DL1,(o).Dtot,(o).runmax_links,(o).runmax_digits,exp((o).runmax_logprice),exp((o).maxlog),(o).maxlog,exp(CSTAR*(o).DL1)); }while(0)
  REP("most Psi-steps:",best_steps,a_steps);
  REP("most L1 links:",best_L1,a_L1);
  REP("most L1 digits:",best_DL1,a_DL1);
  REP("longest consecutive L1 run:",best_run,a_run);
  REP("largest excursion:",best_max,a_max);
  printf("  descents that are not multiplicative: %lld\n",nonmult_scan);
  printf("  budget check: max over orbits and times of ln(m_t/m_0) - c* D_L1(t) = %.6f (<= 0 required) at m=",worst_budget); pr128(stdout,a_budget); printf("\n");
  fflush(stdout);
}

// ------------------------------------------------------------------------------------------------ sample
static void do_sample(const char*fn, u128 M, int nper, uint64_t seed){
  rng_s[0]=seed^0x9E3779B97F4A7C15ULL; rng_s[1]=seed*0xBF58476D1CE4E5B9ULL+1;
  FILE*f=fopen(fn,"r"); if(!f){perror(fn);exit(1);}
  char pbuf[128]; int e;
  printf("sample: m = (p + 3^j W)/2^e <= "); pr128(stdout,M); printf(", %d random W per (point, j)\n",nper);
  double glob_max=0; u128 glob_arg=0; int glob_L1=0; u128 glob_L1arg=0; double glob_budget=-1e300; long long glob_cnt=0, glob_nodesc=0;
  while(fscanf(f,"%127s %d",pbuf,&e)==2){
    u128 p=0; for(char*c=pbuf;*c;c++) p=p*10+(*c-'0');
    u128 twoe=(u128)1<<e;
    double best=0; u128 barg=0; int bestL1=0, bestrun=0; double bestrunp=0; long long cnt=0, nodesc=0; int jmax=0;
    for(int j=1;j<78 && P3[j] < M*twoe; j++){
      // need W = -p * 3^{-j} mod 2^e, W >= 1, 3 !| W, m=(p+3^j W)/2^e <= M
      u128 inv3=1; { // inverse of 3^j mod 2^e (e<=100): Newton iteration mod 2^128 then reduce
        u128 a=P3[j], x=1; for(int it=0;it<7;it++) x=x*(2-a*x); inv3 = (e>=128)? x : (x & (twoe-1)); }
      u128 W0 = e? ((twoe - (p % twoe)) % twoe) * inv3 % twoe : 0; // careful: product may overflow for e>63
      if(e>63){ // recompute with 128-bit modular multiply by bits
        u128 a=(twoe - (p%twoe))%twoe, b=inv3, r=0; while(b){ if(b&1) r=(r+a)&(twoe-1); a=(a<<1)&(twoe-1); b>>=1; } W0=r; }
      u128 Wmax = (M*twoe - p)/P3[j];
      if(Wmax < W0+1) continue;
      u128 span = (Wmax - W0)/twoe; if(span==0) continue;
      jmax=j;
      for(int s=0;s<nper;s++){
        u128 W = W0 + twoe*(rng128()%(span+1)); if(W%3==0) W+=twoe; if(W%3==0) W+=twoe; if(W>Wmax||W==0) continue;
        u128 num=p+P3[j]*W; if(num%twoe) { fprintf(stderr,"integrality bug\n"); exit(2);} u128 m=num/twoe; if(m<2) continue;
        Orbit o; if(!orbit(m,100000,&o)) continue; cnt++;
        if(!o.descended){ nodesc++; continue; }
        if(o.maxlog>best){best=o.maxlog;barg=m;}
        if(o.L1>bestL1) bestL1=o.L1;
        if(o.runmax_links>bestrun){bestrun=o.runmax_links; bestrunp=o.runmax_logprice;}
        if(o.maxlog>glob_max){glob_max=o.maxlog;glob_arg=m;}
        if(o.L1>glob_L1){glob_L1=o.L1;glob_L1arg=m;}
        if(o.maxlog_budget_excess>glob_budget) glob_budget=o.maxlog_budget_excess;
      }
    }
    glob_cnt+=cnt; glob_nodesc+=nodesc;
    printf("  %s/2^%d: j<=%d, %lld orbits (no descent %lld): max L1 links %d, longest L1 run %d (price %.4f), max excursion %.4f at m=",pbuf,e,jmax,cnt,nodesc,bestL1,bestrun,exp(bestrunp),exp(best));
    pr128(stdout,barg); printf("\n"); fflush(stdout);
  }
  printf("  all points: %lld orbits, %lld without descent; max excursion %.4f at m=",glob_cnt,glob_nodesc,exp(glob_max)); pr128(stdout,glob_arg);
  printf("; max L1 links %d at m=",glob_L1); pr128(stdout,glob_L1arg); printf("; budget excess max %.6f (<= 0 required)\n",glob_budget);
  fclose(f);
}

// ------------------------------------------------------------------------------------------------ dfs
// DFS over the least representatives r (mod 3^n) of Psi-alive classes.  State: r, n, and the image
// x = (2^A r - B)/3^D of r under the completed steps (exact), A, D.  Adding digit d at position n adds
// 2^A d 3^(n-D) to x.
static int NT; static u128 MM; static long long dfs_nodes=0, leaves=0, leaves_le=0, alive_leaf_nodesc=0;
static int top_D=0; static double top_max=0; static u128 top_Darg=0, top_maxarg=0; static int hist_surv[400];
static int top_run=0, top_rundig=0; static double top_runp=0; static u128 top_runarg=0; static int top_L1=0; static u128 top_L1arg=0;
static long long nonmult=0; static int top_DL1=0; static u128 top_DL1arg=0;
static void follow_leaf(u128 r){
  leaves++;
  if(r>MM || r<2) return;
  leaves_le++;
  Orbit o; if(!orbit(r,100000,&o)){ fprintf(stderr,"overflow at leaf\n"); return; }
  if(!o.descended){ alive_leaf_nodesc++; fprintf(stderr,"leaf without descent: "); pr128(stderr,r); fprintf(stderr,"\n"); return; }
  // digits consumed up to descent
  int D=o.Dtot; hist_surv[D<399?D:399]++;
  if(D>top_D){top_D=D; top_Darg=r;}
  if(o.maxlog>top_max){top_max=o.maxlog; top_maxarg=r;}
  if(o.runmax_links>top_run || (o.runmax_links==top_run && o.runmax_logprice>top_runp)){top_run=o.runmax_links; top_rundig=o.runmax_digits; top_runp=o.runmax_logprice; top_runarg=r;}
  if(o.L1>top_L1){top_L1=o.L1; top_L1arg=r;}
  if(o.DL1>top_DL1){top_DL1=o.DL1; top_DL1arg=r;}
  if(!o.mult_at_descent){ nonmult++; fprintf(stderr,"non-multiplicative first descent at leaf m="); pr128(stderr,r);
    // continue the Psi-orbit until the prefix multiplier drops below 1 (multiplicative certificate)
    u128 x=r; int A=0,D=0,steps=0,ok=0; while(steps<100000){ u128 y; Step st; if(!psi(x,&y,&st)) break; steps++; A+= (st.half? KS[st.k-1]+1:KS[st.k-1]); D+=st.k;
      if(A <= FL[D]){ ok=1; break; } x=y; }
    fprintf(stderr,"  multiplicative Psi-certificate after %d steps: %s (A=%d, D=%d)\n",steps,ok?"yes":"NO",A,D); }
}
static void dfs(u128 r, int n, u128 x, int A, int D){
  // known: r mod 3^n (r < 3^n), x = image mod 3^(n-D) (exact for r), cumulative 2^A/3^D >= 1
  dfs_nodes++;
  // apply every step that is now determined
  for(;;){
    int avail = n - D;                   // digits of x known
    if(avail<=0) break;
    u128 t; int half;
    if(x%3==1){ t=x-1; half=0; } else { t=2*x-1; half=1; }
    // valuation known iff t != 0 mod 3^avail
    int k=0; u128 tt=t; while(k<avail && tt%3==0){ tt/=3; k++; }
    if(k==avail) break;                  // undetermined
    // step determined: k digits consumed, needs avail >= k+1 (true since k < avail)
    int a=KS[k-1]+half;
    A+=a; D+=k;
    // multiplicative test 2^A < 3^D: dead
    if(A <= FL[D]) return;
    x = tt<<KS[k-1];
    if(half){ /* tt = (2x-1)/3^k already includes the factor 2 */ }
  }
  if(n==NT){ follow_leaf(r); return; }
  for(int d=0; d<3; d++){
    u128 r2 = r + (u128)d*P3[n];
    if(n==0 && d==0) continue;           // unit classes only
    if(r2 > MM && d>0) { /* least rep already > MM: no rep <= MM in this class */ continue; }
    // x update: x2 = x + 2^A d 3^(n-D) ... only when D <= n (always true here)
    u128 add = ((u128)d * P3[n-D]) << A;   // may overflow for large A: guard
    if(A>100){ fprintf(stderr,"A too large\n"); exit(3); }
    dfs(r2, n+1, x + add, A, D);
  }
}
static void do_dfs(int nt, u128 M){
  NT=nt; MM=M;
  // start: r in {1,2} (n=1), x = r, A=D=0.  (units mod 3)
  dfs(1,1,1,0,0); dfs(2,1,2,0,0);
  printf("dfs over Psi-alive classes to 3^%d, representatives <= ",nt); pr128(stdout,M);
  printf(": nodes %lld, alive leaves %lld, leaves with least rep <= M: %lld, of which not descended within cap: %lld\n",dfs_nodes,leaves,leaves_le,alive_leaf_nodesc);
  printf("  largest digits-to-descent D among them: %d at m=",top_D); pr128(stdout,top_Darg);
  printf("\n  largest excursion among them: %.4f at m=",exp(top_max)); pr128(stdout,top_maxarg);
  printf("\n  longest consecutive L1 run: %d links / %d digits (price %.4f) at m=",top_run,top_rundig,exp(top_runp)); pr128(stdout,top_runarg);
  printf("\n  most L1 links: %d at m=",top_L1); pr128(stdout,top_L1arg);
  printf("\n  most L1 digits: %d at m=",top_DL1); pr128(stdout,top_DL1arg);
  printf("\n  descents that are not multiplicative (2^A >= 3^D at the first m_t < m_0): %lld",nonmult);
  printf("\n  histogram of D (digits consumed until descent):");
  for(int i=0;i<400;i++) if(hist_surv[i]) printf(" %d:%d",i,hist_surv[i]); printf("\n");
}

int main(int argc,char**argv){
  init();
  if(argc<2){ fprintf(stderr,"usage\n"); return 1; }
  if(!strcmp(argv[1],"scan")){
    u128 N=strtoull(argv[2],0,10); u128 st=14, stride=27;
    if(argc>4){ st=strtoull(argv[3],0,10); stride=strtoull(argv[4],0,10); }
    do_scan(N,st,stride);
  } else if(!strcmp(argv[1],"sample")){
    u128 M=strtoull(argv[3],0,10); do_sample(argv[2],M,atoi(argv[4]),strtoull(argv[5],0,10));
  } else if(!strcmp(argv[1],"dfs")){
    do_dfs(atoi(argv[2]), strtoull(argv[3],0,10));
  } else if(!strcmp(argv[1],"orbit")){
    u128 m=0; for(char*c=argv[2];*c;c++) m=m*10+(*c-'0');
    u128 x=m; printf("orbit of "); pr128(stdout,m); printf(":\n");
    for(int t=0;t<200;t++){ u128 y; Step st; if(!psi(x,&y,&st)) break; printf("  %s k=%d -> ",st.half?"1/2":"1  ",st.k); pr128(stdout,y);
      printf("   ratio %.6f\n",(double)((long double)y/(long double)m)); x=y; if(y<m) break; }
  }
  return 0;
}
