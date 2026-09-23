// collatz_procgen_20260922_mirror_loops_dp.c
// Q1 mirror of loops_height.c / loops_dp.c: E-loops through -1 in the NEGATIVE E-graph.
//   Forward E-moves on negative integers: v -> 3v+1 (always), v -> v/2 (v even).
//   A loop through -1 with a multiplications and K halvings satisfies B = 3^a - 2^K,
//   B = sum_{i<a} 3^(a-1-i) 2^(K_i)  (K_i = halvings before the i-th multiplication).
//   Lower bound (mirror of Prop 1.2): a >= 2 implies 3^a > 2^(K+1), i.e. a >= a0(K) := ceil((K+1) log_3 2).
// DP layered by halvings.  Layer b holds the states w reached right after the b-th halving
// (layer 0 = {-1}); A_b(w) = min #multiplications, H_b(w) = among those, min height
// (height = max |multiplication point|, i.e. max |v| over the v at which 3v+1 is applied).
// Transition from state u: t >= 0 multiplications (t has the parity of u) then one halving.
// Cap: every multiplication point |v| <= HCAP (so states |w| <= (3 HCAP - 1)/2).
// Pruning (exact for loops with a <= a0(K)+SLACK, see the note): a state w at layer L with
//   A_L(w) >= (L+1) log_3 2 + log_3|w| + 1 + SLACK  cannot lie on such a loop (continuation needs
//   3^{a''}|w| >= 2^{K-L}).  SLACK < 0 disables pruning.
// Lexicographic (min a, then min height) is exact at every K with amin(K) = a0(K) (an a0-loop is
// a-minimal at every intermediate state, else a loop with a < a0(K) would exist).
// Usage: mirror_loops_dp KMAX HCAP SLACK [RECON_KMAX [RECON_KMIN]]
//   prints per K: a0, amin (within cap), ratio 3^amin/2^K, eta(K+1), H(K), #live states.
//   RECON_KMAX>0: keeps every layer's A-array (2 bytes/state/layer) and prints one a-minimal loop
//   for each RECON_KMIN <= K <= RECON_KMAX as a forward word (M^x H^y ...).  With HCAP = H(K) the
//   printed loop has minimal height.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#define NONE 0xFFFF
static long double TH;   // log_3 2
static int a0_of(int K){ // least a with 3^a > 2^(K+1)
  long double v=(long double)(K+1)*TH; int a=(int)ceill(v); return a; }
int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage: %s KMAX HCAP SLACK [RECON_KMAX [RECON_KMIN]]\n",argv[0]);return 1;}
  int KMAX=atoi(argv[1]); int64_t HCAP=atoll(argv[2]); int SLACK=atoi(argv[3]);
  int RK=argc>4?atoi(argv[4]):0; int RKMIN=argc>5?atoi(argv[5]):1;
  TH=logl(2.0L)/logl(3.0L);
  int64_t N=(3*HCAP)/2+2;
  uint16_t *cA=malloc(2*N),*nA=malloc(2*N); uint32_t *cH=malloc(4*N),*nH=malloc(4*N);
  uint16_t **LA=NULL;
  if(!cA||!nA||!cH||!nH){fprintf(stderr,"oom\n");return 1;}
  if(RK>0){ LA=malloc(sizeof(uint16_t*)*(RK+1)); for(int b=0;b<=RK;b++){ LA[b]=malloc(2*N); if(!LA[b]){fprintf(stderr,"oom recon\n");return 1;} } }
  memset(cA,255,2*N); memset(cH,0,4*N); cA[1]=0; cH[1]=0;
  if(RK>0) memcpy(LA[0],cA,2*N);
  int64_t *minW=malloc(sizeof(int64_t)*70000);
  printf("# mirror_loops_dp KMAX=%d HCAP=%lld SLACK=%d  (theta=log_3 2)\n",KMAX,(long long)HCAP,SLACK);
  for(int b=0;b<KMAX;b++){
    int L=b+1;  // layer of the new states
    // minW[a] = least |w| kept at layer L with multiplication count a (permissive by 1e-9)
    long double base=(long double)(L+1)*TH+1.0L+(long double)SLACK;
    for(int a=0;a<70000;a++){
      if(SLACK<0){ minW[a]=0; continue; }
      long double T=(long double)a-base;
      if(T<0){ minW[a]=0; continue; }
      long double w=powl(3.0L,T)*(1.0L-1e-9L);
      minW[a]= (w>(long double)N)? N+1 : (int64_t)floorl(w);
    }
    memset(nA,255,2*N);
    long long live=0;
    for(int64_t i=1;i<N;i++){
      uint16_t a=cA[i]; if(a==NONE) continue; live++;
      uint32_t h=cH[i];
      int64_t x=-i; int aa=a; int64_t hh=h;
      if(x&1){ if(-x>HCAP) continue; if(-x>hh) hh=-x; x=3*x+1; aa++; }
      for(;;){
        int64_t w=-(x/2);   // |x/2|, x even and negative
        if(aa<70000 && w>minW[aa]){   // keep iff |w| > 3^(a-(L+1)theta-1-SLACK) (permissive)
          if(aa<nA[w] || (aa==nA[w] && hh<nH[w])){ nA[w]=(uint16_t)aa; nH[w]=(uint32_t)hh; }
        }
        // two more multiplications: mult points x and 3x+1
        int64_t m2=-(3*x+1); if(-x>HCAP || m2>HCAP) break;
        if(m2>hh) hh=m2;
        x=9*x+4; aa+=2;
      }
    }
    int K=L; int a0=a0_of(K);
    long double eta=(long double)a0-(long double)(K+1)*TH;
    if(nA[1]==NONE) printf("K=%5d a0=%4d amin=none(<=cap/slack) eta=%.7Lf live=%lld\n",K,a0,eta,live);
    else {
      long double r=exp2l((long double)nA[1]*log2l(3.0L)-(long double)K);  // 3^amin/2^K without overflow
      printf("K=%5d a0=%4d amin=%4d ratio=%.6Lf eta=%.7Lf H=%u live=%lld %s\n",K,a0,nA[1],r,eta,nH[1],live,
        nA[1]==a0?"=a0":(nA[1]<a0?"<a0":">a0"));
    }
    fflush(stdout);
    uint16_t*t1=cA;cA=nA;nA=t1; uint32_t*t2=cH;cH=nH;nH=t2;
    if(RK>0 && L<=RK) memcpy(LA[L],cA,2*N);
    if(RK>0 && L==RK) break;
  }
  if(RK>0){
    // reconstruction: for each K in [RKMIN,RK], backtrack from (K,-1)
    for(int K=RKMIN;K<=RK;K++){
      if(LA[K][1]==NONE){ printf("RECON K=%d none\n",K); continue; }
      int nmoves=0; char *word=malloc(16*(K+1)*64); word[0]=0; int64_t w=-1; int A=LA[K][1]; int ok=1; int64_t hmax=0;
      // collect runs in reverse, then print forward
      int *tlist=malloc(sizeof(int)*(K+1)); int64_t *ulist=malloc(sizeof(int64_t)*(K+1));
      for(int L=K;L>=1;L--){
        // predecessors u of w at layer L-1: 2w = M^t(u)
        int64_t x=2*w; int found=0;
        for(int t=0;t<=60;t++){
          // u = M^{-t}(x); compute by undoing t multiplications
          int64_t u=x; int good=1; int64_t mx=0;
          for(int s=0;s<t;s++){ int64_t y=u-1; if(y%3!=0){good=0;break;} u=y/3; if(-u>mx) mx=-u; }
          if(!good) break;             // further t also fail (need integrality at each step)
          if(t>0 && mx>HCAP) continue;
          if(((u&1)!=0) != ((t&1)!=0)) continue;   // parity of u must equal parity of t
          if(-u>=N || -u<1) continue;
          uint16_t au=LA[L-1][-u];
          if(au!=NONE && au+t==A){ tlist[L]=t; ulist[L]=u; A=au; w=u; found=1; break; }
        }
        if(!found){ ok=0; break; }
      }
      if(!ok || w!=-1){ printf("RECON K=%d FAILED\n",K); free(tlist); free(ulist); free(word); continue; }
      // forward word: for L=1..K: M^{t_L} H ; compute height and verify
      int64_t v=-1; int aa=0; char buf[64]; int len=0; char last=0; int run=0;
      #define FLUSH() do{ if(run>0){ len+=sprintf(word+len,"%c%d ",last,run); } }while(0)
      for(int L=1;L<=K;L++){
        for(int s=0;s<tlist[L];s++){ if(-v>hmax) hmax=-v; v=3*v+1; aa++; if(last!='M'){FLUSH(); last='M'; run=0;} run++; }
        if(v&1){ ok=0; break; } v/=2; if(last!='H'){FLUSH(); last='H'; run=0;} run++;
      }
      FLUSH();
      printf("RECON K=%d a=%d H=%lld end=%lld %s word: %s\n",K,aa,(long long)hmax,(long long)v,(ok&&v==-1)?"OK":"BAD",word);
      free(tlist); free(ulist); free(word);
    }
  }
  return 0;
}
