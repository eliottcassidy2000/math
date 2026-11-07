/* rc_lf.c - kind-pasteur-2026-11-07-S15 (tests R2-R5). Exact two-root-change word multiplicity M_K(N) for odd N <= X
 * (same enumeration as the S14 rc_census.c; 16-bit counters, exact for X <= 2^24 where max M_2 = 3792), then the mean
 * of M_K over odd N in [LO, HI) by class mod 27, 5, 7, 11, 13, 17, 23, 31, and the overall mean.
 * usage: rc_lf X K LO HI
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
typedef long long ll;
static ll X; static int K; static unsigned short *cnt;
static inline void mark(ll v){ if(v>=1&&v<=X&&(v&1)) __atomic_fetch_add(&cnt[v>>1],(unsigned short)1,__ATOMIC_RELAXED); }
typedef struct { ll h,R,e,d; } tailctx;
static void tdfs(const tailctx*c, ll m, ll q){
  if (c->h+(c->R-c->d)*m > X) return;
  mark(c->h+c->R*m-c->d*q-c->e*(m&1));
  tdfs(c,2*m,q); tdfs(c,2*m+1,q+m+(m&1));
}
static void tail(ll h, ll R, ll e, ll d){ if(h>X) return; tailctx c={h,R,e,d}; mark(h); tdfs(&c,1,0); }
static void block(ll D, ll E, ll U, ll d, int j){
  if (j < K){ ll D2=D-E+U, E2=D-E, U2=D+U, d2=D; if (D2<=X){ tail(D2, U2, E2, d2); block(D2,E2,U2,d2,j+1); } }
  if (j < K){ if (2*U <= X) block(D,0,2*U,d,j); if (2*U-d <= X) block(D-E+U,d,2*U-d,d,j); }
}
typedef struct { int type,j; ll D,E,U,d; } task;
static task *T; static size_t nT=0, capT=0; static ll USPLIT=512;
static void push(int type,int j,ll D,ll E,ll U,ll d){ if(nT==capT){capT=capT?2*capT:1024; T=realloc(T,capT*sizeof(task));} T[nT++]=(task){type,j,D,E,U,d}; }
static void expand(ll D, ll E, ll U, ll d, int j){
  if (U >= USPLIT){ push(1,j,D,E,U,d); return; }
  if (j < K){ ll D2=D-E+U, E2=D-E, U2=D+U, d2=D; if (D2<=X){ push(0,j+1,D2,E2,U2,d2); expand(D2,E2,U2,d2,j+1); } }
  if (j < K){ if (2*U <= X) expand(D,0,2*U,d,j); if (2*U-d <= X) expand(D-E+U,d,2*U-d,d,j); }
}
int main(int argc,char**argv){
  if(argc<5){fprintf(stderr,"usage: rc_lf X K LO HI\n");return 1;}
  X=atoll(argv[1]); K=atoi(argv[2]); ll LO=atoll(argv[3]), HI=atoll(argv[4]);
  cnt=calloc((X>>1)+2,2);
  push(0,0,1,0,2,1); expand(1,0,2,1,0);
  #pragma omp parallel for schedule(dynamic,1)
  for (ll i=0;i<(ll)nT;i++){ task*k=&T[i]; if(k->type==0) tail(k->D,k->U,k->E,k->d); else block(k->D,k->E,k->U,k->d,k->j); }
  int mods[8]={27,5,7,11,13,17,23,31};
  double all=0, nall=0; int maxm=0;
  for (ll v=LO|1; v<HI; v+=2){ all+=cnt[v>>1]; nall++; if(cnt[v>>1]>maxm) maxm=cnt[v>>1]; }
  printf("rc_lf X=%lld K=%d window [%lld,%lld): odd N %.0f, overall mean %.4f, max M %d\n",X,K,LO,HI,nall,all/nall,maxm);
  for (int t=0;t<8;t++){ int q=mods[t]; double *s=calloc(q,sizeof(double)), *n=calloc(q,sizeof(double));
    for (ll v=LO|1; v<HI; v+=2){ s[v%q]+=cnt[v>>1]; n[v%q]++; }
    printf("mod %d class means:",q); for(int r=0;r<q;r++) printf(" %.4f", s[r]/n[r]); printf("\n");
    free(s); free(n); }
  return 0;
}
