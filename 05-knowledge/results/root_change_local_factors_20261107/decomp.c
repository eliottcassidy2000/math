/* decomp.c - kind-pasteur-2026-11-07-S15. Exact decomposition of the two-root-change word multiplicity M_K(N) by the
 * 3-adic confinement level of the seed that produced each representation.
 * Seeds: the root t and every state right after a root change J (levels 1..K). A seed s = (h,R,e,d) produces the
 * values H_s(m) = h + R m - d Q(m) - e (m mod 2), m >= 0 (m = 0 is the word ending in J, or the empty word for t).
 * Confinement lemma: H_s(m) = h - e (m mod 2) (mod gcd(R,d)). For a J-seed, gcd(R,d) = gcd(D,U) of the pre-J state.
 * For odd N in the window [LO, HI) this program accumulates, for each confinement level j = v3(gcd(R,d)) (capped at
 * JMAX) and parity b = m mod 2, the number of representations in each class N mod 3^6, and also the total by
 * v3(N), N mod 5 and N mod 7 split by v5(gcd), v7(gcd) confinement. Output is deterministic text.
 * usage: decomp X K LO HI   (enumerates words with values <= X; only values in [LO,HI) are counted)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
typedef long long ll;
#define JMAX 6
#define M3 729            /* 3^6 */
static ll X, LO, HI; static int K;
typedef struct { ll c3[JMAX+1][2][M3]; ll c5[3][2][5]; ll c7[3][2][7]; } acc_t;
static acc_t *ACC;        /* one per thread */
static int v3(ll g){ int j=0; while(g%3==0 && j<JMAX){ g/=3; j++; } return j; }
static int vp(ll g, int p){ int j=0; while(g%p==0 && j<2){ g/=p; j++; } return j; }
static ll gcdll(ll a, ll b){ while(b){ ll t=a%b; a=b; b=t; } return a; }
typedef struct { ll h,R,e,d; int j3,j5,j7; acc_t *A; } ctx;
static inline void rec(const ctx*c, ll v, int b){
  if (v < LO || v >= HI) return;
  c->A->c3[c->j3][b][v % M3]++;
  c->A->c5[c->j5][b][v % 5]++;
  c->A->c7[c->j7][b][v % 7]++;
}
static void tdfs(const ctx*c, ll m, ll q){
  if (c->h + (c->R - c->d)*m > X) return;
  ll v = c->h + c->R*m - c->d*q - c->e*(m&1);
  rec(c, v, (int)(m&1));
  tdfs(c, 2*m, q); tdfs(c, 2*m+1, q + m + (m&1));
}
static void tail(ll h, ll R, ll e, ll d){
  if (h > X) return;
  ctx c; c.h=h; c.R=R; c.e=e; c.d=d; ll g=gcdll(R,d); c.j3=v3(g); c.j5=vp(g,5); c.j7=vp(g,7);
  c.A = &ACC[omp_get_thread_num()];
  rec(&c, h, 0); tdfs(&c, 1, 0);
}
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
int main(int argc, char**argv){
  if (argc<5){ fprintf(stderr,"usage: decomp X K LO HI\n"); return 1; }
  X=atoll(argv[1]); K=atoi(argv[2]); LO=atoll(argv[3]); HI=atoll(argv[4]);
  int nth = omp_get_max_threads(); ACC = calloc(nth, sizeof(acc_t));
  push(0,0,1,0,2,1); expand(1,0,2,1,0);
  #pragma omp parallel for schedule(dynamic,1)
  for (ll i=0;i<(ll)nT;i++){ task*k=&T[i];
    if (k->type==0) tail(k->D, k->U, k->E, k->d); else block(k->D,k->E,k->U,k->d,k->j); }
  acc_t S; memset(&S,0,sizeof S);
  for (int t=0;t<nth;t++){ for(int j=0;j<=JMAX;j++) for(int b=0;b<2;b++) for(int r=0;r<M3;r++) S.c3[j][b][r]+=ACC[t].c3[j][b][r];
    for(int j=0;j<3;j++) for(int b=0;b<2;b++){ for(int r=0;r<5;r++) S.c5[j][b][r]+=ACC[t].c5[j][b][r]; for(int r=0;r<7;r++) S.c7[j][b][r]+=ACC[t].c7[j][b][r]; } }
  /* number of odd N in [LO,HI) per class mod 729, 5, 7 */
  static ll n3[M3]; ll n5[5]={0}, n7[7]={0};
  for (ll v=(LO|1); v<HI; v+=2){ n3[v%M3]++; n5[v%5]++; n7[v%7]++; }
  printf("decomp X=%lld K=%d window [%lld,%lld)\n", X, K, LO, HI);
  /* mean multiplicity contributed by confinement level j, per class r mod 3^6, aggregated: by v3(r) and by r mod 9 */
  for (int j=0;j<=JMAX;j++) for (int b=0;b<2;b++){
    double byv[7]={0}, nv[7]={0}; double by9[9]={0}, n9[9]={0}; ll tot=0;
    for (int r=0;r<M3;r++){ if(!n3[r]) continue; int k=0; int rr=r; while(rr%3==0 && k<6 && rr){ rr/=3; k++; } if(r==0) k=6;
      byv[k]+=S.c3[j][b][r]; nv[k]+=n3[r]; by9[r%9]+=S.c3[j][b][r]; n9[r%9]+=n3[r]; tot+=S.c3[j][b][r]; }
    printf("j3=%d b=%d total=%lld  mean by v3(N)=0..6:", j, b, tot);
    for (int k=0;k<7;k++) printf(" %.3f", nv[k]? byv[k]/nv[k] : 0.0);
    printf("  | by N mod 9:"); for (int r=0;r<9;r++) printf(" %.3f", n9[r]? by9[r]/n9[r] : 0.0); printf("\n");
  }
  /* uniformity of the unconfined part (j3 = 0) across all 729 classes */
  { double mn=1e300, mx=0; for (int r=0;r<M3;r++){ if(!n3[r]) continue; double m=(double)(S.c3[0][0][r]+S.c3[0][1][r])/n3[r]; if(m<mn)mn=m; if(m>mx)mx=m; }
    printf("unconfined (j3=0) mean per class mod 729: min %.4f max %.4f ratio %.5f\n", mn, mx, mx/mn); }
  /* totals by v3(N) and overall */
  { double byv[7]={0}, nv[7]={0}; double all=0, nall=0;
    for (int r=0;r<M3;r++){ if(!n3[r]) continue; int k=0; int rr=r; while(rr%3==0 && k<6 && rr){ rr/=3; k++; } if(r==0) k=6;
      double s=0; for(int j=0;j<=JMAX;j++) for(int b=0;b<2;b++) s+=S.c3[j][b][r]; byv[k]+=s; nv[k]+=n3[r]; all+=s; nall+=n3[r]; }
    printf("TOTAL mean M_K by v3(N)=0..6:"); for(int k=0;k<7;k++) printf(" %.3f", byv[k]/nv[k]); printf("  overall %.3f\n", all/nall); }
  for (int p=5;p<=7;p+=2){ printf("p=%d: mean by N mod p (all seeds):", p);
    for (int r=0;r<p;r++){ double s=0; for(int j=0;j<3;j++) for(int b=0;b<2;b++) s+= (p==5? S.c5[j][b][r] : S.c7[j][b][r]); printf(" %.3f", s/(p==5? n5[r]:n7[r])); }
    printf("  | unconfined only:"); for (int r=0;r<p;r++){ double s=0; for(int b=0;b<2;b++) s+= (p==5? S.c5[0][b][r] : S.c7[0][b][r]); printf(" %.3f", s/(p==5? n5[r]:n7[r])); }
    printf("\n"); }
  return 0;
}
