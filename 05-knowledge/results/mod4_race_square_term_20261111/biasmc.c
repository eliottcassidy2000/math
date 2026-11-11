// biasmc.c - Monte Carlo with exact conditional expectation for the mod-4 law one order beyond the censuses.
// Samples N uniformly random LABELLED tournaments R on m vertices (independent fair arcs; splitmix64 seeded by argv),
// and for each computes exactly, as in biasweights.c:
//   H(R), eps_R = (-1)^((H(R)-1)/2);  mu, beta = mu/2 mod 2, ell = a/2 mod 2;  rank(beta), radical test;
//   diff_R = sum over all 2^m out-sets X of (-1)^(q_R(X))   (Gray code; b_R = diff_R / 2^m).
// Then  D_(m+1) = E[eps_R b_R]  (unbiased estimator, per-sample variance E[b_R^2] ~ 2^-m),  rho_m = 2^m E[b_R^2],
// E|b_R|, the rank law of beta_R and the non-zero-bias frequency by rank.
// Usage: biasmc m N seed      Output: one line of sums (parsed by merge_mc.py).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXN 13
static int n; static u64 f[1<<MAXN][MAXN], g[1<<MAXN][MAXN]; static int outm[MAXN], inm[MAXN];
static u64 sm_state;
static inline u64 splitmix(void){ u64 z=(sm_state+=0x9E3779B97F4A7C15ULL); z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31); }
int main(int argc,char**argv){ n=atoi(argv[1]); u64 N=strtoull(argv[2],0,10); sm_state=strtoull(argv[3],0,10)*0x2545F4914F6CDD1DULL+12345;
  int full=(1<<n)-1; static u64 M[MAXN][MAXN];
  u64 cls[MAXN+1][2]; for(int r=0;r<=MAXN;r++){cls[r][0]=cls[r][1]=0;}
  u64 fail_even=0, fail_a=0, fail_c3=0, fail_edges=0; long double S1=0, S2=0, SA=0, SU=0, SQ=0; i64 npos=0;
  for(u64 it=0; it<N; it++){
    for(int i=0;i<n;i++){outm[i]=0;inm[i]=0;}
    u64 bits=0; int nb=0;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){ if(nb==0){bits=splitmix();nb=64;} if(bits&1){outm[i]|=1<<j;inm[j]|=1<<i;} else {outm[j]|=1<<i;inm[i]|=1<<j;} bits>>=1; nb--; }
    for(int S=1;S<=full;S++){ for(int v=0;v<n;v++){ if(!(S>>v&1)){f[S][v]=g[S][v]=0;continue;} int R=S&~(1<<v);
        if(R==0){f[S][v]=g[S][v]=1;continue;} u64 a=0,b=0; int mm=R&inm[v]; while(mm){int w=__builtin_ctz(mm);mm&=mm-1;a+=f[R][w];}
        mm=R&outm[v]; while(mm){int w=__builtin_ctz(mm);mm&=mm-1;b+=g[R][w];} f[S][v]=a;g[S][v]=b; } }
    u64 H=0; for(int v=0;v<n;v++) H+=f[full][v];
    for(int j=0;j<n;j++)for(int i=0;i<n;i++){ if(i==j){M[j][i]=0;continue;} u64 a=0; int rest=full&~(1<<i)&~(1<<j);
      for(int S0=rest;;S0=(S0-1)&rest){int S=S0|(1<<j); a+=f[S][j]*g[full^S][i]; if(S0==0)break;} M[j][i]=a; }
    int rows[MAXN], ell[MAXN]; int edges=0;
    for(int i=0;i<n;i++){ rows[i]=0; for(int j=0;j<n;j++){ if(i==j)continue; u64 mu=M[i][j]+M[j][i]; if(mu&1)fail_even++; if((mu>>1)&1){ rows[i]|=1<<j; if(j>i) edges++; } }
      i64 a=(i64)g[full][i]-(i64)f[full][i]; for(int j=0;j<n;j++) if(j!=i) a+=(i64)M[j][i]; if(a&1) fail_a++; ell[i]=(int)((a>>1)&1); }
    if(!(edges&1)) fail_edges++;
    int rr[MAXN], id[MAXN]; for(int i=0;i<n;i++){rr[i]=rows[i]; id[i]=1<<i;} int rk=0;
    for(int c=0;c<n;c++){ int p=-1; for(int i=rk;i<n;i++) if(rr[i]>>c&1){p=i;break;} if(p<0)continue;
      int t=rr[p];rr[p]=rr[rk];rr[rk]=t; t=id[p];id[p]=id[rk];id[rk]=t;
      for(int i=0;i<n;i++) if(i!=rk && (rr[i]>>c&1)){ rr[i]^=rr[rk]; id[i]^=id[rk]; } rk++; }
    int radzero=1;
    for(int i=rk;i<n;i++){ int X=id[i]; int q=0; for(int a=0;a<n;a++) if(X>>a&1){ q^=ell[a]; q^=__builtin_popcount(rows[a]&X&((1<<a)-1))&1; } if(q){radzero=0;} }
    i64 diff=1; { int X=0, q=0; for(int t=1;t<=full;t++){ int i=__builtin_ctz(t); q^=ell[i]^(__builtin_popcount(rows[i]&X&~(1<<i))&1); X^=1<<i; diff += q?-1:1; } }
    i64 ad=diff<0?-diff:diff; i64 pred=1LL<<(n-rk/2);
    if((rk&1) || !((ad==0&&!radzero)||(ad==pred&&radzero))) fail_c3++;
    cls[rk][ad!=0]++;
    i64 eps=(((H-1)>>1)&1)?-1:1; long double b=(long double)diff/(long double)(1LL<<n);
    S1+=eps*b; S2+=b*b; SA+=(ad/(long double)(1LL<<n)); SU+=b; SQ+=(eps*b)*(eps*b); if(eps>0) npos++;
  }
  printf("m %d N %llu failures(mu odd, a odd, C3, edges) %llu %llu %llu %llu S1 %.12Le S2 %.12Le SA %.12Le SU %.12Le npos %lld", n, N, fail_even, fail_a, fail_c3, fail_edges, S1, S2, SA, SU, npos);
  for(int r=0;r<=n;r++) printf(" r%d %llu %llu", r, cls[r][0], cls[r][1]);
  printf("\n"); return 0; }
