/* verify_chamber.c - kind-pasteur-2026-11-07-S14, THM-4904 checks (independent of rc_census.c: no closed-form tail,
 * every word is built letter by letter from the rooted edge and stored as a string).
 *
 * Physical moves on s = (D,E,U,d): O = (D,0,2U,d), Z = (D-E+U,d,2U-d,d), J = (D-E+U,D-E,D+U,D); t = (1,0,2,1).
 * F(s) = D-E+U is the value of the word extended by Z (equivalently by J).
 * M_K(N) = #words with <= K letters J, ending in Z or J, with value N; plus 1 for the empty word at N = 1.
 *
 * Checks, for K = 0,1,2 and every odd N <= X:
 *  (A) adjoint: F(t.w) = F(t.w^dagger) for every word w (any number of J) of length <= L  (dagger = reverse, O<->Z, J fixed)
 *  (B) decomposition: M_K(N) = [N=1] + 2 A_(K-1)(N) + B_K(N)  (A_j: words w with <= j J and F(t.w)=N; B_K: exactly K J)
 *  (C) parity law: M_K(N) = [N=1] + Fix_K(N) (mod 2), Fix_K(N) = #{w : w = w^dagger, exactly K J, F(t.w) = N}
 *  (D) chamber: if M_K(N) = 1 and N > 1, the unique word is wZ with w = w^dagger, exactly K J, w not starting or ending in J
 * Prints, per K, the singleton values (M_K = 1) up to X with their words, and window histograms for comparison with
 * rc_census.c.   usage: verify_chamber X L
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef long long ll;
static ll X; static int KMAX = 2;
static int *M[3], *A[3], *B[3], *FX[3];      /* indexed by N (odd N only used) */
static char wbuf[256];
static char *sing_word[3]; static ll *sing_val[3]; static int nsing[3];
static ll states=0, phys_fail=0, growth_fail=0; static double c0;  /* physical inequalities, growth lemma */
static char **witness[3];                     /* one representing prefix w (for B_K words) per N, to print singletons */

static int is_dagger_fixed(const char *w, int n){
  for (int i=0;i<n;i++){ char a=w[i], b=w[n-1-i]; char sb = b=='O'?'Z':(b=='Z'?'O':'J'); if (a!=sb) return 0; }
  return 1;
}
/* enumerate prefixes w (as states) with F(t.w) <= X; j = #J in w */
static void dfs(ll D, ll E, ll U, ll d, int j, int len){
  ll F = D-E+U;
  if (F > X) return;                          /* F strictly increases along every extension */
  states++;
  if (!(D-E >= d && d >= 1 && D <= U && U >= 2*d)) phys_fail++;                 /* Session200 section 1 */
  { int noz = len - j; double lb = 2.0*pow(c0, j+1)*ldexp(1.0, noz);            /* U >= 2 c0^(#J+1) 2^(#O/Z) */
    if ((double)U < lb*(1.0-1e-12)) growth_fail++; }
  if (F & 1){
    for (int K=0;K<=KMAX;K++){
      if (j <= K-1) A[K][F]++;                /* wJ and wZ both counted in M_K */
      if (j == K){ B[K][F]++; if (is_dagger_fixed(wbuf,len)) FX[K][F]++;
        if (!witness[K][F]){ witness[K][F]=malloc(len+2); memcpy(witness[K][F],wbuf,len); witness[K][F][len]='Z'; witness[K][F][len+1]=0; } }
    }
  }
  if (len >= 250) { fprintf(stderr,"word too long\n"); exit(1); }
  if (j < KMAX){ wbuf[len]='J'; dfs(D-E+U, D-E, D+U, D, j+1, len+1); }
  wbuf[len]='O'; dfs(D, 0, 2*U, d, j, len+1);
  wbuf[len]='Z'; dfs(D-E+U, d, 2*U-d, d, j, len+1);
}
/* value of t.w followed by Z, for an arbitrary word */
static ll Fword(const char *w, int n){ ll D=1,E=0,U=2,d=1; for(int i=0;i<n;i++){ ll D2,E2,U2,d2;
  if(w[i]=='O'){D2=D;E2=0;U2=2*U;d2=d;} else if(w[i]=='Z'){D2=D-E+U;E2=d;U2=2*U-d;d2=d;} else {D2=D-E+U;E2=D-E;U2=D+U;d2=D;}
  D=D2;E=E2;U=U2;d=d2; } return D-E+U; }
static ll adj_checked=0, adj_fail=0;
static void adjoint_check(char *w, int n, int L){
  if (n>0){ char v[64]; for(int i=0;i<n;i++){ char b=w[n-1-i]; v[i]= b=='O'?'Z':(b=='Z'?'O':'J'); }
    adj_checked++; if (Fword(w,n)!=Fword(v,n)) { adj_fail++; if(adj_fail<5){ w[n]=0; fprintf(stderr,"adjoint FAIL %s\n",w);} } }
  if (n==L) return; const char *al="OZJ"; for(int a=0;a<3;a++){ w[n]=al[a]; adjoint_check(w,n+1,L); }
}
int main(int argc, char **argv){
  if (argc<3){ fprintf(stderr,"usage: verify_chamber X L\n"); return 1; }
  X = atoll(argv[1]); int L = atoi(argv[2]);
  c0 = 1.0; for (int i=0;i<2000;i++) c0 *= 1.0 - 0.25*pow(2.0/3.0, i);
  char w[64]; adjoint_check(w,0,L);
  printf("(A) adjoint F(t.w)=F(t.w^dagger): %lld words of length 1..%d, failures %lld\n", adj_checked, L, adj_fail);
  for (int K=0;K<=KMAX;K++){ A[K]=calloc(X+2,sizeof(int)); B[K]=calloc(X+2,sizeof(int)); FX[K]=calloc(X+2,sizeof(int)); M[K]=calloc(X+2,sizeof(int)); witness[K]=calloc(X+2,sizeof(char*)); }
  dfs(1,0,2,1,0,0);
  printf("physical inequalities D-E>=d>=1, D<=U, U>=2d on %lld visited states: failures %lld\n", states, phys_fail);
  printf("growth lemma U >= 2 c0^(#J+1) 2^(#O/Z), c0 = %.9f: failures %lld\n", c0, growth_fail);
  int bad = (phys_fail||growth_fail);
  for (int K=0;K<=KMAX;K++){
    ll parity_fail=0, chamber_fail=0, ns=0;
    for (ll N=1;N<=X;N+=2){
      M[K][N] = (N==1) + 2*A[K][N] + B[K][N];             /* (B) by construction of the count */
      if ( ((M[K][N] - (N==1) - FX[K][N]) & 1) ) parity_fail++;          /* (C) */
      if (M[K][N]==1 && N>1){ ns++;
        /* (D): the unique word must be a B_K word, dagger-fixed, not starting/ending with J */
        const char *ws = witness[K][N]; int n = ws? (int)strlen(ws)-1 : -1;
        int ok = A[K][N]==0 && B[K][N]==1 && FX[K][N]==1 && n>=0 && (n==0 || (ws[0]!='J' && ws[n-1]!='J'));
        if (!ok) chamber_fail++;
      }
    }
    ll nch=0; for (ll N=1;N<=X;N+=2) nch += FX[K][N];
    printf("K=%d: (C) parity failures %lld; singletons N>1: %lld; (D) chamber failures %lld; chamber words (dagger-fixed, exactly K J, value <= X): %lld = %.4f sqrt(X)\n",K,parity_fail,ns,chamber_fail,nch,nch/sqrt((double)X));
    if (parity_fail||chamber_fail) bad=1;
    printf("K=%d singleton values (first 40) with word w (value of wZ): ",K); int shown=0;
    for (ll N=3;N<=X && shown<40;N+=2) if (M[K][N]==1){ printf(" %lld:%s",N,witness[K][N]); shown++; }
    printf("\n");
    for (ll lo=1; lo<X; lo*=2){ ll hi = lo*2>X? X+1 : lo*2; ll c[9]={0}; int mn=1<<30; ll n=0;
      for (ll N=lo|1; N<hi; N+=2){ int m=M[K][N]>255?255:M[K][N]; n++; if(m<mn)mn=m; if(m<=8)c[m]++; }
      printf("  K=%d [%lld,%lld): n=%lld min=%d c0..8=",K,lo,hi,n,mn); for(int m=0;m<=8;m++) printf("%s%lld",m?",":"",c[m]); printf("\n"); }
  }
  printf("RESULT %s\n", (bad||adj_fail)? "FAIL" : "PASS");
  return (bad||adj_fail)?1:0;
}
