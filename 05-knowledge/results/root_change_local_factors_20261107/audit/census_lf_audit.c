/* census_lf_audit.c - independent auditor's census for the S15 packet (THM-4916 / HYP-9506).
 *
 * Written from the definitions only (it does not reuse rc_lf.c / decomp.c / census_audit.c):
 *   words over {O,Z,J} act on t = (D,E,U,d) = (1,0,2,1) by
 *     O(D,E,U,d) = (D,0,2U,d),  Z = (D-E+U,d,2U-d,d),  J = (D-E+U,D-E,D+U,D);
 *   value(word) = final D, F(s) = D - E + U, so value(wZ) = value(wJ) = F(t.w);
 *   M_2(N) = [N=1] + #{wJ : #J(w) <= 1, F(t.w) = N} + #{wZ : #J(w) <= 2, F(t.w) = N}.
 * Enumeration: plain letter-by-letter DFS over prefixes w with #J(w) <= 2. Pruning uses only the fact that F strictly
 * increases along every letter (F(sO)-F(s) = U+E, F(sZ)-F(s) = 2U-2d, F(sJ)-F(s) = U+E, with E >= 0 and U > d >= 1;
 * these inequalities are asserted at every node when CHECK is defined).
 *
 * Decomposition (task 4): every counted word is (seed)(tail), seed = t or the state right after the LAST J, tail = the
 * O/Z letters after it, read as m = sum b_i 2^i with b_0 the first tail letter (O=0, Z=1); m = 0 is the empty tail
 * (word ending in J, or the empty word). Seed confinement level j3 = v3(gcd(U_s, d_s)) capped at 6 (for a J-seed this
 * is v3(gcd(D,U)) of the pre-J state). Q(m) mod 3 is tracked incrementally by Q(m + 2^L) = Q(m) + b_{L-1} + rev_L(m)
 * (checked against the recursive definition in identities_audit.py).
 *
 * usage: census_lf_audit X LO HI nthreads [lowcap]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
typedef long long ll;
typedef unsigned long long ull;

static ll X, LO, HI; static int LOWCAP = 12;
#define M729 729
typedef struct {
  unsigned int *cnt;          /* M_2 counts for odd N <= X, index N>>1 */
  unsigned int *ucnt;         /* per-N unconfined (j3 = 0, root included) counts in the window */
  ull c3[7][2][M729];         /* seed level j3, tail parity b, N mod 729 (window only) */
  ull c5[3][2][5], c7[3][2][7];
  ull jz_tc[2][M729];         /* J-seeds with j3 = 0, odd tails: [tail-confined mod 3][N mod 729] */
  ull jz_even[M729];          /* J-seeds with j3 = 0, even tails (m even, including m = 0) */
  ull root[2][M729];          /* root seed t, by tail parity */
  ull words;                  /* all counted words with value <= X (excluding the empty word) */
} acc_t;
static acc_t *ACC;

typedef struct {
  ll D, E, U, d; int j;       /* state and number of J in the prefix */
  int j3, j5, j7, isroot;     /* current seed */
  int L, b0, last, m3, Q3, rev3, pw3;  /* tail bookkeeping (mod 3) */
} node;

static inline int vcap(ll x, int p, int cap){ int k = 0; while (k < cap && x % p == 0){ x /= p; k++; } return k; }
static inline int imin(int a, int b){ return a < b ? a : b; }

static inline void record(acc_t *A, ll N, int j3, int j5, int j7, int isroot, int b, int tc){
  A->cnt[N >> 1]++; A->words++;
  if (N < LO || N >= HI) return;
  int r = (int)(N % M729);
  A->c3[j3][b][r]++; A->c5[j5][b][N % 5]++; A->c7[j7][b][N % 7]++;
  if (j3 == 0) A->ucnt[(N - LO) >> 1]++;
  if (isroot) A->root[b][r]++;
  else if (j3 == 0){ if (b) A->jz_tc[tc][r]++; else A->jz_even[r]++; }
}

/* record the counted words at node w, then recurse */
static void visit(acc_t *A, const node *w){
  ll D = w->D, E = w->E, U = w->U, d = w->d;
#ifdef CHECK
  if (!(E >= 0 && U > d && d >= 1 && D - E >= d && U >= 2*d)){ fprintf(stderr, "phys violated\n"); exit(2); }
#endif
  ll N = D - E + U;                    /* = F(t.w) <= X guaranteed by caller */
  /* word wZ */
  {
    int b = (w->L == 0) ? 1 : w->b0;
    int tc = 0;
    if (!w->isroot && b == 1){
      int mZ = (w->m3 + w->pw3) % 3;
      int QZ = (w->Q3 + (w->L ? (w->last + w->rev3) : 0)) % 3;
      tc = ((mZ + 1) % 3 == 0) && ((QZ + 1) % 3 == 0);
    }
    record(A, N, w->j3, w->j5, w->j7, w->isroot, b, tc);
  }
  /* word wJ */
  if (w->j <= 1){
    int nj3 = imin(imin(vcap(D,3,6), vcap(U,3,6)), 6);
    int nj5 = imin(vcap(D,5,2), vcap(U,5,2));
    int nj7 = imin(vcap(D,7,2), vcap(U,7,2));
    record(A, N, nj3, nj5, nj7, 0, 0, 0);
  }
  /* children */
  node c;
  /* O */
  if (D + 2*U <= X){
    c = *w; c.D = D; c.E = 0; c.U = 2*U; c.d = d;
    c.b0 = (w->L == 0) ? 0 : w->b0; c.last = 0; c.rev3 = (2*w->rev3) % 3; c.L = w->L + 1; c.pw3 = (2*w->pw3) % 3;
    visit(A, &c);
  }
  /* Z */
  if (D - E + 3*U - 2*d <= X){
    c = *w; c.D = D - E + U; c.E = d; c.U = 2*U - d; c.d = d;
    c.m3 = (w->m3 + w->pw3) % 3; c.Q3 = (w->Q3 + (w->L ? (w->last + w->rev3) : 0)) % 3;
    c.b0 = (w->L == 0) ? 1 : w->b0; c.last = 1; c.rev3 = (2*w->rev3 + 1) % 3; c.L = w->L + 1; c.pw3 = (2*w->pw3) % 3;
    visit(A, &c);
  }
  /* J */
  if (w->j < 2 && D + 2*U <= X){
    c.D = D - E + U; c.E = D - E; c.U = D + U; c.d = D; c.j = w->j + 1;
    c.j3 = imin(imin(vcap(D,3,6), vcap(U,3,6)), 6); c.j5 = imin(vcap(D,5,2), vcap(U,5,2)); c.j7 = imin(vcap(D,7,2), vcap(U,7,2));
    c.isroot = 0; c.L = 0; c.b0 = 0; c.last = 0; c.m3 = 0; c.Q3 = 0; c.rev3 = 0; c.pw3 = 1;
    visit(A, &c);
  }
}

/* frontier generation: nodes at depth DEP are collected as tasks; shallower nodes are visited sequentially */
static node *T; static size_t nT = 0, capT = 0;
static void collect(acc_t *A, const node *w, int depth, int DEP){
  if (depth == DEP){ if (nT == capT){ capT = capT ? 2*capT : 4096; T = realloc(T, capT*sizeof(node)); } T[nT++] = *w; return; }
  ll D = w->D, E = w->E, U = w->U, d = w->d;
  ll N = D - E + U;
  { int b = (w->L == 0) ? 1 : w->b0; int tc = 0;
    if (!w->isroot && b == 1){ int mZ = (w->m3 + w->pw3) % 3; int QZ = (w->Q3 + (w->L ? (w->last + w->rev3) : 0)) % 3;
      tc = ((mZ + 1) % 3 == 0) && ((QZ + 1) % 3 == 0); }
    record(A, N, w->j3, w->j5, w->j7, w->isroot, b, tc); }
  if (w->j <= 1){ int nj3 = imin(vcap(D,3,6), vcap(U,3,6)); record(A, N, nj3, imin(vcap(D,5,2), vcap(U,5,2)), imin(vcap(D,7,2), vcap(U,7,2)), 0, 0, 0); }
  node c;
  if (D + 2*U <= X){ c = *w; c.E = 0; c.U = 2*U; c.b0 = (w->L == 0) ? 0 : w->b0; c.last = 0; c.rev3 = (2*w->rev3) % 3; c.L = w->L + 1; c.pw3 = (2*w->pw3) % 3; collect(A, &c, depth+1, DEP); }
  if (D - E + 3*U - 2*d <= X){ c = *w; c.D = D - E + U; c.E = d; c.U = 2*U - d;
    c.m3 = (w->m3 + w->pw3) % 3; c.Q3 = (w->Q3 + (w->L ? (w->last + w->rev3) : 0)) % 3;
    c.b0 = (w->L == 0) ? 1 : w->b0; c.last = 1; c.rev3 = (2*w->rev3 + 1) % 3; c.L = w->L + 1; c.pw3 = (2*w->pw3) % 3; collect(A, &c, depth+1, DEP); }
  if (w->j < 2 && D + 2*U <= X){ c.D = D - E + U; c.E = D - E; c.U = D + U; c.d = D; c.j = w->j + 1;
    c.j3 = imin(vcap(D,3,6), vcap(U,3,6)); c.j5 = imin(vcap(D,5,2), vcap(U,5,2)); c.j7 = imin(vcap(D,7,2), vcap(U,7,2));
    c.isroot = 0; c.L = 0; c.b0 = 0; c.last = 0; c.m3 = 0; c.Q3 = 0; c.rev3 = 0; c.pw3 = 1; collect(A, &c, depth+1, DEP); }
}

static int v3cap6(ll x){ if (x == 0) return 6; return vcap(x, 3, 6); }

int main(int argc, char **argv){
  if (argc < 5){ fprintf(stderr, "usage: census_lf_audit X LO HI nthreads [lowcap]\n"); return 1; }
  X = atoll(argv[1]); LO = atoll(argv[2]); HI = atoll(argv[3]); int nth = atoi(argv[4]); if (argc > 5) LOWCAP = atoi(argv[5]);
  if (HI > X + 1) HI = X + 1;
  omp_set_num_threads(nth);
  ACC = calloc(nth + 1, sizeof(acc_t));
  size_t NC = (size_t)(X >> 1) + 2, NW = (size_t)((HI - LO) >> 1) + 2;
  for (int t = 0; t <= nth; t++){ ACC[t].cnt = calloc(NC, sizeof(unsigned int)); ACC[t].ucnt = calloc(NW, sizeof(unsigned int));
    if (!ACC[t].cnt || !ACC[t].ucnt){ fprintf(stderr, "oom\n"); return 1; } }
  node root = { 1, 0, 2, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };
  int DEP = 12;
  collect(&ACC[nth], &root, 0, DEP);   /* shallow part, sequential, into the extra accumulator */
  #pragma omp parallel for schedule(dynamic, 64)
  for (ll i = 0; i < (ll)nT; i++) visit(&ACC[omp_get_thread_num()], &T[i]);
  /* reduce into ACC[nth] */
  acc_t *S = &ACC[nth];
  for (int t = 0; t < nth; t++){
    for (size_t i = 0; i < NC; i++) S->cnt[i] += ACC[t].cnt[i];
    for (size_t i = 0; i < NW; i++) S->ucnt[i] += ACC[t].ucnt[i];
    for (int a = 0; a < 7; a++) for (int b = 0; b < 2; b++) for (int r = 0; r < M729; r++) S->c3[a][b][r] += ACC[t].c3[a][b][r];
    for (int a = 0; a < 3; a++) for (int b = 0; b < 2; b++){ for (int r = 0; r < 5; r++) S->c5[a][b][r] += ACC[t].c5[a][b][r]; for (int r = 0; r < 7; r++) S->c7[a][b][r] += ACC[t].c7[a][b][r]; }
    for (int a = 0; a < 2; a++) for (int r = 0; r < M729; r++){ S->jz_tc[a][r] += ACC[t].jz_tc[a][r]; S->root[a][r] += ACC[t].root[a][r]; }
    for (int r = 0; r < M729; r++) S->jz_even[r] += ACC[t].jz_even[r];
    S->words += ACC[t].words;
  }
  S->cnt[0] += 1;   /* the empty word, value 1 */
  unsigned int *cnt = S->cnt;
  printf("census_lf_audit X=%lld window [%lld,%lld) threads %d tasks %zu counted words (excl. empty) %llu\n", X, LO, HI, nth, nT, S->words);
  /* misses below X */
  { ll nm = 0, last = 0; printf("misses (odd N <= X with M_2 = 0):"); for (ll v = 1; v <= X; v += 2) if (!cnt[v >> 1]){ nm++; last = v; printf(" %lld", v); }
    printf("\nmiss count %lld largest %lld\n", nm, last); }
  /* low values above 2^22 */
  { printf("low values (2^22 < N <= X, M_2 <= %d):", LOWCAP); ll nl = 0; for (ll v = (1LL << 22) + 1; v <= X; v += 2) if (cnt[v >> 1] <= (unsigned)LOWCAP){ nl++; printf(" %lld:%u", v, cnt[v >> 1]); }
    printf("\nlow count %lld\n", nl); }
  /* dyadic window means and minima */
  for (int k = 3; (1LL << (k+1)) <= X + 1; k++){ ll a = 1LL << k, b = 1LL << (k+1); double s = 0; ll n = 0; unsigned mn = ~0u, mx = 0;
    for (ll v = a | 1; v < b; v += 2){ s += cnt[v >> 1]; n++; if (cnt[v >> 1] < mn) mn = cnt[v >> 1]; if (cnt[v >> 1] > mx) mx = cnt[v >> 1]; }
    printf("dyadic [2^%d,2^%d): mean %.4f min %u max %u\n", k, k+1, s/n, mn, mx); }
  /* window statistics */
  double all = 0, nall = 0; unsigned maxm = 0;
  for (ll v = LO | 1; v < HI; v += 2){ all += cnt[v >> 1]; nall++; if (cnt[v >> 1] > maxm) maxm = cnt[v >> 1]; }
  printf("window odd N %.0f overall mean %.6f max M_2 %u\n", nall, all/nall, maxm);
  int mods[8] = {27, 5, 7, 11, 13, 17, 23, 31};
  for (int t = 0; t < 8; t++){ int q = mods[t]; double *s = calloc(q, sizeof(double)), *n = calloc(q, sizeof(double));
    for (ll v = LO | 1; v < HI; v += 2){ s[v % q] += cnt[v >> 1]; n[v % q]++; }
    printf("mod %d class means:", q); for (int r = 0; r < q; r++) printf(" %.6f", s[r]/n[r]); printf("\n"); free(s); free(n); }
  { double sv[7] = {0}, nv[7] = {0}; for (ll v = LO | 1; v < HI; v += 2){ int k = v3cap6(v); sv[k] += cnt[v >> 1]; nv[k]++; }
    printf("mean M_2 by v3(N)=0..6(+):"); for (int k = 0; k < 7; k++) printf(" %.4f", sv[k]/nv[k]); printf("\n"); }
  /* decomposition */
  static ll n729[M729]; ll n5[5] = {0}, n7[7] = {0};
  for (ll v = LO | 1; v < HI; v += 2){ n729[v % M729]++; n5[v % 5]++; n7[v % 7]++; }
  for (int a = 0; a < 7; a++) for (int b = 0; b < 2; b++){
    double byv[7] = {0}, nv[7] = {0}, by9[9] = {0}, n9[9] = {0}; ull tot = 0;
    for (int r = 0; r < M729; r++){ if (!n729[r]) continue; int k = v3cap6(r); byv[k] += S->c3[a][b][r]; nv[k] += n729[r]; by9[r % 9] += S->c3[a][b][r]; n9[r % 9] += n729[r]; tot += S->c3[a][b][r]; }
    printf("j3=%d b=%d total=%llu mean by v3(N)=0..6:", a, b, tot); for (int k = 0; k < 7; k++) printf(" %.3f", nv[k] ? byv[k]/nv[k] : 0.0);
    printf(" | by N mod 9:"); for (int r = 0; r < 9; r++) printf(" %.3f", n9[r] ? by9[r]/n9[r] : 0.0); printf("\n");
  }
  /* confinement exactness check: seeds with j3 = a >= 1 and odd tails must have N = 0 mod 3^a */
  { ull bad = 0; for (int a = 1; a < 7; a++){ ll mod = 1; for (int i = 0; i < a; i++) mod *= 3; for (int r = 0; r < M729; r++) if (r % mod) bad += S->c3[a][1][r]; }
    printf("odd-tail mass of confined seeds (j3=a>=1) outside N = 0 mod 3^a: %llu\n", bad); }
  { double mn = 1e300, mx = 0; int rmn = -1, rmx = -1; for (int r = 0; r < M729; r++){ if (!n729[r]) continue; double m = (double)(S->c3[0][0][r] + S->c3[0][1][r]) / n729[r]; if (m < mn){ mn = m; rmn = r; } if (m > mx){ mx = m; rmx = r; } }
    printf("unconfined (j3=0, root included) mean per class mod 729: min %.4f (r=%d) max %.4f (r=%d) ratio %.5f\n", mn, rmn, mx, rmx, mx/mn); }
  /* permutation null for the 729-class max/min ratio: random relabelling of the window's N into groups of the same sizes */
  { ll nw = (HI - (LO | 1) + 1) / 2; unsigned int *u = S->ucnt; ull seed = 88172645463325252ULL;
    double worst = 0, avg = 0; int R = 20;
    ll *perm = malloc(nw * sizeof(ll)); for (ll i = 0; i < nw; i++) perm[i] = i;
    for (int rep = 0; rep < R; rep++){
      for (ll i = nw - 1; i > 0; i--){ seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; ll k = (ll)(seed % (ull)(i + 1)); ll t = perm[i]; perm[i] = perm[k]; perm[k] = t; }
      double s[M729] = {0}; ll n[M729] = {0};
      for (ll i = 0; i < nw; i++){ ll v = (LO | 1) + 2*i; int r = (int)(v % M729); s[r] += u[perm[i]]; n[r]++; }
      double mn = 1e300, mx = 0; for (int r = 0; r < M729; r++){ if (!n[r]) continue; double m = s[r]/n[r]; if (m < mn) mn = m; if (m > mx) mx = m; }
      avg += mx/mn; if (mx/mn > worst) worst = mx/mn;
    }
    printf("permutation null (unconfined per-N counts shuffled across the window, %d reps): mean max/min %.5f, worst %.5f\n", R, avg/R, worst);
    free(perm); }
  /* tail-confinement test for J-seeds with j3 = 0 */
  { const char *nm[4] = {"J-seed j3=0 odd tail, tail-confined mod 3", "J-seed j3=0 odd tail, not tail-confined", "J-seed j3=0 even tail", "root seed (no J), both parities"};
    for (int g = 0; g < 4; g++){ double s3[3] = {0}, n3[3] = {0}; ull tot = 0;
      for (int r = 0; r < M729; r++){ ull x = g == 0 ? S->jz_tc[1][r] : g == 1 ? S->jz_tc[0][r] : g == 2 ? S->jz_even[r] : S->root[0][r] + S->root[1][r];
        s3[r % 3] += x; n3[r % 3] += n729[r]; tot += x; }
      printf("%s: total %llu, mean by N mod 3: %.4f %.4f %.4f, ratio class0/avg(nonzero) %.5f\n", nm[g], tot, s3[0]/n3[0], s3[1]/n3[1], s3[2]/n3[2], (s3[0]/n3[0]) / (0.5*(s3[1]/n3[1] + s3[2]/n3[2]))); }
  }
  /* by v3 total (decomp.c's TOTAL line) and p = 5, 7 lines */
  { double byv[7] = {0}, nv[7] = {0}; for (int r = 0; r < M729; r++){ if (!n729[r]) continue; int k = v3cap6(r); double s = 0; for (int a = 0; a < 7; a++) for (int b = 0; b < 2; b++) s += S->c3[a][b][r]; byv[k] += s; nv[k] += n729[r]; }
    printf("TOTAL mean by v3(N)=0..6:"); for (int k = 0; k < 7; k++) printf(" %.3f", byv[k]/nv[k]); printf("\n"); }
  for (int p = 5; p <= 7; p += 2){ printf("p=%d mean by N mod p (all seeds):", p);
    for (int r = 0; r < p; r++){ double s = 0; for (int a = 0; a < 3; a++) for (int b = 0; b < 2; b++) s += (p == 5 ? S->c5[a][b][r] : S->c7[a][b][r]); printf(" %.3f", s / (p == 5 ? n5[r] : n7[r])); }
    printf(" | unconfined only:"); for (int r = 0; r < p; r++){ double s = 0; for (int b = 0; b < 2; b++) s += (p == 5 ? S->c5[0][b][r] : S->c7[0][b][r]); printf(" %.3f", s / (p == 5 ? n5[r] : n7[r])); }
    printf("\n"); }
  return 0;
}
