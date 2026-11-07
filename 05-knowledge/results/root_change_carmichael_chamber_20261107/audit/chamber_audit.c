/* chamber_audit.c -- independent auditor code for THM-4904 (1)-(5) (audit of kind-pasteur-2026-11-07-S14).
 * Written from the definitions only; it shares no code with rc_census.c, verify_chamber.c or literal_witnesses.py.
 *
 * Moves on s = (D,E,U,d), from t = (1,0,2,1):
 *   O = (D,0,2U,d),  Z = (D-E+U,d,2U-d,d),  J = (D-E+U,D-E,D+U,D).
 * Census definition (PROOF section 0): M_K(N) = #{u : u empty or u ends in Z or J, #J(u) <= K, D(t.u) = N}.
 * The DFS runs over prefixes w (any last letter) with F(t.w) <= X, F = D-E+U; the counted words are u = wZ
 * (needs #J(w) <= K) and u = wJ (needs #J(w)+1 <= K), and their values are read off as the D coordinate of the
 * child state (not from the formula F), so the twin decomposition is tested rather than assumed.
 *
 * Checks (all exact integer checks unless stated):
 *  [adj]   F(t.w) = F(t.w^dagger) for every word of length <= LADJ (any number of J), by carrying the row vector
 *          f*G_(w^dagger) (f = (1,-1,1,0)) and updating it on the right by sigma(letter);
 *  [Q]     Session199 eq. (3): O^T Q = Q Z, Z^T Q = Q O, J^T Q = Q J;
 *  [phys]  D-E >= d >= 1, D <= U, U >= 2d, E >= 0, d >= E, D-E+d <= U, D,d odd, E = U mod 2, on every visited state;
 *  [incr]  F(child) > F(parent) for each of the three children of every visited state (justifies pruning);
 *  [grow]  U >= 2 c0^(#J+1) 2^n (THM-4904 (4), long double) and the stronger W = U-d >= 2^n (Session200
 *          root-weighted-height section 1), n = number of letters O, Z;
 *  [T1]    M_K = [N=1] + 2 A_(K-1) + B_K for every odd N <= X;
 *  [T2]    Fix_K computed twice: (a) testing each enumerated prefix for w = w^dagger, (b) constructing every
 *          A c A^dagger directly from A; then M_K = [N=1] + Fix_K (mod 2);
 *  [T3]    every singleton N > 1 (M_K(N) = 1): the unique counted word is wZ with #J(w) = K, w = w^dagger,
 *          N not in F_<=K-1, w not beginning/ending with J (K >= 1), K = 2 shape B J C C^dagger J B^dagger, B != empty;
 *          extra (JZ = OJ identity): w contains neither the factor JZ nor OJ.
 * usage: chamber_audit X KMAX LADJ [windows_from_k]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

typedef long long ll;
typedef struct { ll D, E, U, d; } st;

static inline st mO(st s){ st r = { s.D, 0, 2*s.U, s.d }; return r; }
static inline st mZ(st s){ st r = { s.D - s.E + s.U, s.d, 2*s.U - s.d, s.d }; return r; }
static inline st mJ(st s){ st r = { s.D - s.E + s.U, s.D - s.E, s.D + s.U, s.D }; return r; }
static inline ll FF(st s){ return s.D - s.E + s.U; }
static inline st mv(st s, char c){ return c == 'O' ? mO(s) : (c == 'Z' ? mZ(s) : mJ(s)); }
static inline char sg(char c){ return c == 'O' ? 'Z' : (c == 'Z' ? 'O' : 'J'); }

static ll X; static int KMAX;
static uint32_t *Mdir[3];   /* direct census count of counted words by value, index N>>1 */
static uint32_t *Bex[3];    /* prefixes with exactly j letters J and F(t.w) = N */
static uint32_t *FixA[3];   /* dagger-fixed prefixes, exactly j letters J (test on enumerated prefixes) */
static uint32_t *FixB[3];   /* dagger-fixed prefixes built as A c A^dagger */
static ll nstates = 0, fail_phys = 0, fail_incr = 0, fail_grow = 0, fail_W = 0, fail_par = 0;
static long double c0, minratio_c0 = 1e30L, minratio_W = 1e30L;
static char wb[512];

static int dagger_fixed(const char *w, int n){
  for (int i = 0; i < n; i++) if (w[i] != sg(w[n-1-i])) return 0;
  return 1;
}

static void dfs(st s, int j, int len){
  ll F = FF(s);
  if (F > X) return;
  nstates++;
  /* physical inequalities (positive coordinates) */
  if (!(s.D - s.E >= s.d && s.d >= 1 && s.D <= s.U && s.U >= 2*s.d && s.E >= 0 && s.d >= s.E && s.D - s.E + s.d <= s.U))
    fail_phys++;
  if (!((s.D & 1) && (s.d & 1) && (((s.E - s.U) & 1) == 0))) fail_par++;
  /* growth */
  int n = len - j;
  { long double lb = 2.0L * powl(c0, (long double)(j + 1)) * ldexpl(1.0L, n);
    long double r = (long double)s.U / lb; if (r < minratio_c0) minratio_c0 = r;
    if (r < 1.0L - 1e-15L) fail_grow++;
    ll W = s.U - s.d;
    if (n >= 62 || W < (1LL << n)) fail_W++;
    long double rw = (long double)W / ldexpl(1.0L, n); if (rw < minratio_W) minratio_W = rw; }
  st cO = mO(s), cZ = mZ(s), cJ = mJ(s);
  if (!(FF(cO) > F && FF(cZ) > F && FF(cJ) > F)) fail_incr++;
  Bex[j][F >> 1]++;
  if (dagger_fixed(wb, len)) FixA[j][F >> 1]++;
  for (int K = 0; K <= KMAX; K++){
    if (j <= K)     { ll v = cZ.D; if (v <= X) Mdir[K][v >> 1]++; }
    if (j + 1 <= K) { ll v = cJ.D; if (v <= X) Mdir[K][v >> 1]++; }
  }
  if (len >= 500) { fprintf(stderr, "word too long\n"); exit(2); }
  wb[len] = 'O'; dfs(cO, j, len + 1);
  wb[len] = 'Z'; dfs(cZ, j, len + 1);
  if (j < KMAX){ wb[len] = 'J'; dfs(cJ, j + 1, len + 1); }
}

/* chamber words A c A^dagger, #J(A) = floor(K/2), c = J iff K odd; prune on F(t.A) > X and 4^n(A) >= X
   (W >= 2^n gives F(t.A c A^dagger) > 2^(2 n(A))) */
static char ab[256]; static ll chamber_words[3];
static void chamber(st s, int jA, int len){
  ll F = FF(s);
  if (F > X) return;
  int nA = len - jA;
  if (2*nA < 62 && (1LL << (2*nA)) >= X) return;
  for (int K = 0; K <= KMAX; K++){
    if (K/2 != jA) continue;
    st q = s; int L = len;
    if (K & 1) q = mJ(q);
    for (int i = len - 1; i >= 0; i--) q = mv(q, sg(ab[i]));
    (void)L;
    ll v = FF(q);
    if (v <= X){ FixB[K][v >> 1]++; chamber_words[K]++; }
  }
  if (len >= 250) return;
  ab[len] = 'O'; chamber(mO(s), jA, len + 1);
  ab[len] = 'Z'; chamber(mZ(s), jA, len + 1);
  if (jA < KMAX/2){ ab[len] = 'J'; chamber(mJ(s), jA + 1, len + 1); }
}

/* second pass: collect the counted words of singleton values */
typedef struct { ll N; int K; char w[200]; } sing_t;
static sing_t *SL; static ll nSL = 0, capSL = 0;
static void dfs2(st s, int j, int len){
  ll F = FF(s);
  if (F > X) return;
  st cO = mO(s), cZ = mZ(s), cJ = mJ(s);
  for (int K = 0; K <= KMAX; K++){
    for (int e = 0; e < 2; e++){
      if (e == 0 && !(j <= K)) continue;
      if (e == 1 && !(j + 1 <= K)) continue;
      ll v = (e == 0) ? cZ.D : cJ.D;
      if (v > 1 && v <= X && Mdir[K][v >> 1] == 1){
        if (nSL == capSL){ capSL = capSL ? 2*capSL : 4096; SL = realloc(SL, capSL * sizeof(sing_t)); }
        SL[nSL].N = v; SL[nSL].K = K; memcpy(SL[nSL].w, wb, len); SL[nSL].w[len] = e ? 'J' : 'Z'; SL[nSL].w[len+1] = 0;
        nSL++;
      }
    }
  }
  wb[len] = 'O'; dfs2(cO, j, len + 1);
  wb[len] = 'Z'; dfs2(cZ, j, len + 1);
  if (j < KMAX){ wb[len] = 'J'; dfs2(cJ, j + 1, len + 1); }
}
static int cmp_sing(const void *a, const void *b){
  const sing_t *x = a, *y = b;
  if (x->K != y->K) return x->K - y->K;
  return (x->N > y->N) - (x->N < y->N);
}

/* adjoint check by right multiplication of a row vector */
static const ll MO[4][4] = {{1,0,0,0},{0,0,0,0},{0,0,2,0},{0,0,0,1}};
static const ll MZ[4][4] = {{1,-1,1,0},{0,0,0,1},{0,0,2,-1},{0,0,0,1}};
static const ll MJ[4][4] = {{1,-1,1,0},{1,-1,0,0},{1,0,1,0},{1,0,0,0}};
static ll adj_n = 0, adj_fail = 0; static int LADJ;
static void adjdfs(st s, const ll r[4], int len){
  if (len > 0){
    adj_n++;
    ll lhs = FF(s);
    ll rhs = r[0]*1 + r[1]*0 + r[2]*2 + r[3]*1;     /* r . t */
    if (lhs != rhs) adj_fail++;
  }
  if (len == LADJ) return;
  const char *al = "OZJ";
  for (int a = 0; a < 3; a++){
    char c = al[a], cs = sg(c);
    const ll (*M)[4] = cs == 'O' ? MO : (cs == 'Z' ? MZ : MJ);
    ll r2[4];
    for (int k = 0; k < 4; k++){ r2[k] = 0; for (int i = 0; i < 4; i++) r2[k] += r[i] * M[i][k]; }
    adjdfs(mv(s, c), r2, len + 1);
  }
}
static int checkQ(void){
  static const ll Q[4][4] = {{0,0,0,1},{0,1,0,-1},{0,0,1,-1},{1,-1,-1,1}};
  const ll (*L[3])[4] = {MO, MZ, MJ}; const ll (*R[3])[4] = {MZ, MO, MJ};
  int ok = 1;
  for (int k = 0; k < 3; k++)
    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++){
      ll x = 0, y = 0;
      for (int i = 0; i < 4; i++){ x += L[k][i][a] * Q[i][b]; y += Q[a][i] * R[k][i][b]; }
      if (x != y) ok = 0;
    }
  /* the matrices must reproduce the moves */
  st s = {7, 3, 11, 5}; st o = mO(s), z = mZ(s), jj = mJ(s); ll v[4] = {7,3,11,5};
  const ll (*Ms[3])[4] = {MO, MZ, MJ}; st outs[3] = {o, z, jj};
  for (int k = 0; k < 3; k++){ ll w[4]; for (int a = 0; a < 4; a++){ w[a] = 0; for (int i = 0; i < 4; i++) w[a] += Ms[k][a][i]*v[i]; }
    if (w[0] != outs[k].D || w[1] != outs[k].E || w[2] != outs[k].U || w[3] != outs[k].d) ok = 0; }
  return ok;
}

int main(int argc, char **argv){
  if (argc < 4){ fprintf(stderr, "usage: chamber_audit X KMAX LADJ\n"); return 1; }
  X = atoll(argv[1]); KMAX = atoi(argv[2]); LADJ = atoi(argv[3]);
  if (KMAX < 0 || KMAX > 2){ fprintf(stderr, "KMAX in 0..2\n"); return 1; }
  c0 = 1.0L; for (int i = 0; i < 4000; i++) c0 *= 1.0L - 0.25L * powl(2.0L/3.0L, (long double)i);
  printf("c0 = prod_(i>=0) (1 - (1/4)(2/3)^i) = %.12Lf\n", c0);
  printf("[Q] Session199 identities O^T Q = Q Z, Z^T Q = Q O, J^T Q = Q J and matrix/move agreement: %s\n", checkQ() ? "OK" : "FAIL");
  if (LADJ > 0){ st t = {1,0,2,1}; ll f[4] = {1,-1,1,0}; adjdfs(t, f, 0);
    printf("[adj] F(t.w) = F(t.w^dagger) on all %lld words of length 1..%d (any number of J): failures %lld\n", adj_n, LADJ, adj_fail); }
  ll H = (X >> 1) + 2;
  for (int K = 0; K <= KMAX; K++){
    Mdir[K] = calloc(H, 4); Bex[K] = calloc(H, 4); FixA[K] = calloc(H, 4); FixB[K] = calloc(H, 4);
    if (!Mdir[K] || !Bex[K] || !FixA[K] || !FixB[K]){ fprintf(stderr, "oom\n"); return 1; }
    Mdir[K][0]++;   /* the empty word, value D(t) = 1 */
  }
  st t = {1, 0, 2, 1};
  dfs(t, 0, 0);
  printf("visited prefixes w with #J(w) <= %d and F(t.w) <= %lld: %lld\n", KMAX, X, nstates);
  printf("[phys] failures %lld; [parity of coordinates] failures %lld; [incr] failures %lld\n", fail_phys, fail_par, fail_incr);
  printf("[grow] THM-4904 (4) U >= 2 c0^(#J+1) 2^n: failures %lld (min ratio %.6Lf); Session200 W=U-d >= 2^n: failures %lld (min ratio %.6Lf)\n",
         fail_grow, minratio_c0, fail_W, minratio_W);
  chamber(t, 0, 0);
  for (int K = 0; K <= KMAX; K++){
    ll t1 = 0, fixdiff = 0, par = 0, sing = 0, nfix = 0, oddM = 0, oddM15 = 0;
    for (ll i = 0; i < (X + 1) / 2; i++){
      ll N = 2*i + 1;
      ll A = 0; for (int jj = 0; jj <= K - 1; jj++) A += Bex[jj][i];
      ll pred = (N == 1) + 2*A + Bex[K][i];
      if (pred != Mdir[K][i]) t1++;
      if (FixA[K][i] != FixB[K][i]) fixdiff++;
      nfix += FixB[K][i];
      if (((Mdir[K][i] - (N == 1) - FixB[K][i]) & 1) != 0) par++;
      if (N > 1 && Mdir[K][i] == 1) sing++;
      if (N >= 3 && (Mdir[K][i] & 1)){ oddM++; if (N > 32768) oddM15++; }
    }
    printf("K=%d: [T1] mismatches %lld; [T2] Fix(test) vs Fix(construct) mismatches %lld; parity-law failures %lld;"
           " singletons N>1: %lld; chamber words with value <= X: %lld (%.4f sqrt X); odd M_K at N>=3: %lld (above 2^15: %lld)\n",
           K, t1, fixdiff, par, sing, chamber_words[K], nfix / sqrt((double)X), oddM, oddM15);
  }
  /* singletons and their words */
  dfs2(t, 0, 0);
  qsort(SL, nSL, sizeof(sing_t), cmp_sing);
  for (int K = 0; K <= KMAX; K++){
    ll bad = 0, cnt = 0, badJZ = 0; ll prevN = -1; int dupl = 0;
    printf("K=%d singleton values N>1 (first 60) with the counted word u = wZ:", K);
    int shown = 0;
    for (ll q = 0; q < nSL; q++){
      if (SL[q].K != K) continue;
      if (SL[q].N == prevN){ dupl++; continue; }
      prevN = SL[q].N; cnt++;
      const char *u = SL[q].w; int L = (int)strlen(u); int n = L - 1;  /* prefix w = u[0..n-1] */
      int nj = 0; for (int i = 0; i < n; i++) nj += (u[i] == 'J');
      ll i2 = SL[q].N >> 1; ll Aprev = 0; for (int jj = 0; jj <= K - 1; jj++) Aprev += Bex[jj][i2];
      int ok = (u[L-1] == 'Z') && nj == K && dagger_fixed(u, n) && Aprev == 0 && Bex[K][i2] == 1;
      if (K >= 1) ok = ok && n > 0 && u[0] != 'J' && u[n-1] != 'J';
      if (K == 2){ /* B J C C^dagger J B^dagger with B nonempty: first J at position p >= 1 */
        int p = -1; for (int i = 0; i < n; i++) if (u[i] == 'J'){ p = i; break; }
        ok = ok && p >= 1; }
      if (!ok) bad++;
      for (int i = 0; i + 1 < n; i++) if ((u[i] == 'J' && u[i+1] == 'Z') || (u[i] == 'O' && u[i+1] == 'J')) { badJZ++; break; }
      if (shown < 60){ printf(" %lld:%s", SL[q].N, u); shown++; }
    }
    printf("\nK=%d singletons %lld; [T3] failures %lld; duplicate words for a singleton value %d; singleton prefixes containing JZ or OJ: %lld\n",
           K, cnt, bad, dupl, badJZ);
  }
  /* window tables (format comparable to VERIFY_CHAMBER_2p20.txt and the census tables) */
  for (int K = 0; K <= KMAX; K++){
    for (ll lo = 1; lo < X; lo *= 2){
      ll hi = lo * 2 > X ? X + 1 : lo * 2; ll c[9] = {0}; ll n = 0, mn = 1LL << 40, sat = 0; double s = 0, sx = 0;
      for (ll N = lo | 1; N < hi; N += 2){ ll m = Mdir[K][N >> 1]; ll ms = m > 255 ? 255 : m; n++; s += ms; sx += m; if (ms < mn) mn = ms; if (ms <= 8) c[ms]++; if (ms == 255) sat++; }
      printf("  K=%d [%lld,%lld): n=%lld min=%lld mean255=%.2f exactmean=%.2f sat255=%lld c0..8=", K, lo, hi, n, mn, s / n, sx / n, sat);
      for (int m = 0; m <= 8; m++) printf("%s%lld", m ? "," : "", c[m]);
      printf("\n");
    }
  }
  int pass = !fail_phys && !fail_par && !fail_incr && !fail_grow && !fail_W && !adj_fail && checkQ();
  printf("STRUCTURAL RESULT %s\n", pass ? "PASS" : "FAIL");
  return 0;
}
