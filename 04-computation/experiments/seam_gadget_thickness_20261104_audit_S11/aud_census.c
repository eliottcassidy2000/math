/* aud_census.c -- audit S11: census of SEAM transfer matrices of marked gadgets (G, u, z) by exact order.
   Independent of the audited gadget_matrices.c:
     * classes of order <= 7 by vertex extension + canonical form = minimum code over ALL n! relabellings (no pruning);
     * order 8 WITHOUT canonical forms: every order-8 tournament is isomorphic to an extension of an order-7 class
       representative by one vertex (128 patterns), so the 58,368 extensions carry every order-8 matrix;
     * matrices by route II = the (A) row statistics counted by subset DP on G alone (for every order <= 8), and for
       orders <= 7 also by route I = prefix fit with the prefixes singleton (1,1,1), TT_2 at its sink (1,0,2), TT_3 at
       its sink (1,0,4) (different prefixes from the audited code), counts by subset DP on the seamed tournament;
       route I vs route II is compared on every (class, u, z) of order <= 7 and on a random sample of order 8.
   Output: per-order distinct matrices to OUTDIR/aud_matsN.txt (9 ints per line, row-major) and a summary on stdout.
   Usage: aud_census NMAX OUTDIR                                                                                  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll; typedef unsigned long long ull;
#define MV 12
static unsigned long long rng_s = 0x9E3779B97F4A7C15ULL;
static unsigned long long rnd(void){ rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }
/* ---------------- canonical forms ---------------- */
static int nperm; static unsigned char (*perms)[8];
static void genperms(int n){ int idx[8]; for (int i = 0; i < n; i++) idx[i] = i; nperm = 0;
  while (1) { for (int i = 0; i < n; i++) perms[nperm][i] = idx[i]; nperm++;
    int i = n - 2; while (i >= 0 && idx[i] > idx[i + 1]) i--; if (i < 0) break;
    int j = n - 1; while (idx[j] < idx[i]) j--; int t = idx[i]; idx[i] = idx[j]; idx[j] = t;
    for (int l = i + 1, r = n - 1; l < r; l++, r--) { t = idx[l]; idx[l] = idx[r]; idx[r] = t; } } }
static ull encode(int n, int A[MV][MV]){ ull m = 0; int b = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (A[i][j]) m |= 1ULL << b; b++; } return m; }
static void decode(int n, ull m, int A[MV][MV]){ int b = 0; for (int i = 0; i < n; i++) { A[i][i] = 0; for (int j = i + 1; j < n; j++) { int e = (m >> b) & 1; b++; A[i][j] = e; A[j][i] = 1 - e; } } }
static ull canon(int n, int A[MV][MV]){ ull best = ~0ULL;
  for (int q = 0; q < nperm; q++) { unsigned char *p = perms[q]; ull m = 0; int b = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (A[p[i]][p[j]]) m |= 1ULL << b; b++; }
    if (m < best) best = m; }
  return best; }
/* ---------------- hash set of ull (class codes) ---------------- */
typedef struct { ull *t; size_t cap, cnt; } uset;
static int uset_ins(uset *s, ull k){ k += 1; size_t i = (size_t)((k * 0x9E3779B97F4A7C15ULL) >> 17) & (s->cap - 1);
  while (s->t[i]) { if (s->t[i] == k) return 0; i = (i + 1) & (s->cap - 1); } s->t[i] = k; s->cnt++; return 1; }
/* ---------------- hash set of matrices ---------------- */
typedef struct { ll (*m)[9]; size_t n, mcap; int *tab; size_t cap; } mset;
static ull mhash(const ll *a){ ull h = 1469598103934665603ULL; for (int i = 0; i < 9; i++) { h ^= (ull)a[i] + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2); h *= 1099511628211ULL; } return h; }
static void mset_init(mset *s, size_t cap){ s->cap = cap; s->tab = malloc(sizeof(int) * cap); memset(s->tab, -1, sizeof(int) * cap); s->n = 0; s->mcap = cap / 2; s->m = malloc(sizeof(ll) * 9 * s->mcap); }
static int mset_ins(mset *s, const ll *a){ size_t i = mhash(a) & (s->cap - 1);
  while (s->tab[i] >= 0) { if (!memcmp(s->m[s->tab[i]], a, sizeof(ll) * 9)) return 0; i = (i + 1) & (s->cap - 1); }
  if (s->n >= s->mcap) { fprintf(stderr, "mset full\n"); exit(1); }
  memcpy(s->m[s->n], a, sizeof(ll) * 9); s->tab[i] = (int)s->n; s->n++; return 1; }
static int mset_has(mset *s, const ll *a){ size_t i = mhash(a) & (s->cap - 1);
  while (s->tab[i] >= 0) { if (!memcmp(s->m[s->tab[i]], a, sizeof(ll) * 9)) return 1; i = (i + 1) & (s->cap - 1); } return 0; }
/* ---------------- route II: (A) statistics by subset DP on G ---------------- */
static int G[MV][MV];
static ll f[1 << 8][8], gs[1 << 8][8], Hm[1 << 8], px[8][1 << 8][8], C2[1 << 8];
static void gadget_tables(int n){
  int full = (1 << n) - 1;
  memset(f, 0, sizeof f); memset(gs, 0, sizeof gs); memset(px, 0, sizeof px);
  for (int v = 0; v < n; v++) { f[1 << v][v] = 1; gs[1 << v][v] = 1; px[v][1 << v][v] = 1; }
  for (int m = 1; m <= full; m++) for (int v = 0; v < n; v++) {
    for (int w = 0; w < n; w++) if (!(m >> w & 1)) {
      if (G[v][w]) { f[m | 1 << w][w] += f[m][v]; for (int x = 0; x < n; x++) px[x][m | 1 << w][w] += px[x][m][v]; }
      if (G[w][v]) gs[m | 1 << w][w] += gs[m][v]; } }
  Hm[0] = 1; for (int m = 1; m <= full; m++) { ll s = 0; for (int v = 0; v < n; v++) s += f[m][v]; Hm[m] = s; }
  for (int Yy = 0; Yy <= full; Yy++) { ll s = 0; for (int U = Yy; ; U = (U - 1) & Yy) { s += Hm[U] * Hm[Yy ^ U]; if (!U) break; } C2[Yy] = s; }
}
static void route2(int n, int u, int z, ll M[9]){
  int full = (1 << n) - 1;
  ll h = Hm[full], Eu = f[full][u], Lu = 0, sz = gs[full][z], Pzu = px[z][full][u], Lzu = 0, Rz = 0, Mzu = 0, K1 = 0, K2 = 0;
  for (int m = 1; m <= full; m++) {
    if (m >> u & 1) Lu += gs[m][u] * Hm[full ^ m];
    if (m >> z & 1) Rz += f[m][z] * Hm[full ^ m];
    if ((m >> z & 1) && !(m >> u & 1)) { Lzu += gs[m][z] * gs[full ^ m][u]; Mzu += f[m][z] * f[full ^ m][u]; }
    if ((m >> u & 1) && (m >> z & 1)) K1 += px[u][m][z] * C2[full ^ m];
  }
  if (u != z) for (int A = 1; A <= full; A++) if ((A >> z & 1) && !(A >> u & 1)) { ll fa = f[A][z]; if (!fa) continue; int rest = full ^ A;
      for (int B = rest; B; B = (B - 1) & rest) if (B >> u & 1) K2 += fa * gs[B][u] * Hm[rest ^ B]; }
  M[0] = h; M[1] = -Eu; M[2] = Lu; M[3] = sz; M[4] = -Pzu; M[5] = Lzu; M[6] = Rz; M[7] = -Mzu; M[8] = K1 + K2;
}
/* ---------------- route I: prefix fit by subset DP on the seamed tournament ---------------- */
static int Y[MV][MV]; static ll dp[1 << MV][MV], ds[1 << MV][MV], Hy[1 << MV];
static void state_dp(int n, int a, ll out[3]){
  int full = (1 << n) - 1;
  for (int m = 0; m <= full; m++) for (int v = 0; v < n; v++) { dp[m][v] = 0; ds[m][v] = 0; }
  for (int v = 0; v < n; v++) dp[1 << v][v] = 1; ds[1 << a][a] = 1;
  for (int m = 1; m <= full; m++) for (int v = 0; v < n; v++) { ll x = dp[m][v], y = ds[m][v]; if (!x && !y) continue;
    for (int w = 0; w < n; w++) if (!(m >> w & 1) && Y[v][w]) { dp[m | 1 << w][w] += x; ds[m | 1 << w][w] += y; } }
  Hy[0] = 1; for (int m = 1; m <= full; m++) { ll s = 0; for (int v = 0; v < n; v++) s += dp[m][v]; Hy[m] = s; }
  ll R = 0; for (int m = 1; m <= full; m++) if (m >> a & 1) R += dp[m][a] * Hy[full ^ m];
  ll S = 0; for (int v = 0; v < n; v++) S += ds[full][v];
  out[0] = Hy[full]; out[1] = S; out[2] = R;
}
static void route1(int gn, int u, int z, ll M[9]){
  static const int pn_[3] = {1, 2, 3}; static const ll V[3][3] = {{1, 1, 1}, {1, 0, 0}, {1, 2, 4}}; /* columns (1,1,1),(1,0,2),(1,0,4) */
  ll adjV[3][3]; ll dV = V[0][0] * (V[1][1] * V[2][2] - V[1][2] * V[2][1]) - V[0][1] * (V[1][0] * V[2][2] - V[1][2] * V[2][0]) + V[0][2] * (V[1][0] * V[2][1] - V[1][1] * V[2][0]);
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) { int r1 = (j + 1) % 3, r2 = (j + 2) % 3, c1 = (i + 1) % 3, c2 = (i + 2) % 3; adjV[i][j] = V[r1][c1] * V[r2][c2] - V[r1][c2] * V[r2][c1]; }
  ll Vo[3][3];
  for (int p = 0; p < 3; p++) { int pn = pn_[p], n = pn + gn, a = pn - 1;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Y[i][j] = 0;
    for (int i = 0; i < pn; i++) for (int j = i + 1; j < pn; j++) Y[i][j] = 1;            /* transitive prefix, exit = its sink */
    for (int i = 0; i < gn; i++) for (int j = 0; j < gn; j++) Y[pn + i][pn + j] = G[i][j];
    for (int i = 0; i < gn; i++) for (int j = 0; j < pn; j++) { Y[pn + i][j] = 1; Y[j][pn + i] = 0; }
    Y[pn + u][a] = 0; Y[a][pn + u] = 1;
    ll o[3]; state_dp(n, pn + z, o); for (int r = 0; r < 3; r++) Vo[r][p] = o[r]; }
  for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) { ll s = 0; for (int t = 0; t < 3; t++) s += Vo[r][t] * adjV[t][c];
    if (s % dV) { fprintf(stderr, "non-integer fit\n"); exit(1); } M[3 * r + c] = s / dV; }
}
static int cmp9(const void *a, const void *b){ const ll *x = a, *y = b; for (int i = 0; i < 9; i++) { if (x[i] < y[i]) return -1; if (x[i] > y[i]) return 1; } return 0; }
int main(int argc, char **argv){
  int NMAX = argc > 1 ? atoi(argv[1]) : 7; const char *outdir = argc > 2 ? argv[2] : ".";
  perms = malloc(sizeof(*perms) * 40320);
  ull *cls[9]; int ncls[9]; cls[1] = malloc(sizeof(ull)); cls[1][0] = 0; ncls[1] = 1;
  int maxc = NMAX < 7 ? NMAX : 7;
  for (int n = 2; n <= maxc; n++) { genperms(n); uset s; s.cap = 1 << 16; s.t = calloc(s.cap, sizeof(ull)); s.cnt = 0;
    cls[n] = malloc(sizeof(ull) * 4096); ncls[n] = 0;
    for (int c = 0; c < ncls[n - 1]; c++) { int A[MV][MV]; decode(n - 1, cls[n - 1][c], A);
      for (int pat = 0; pat < (1 << (n - 1)); pat++) { for (int i = 0; i < n - 1; i++) { A[i][n - 1] = (pat >> i) & 1; A[n - 1][i] = 1 - A[i][n - 1]; } A[n - 1][n - 1] = 0;
        ull cf = canon(n, A); if (uset_ins(&s, cf)) cls[n][ncls[n]++] = cf; } }
    free(s.t); }
  printf("isomorphism classes by order:"); for (int n = 1; n <= maxc; n++) printf(" %d", ncls[n]); printf("  (A000568: 1 1 2 4 12 56 456)\n");
  mset per[9]; for (int n = 1; n <= NMAX; n++) mset_init(&per[n], n <= 6 ? (1 << 14) : (n == 7 ? (1 << 17) : (1 << 21)));
  ll route_cmp = 0, route_bad = 0;
  for (int n = 1; n <= NMAX; n++) {
    ll pairs = 0;
    if (n <= 7) {
      for (int c = 0; c < ncls[n]; c++) { decode(n, cls[n][c], G); gadget_tables(n);
        for (int u = 0; u < n; u++) for (int z = 0; z < n; z++) { ll M2[9], M1[9]; route2(n, u, z, M2); route1(n, u, z, M1); route_cmp++;
          if (memcmp(M1, M2, sizeof M1)) { route_bad++; if (route_bad < 5) printf("ROUTE MISMATCH n=%d class %d u=%d z=%d\n", n, c, u, z); }
          mset_ins(&per[n], M2); pairs++; } }
    } else { /* n == 8: extensions of the order-7 class representatives */
      ll samp = 0;
      for (int c = 0; c < ncls[7]; c++) { int A[MV][MV]; decode(7, cls[7][c], A);
        for (int pat = 0; pat < 128; pat++) {
          for (int i = 0; i < 7; i++) for (int j = 0; j < 7; j++) G[i][j] = A[i][j];
          for (int i = 0; i < 7; i++) { G[i][7] = (pat >> i) & 1; G[7][i] = 1 - G[i][7]; } G[7][7] = 0;
          gadget_tables(8);
          for (int u = 0; u < 8; u++) for (int z = 0; z < 8; z++) { ll M2[9]; route2(8, u, z, M2); mset_ins(&per[8], M2); pairs++;
            if ((rnd() % 4000) == 0) { ll M1[9]; route1(8, u, z, M1); samp++; route_cmp++; if (memcmp(M1, M2, sizeof M1)) { route_bad++; printf("ROUTE MISMATCH order 8\n"); } } } } }
      printf("order 8: route I vs II checked on a random sample of %lld (G,u,z)\n", samp);
    }
    printf("order %d: (G,u,z) triples %lld, distinct matrices %zu\n", n, pairs, per[n].n); fflush(stdout);
  }
  printf("route I (prefix fit) vs route II ((A) statistics): %lld comparisons, %lld mismatches\n", route_cmp, route_bad);
  /* write per-order files */
  for (int n = 1; n <= NMAX; n++) { qsort(per[n].m, per[n].n, sizeof(ll) * 9, cmp9); char fn[512]; snprintf(fn, sizeof fn, "%s/aud_mats%d.txt", outdir, n);
    FILE *o = fopen(fn, "w"); for (size_t i = 0; i < per[n].n; i++) { ll *m = per[n].m[i]; fprintf(o, "%lld %lld %lld %lld %lld %lld %lld %lld %lld\n", m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]); } fclose(o); }
  /* cumulative statistics, parity checks (E), (F3), P5 */
  mset cum; mset_init(&cum, 1 << 21); mset r2s; mset_init(&r2s, 1 << 21); mset s0; mset_init(&s0, 1 << 21);
  printf("k | per-order distinct | cumulative distinct | S-row-zero per order | S-row-zero cumulative | distinct nonsingular 2x2 restrictions cumulative | per-order S-row-zero: R_z odd, M odd, det2 odd, det2!=M mod 2 | all per-order matrices: det2 != M mod 2 (no S-row hypothesis)\n");
  for (int n = 1; n <= NMAX; n++) {
    ll s0n = 0, rodd = 0, modd = 0, dodd = 0, dneq = 0, gen_neq = 0;
    for (size_t i = 0; i < per[n].n; i++) { ll *m = per[n].m[i]; mset_ins(&cum, m);
      ll d2 = m[0] * m[8] - m[2] * m[6];
      if (((d2 - (-m[7])) & 1)) gen_neq++;
      if (m[3] == 0 && m[4] == 0 && m[5] == 0) { s0n++; mset_ins(&s0, m); if (m[6] & 1) rodd++; if (m[7] & 1) modd++; if (d2 & 1) dodd++; if ((d2 - (-m[7])) & 1) dneq++;
        if (d2 != 0) { ll r[9] = {m[0], m[2], m[6], m[8], 0, 0, 0, 0, 0}; mset_ins(&r2s, r); } } }
    printf("%d | %zu | %zu | %lld | %zu | %zu | %lld, %lld, %lld, %lld | %lld\n", n, per[n].n, cum.n, s0n, s0.n, r2s.n, rodd, modd, dodd, dneq, gen_neq);
  }
  /* determinant statistics at orders 3 and 5 */
  for (int n = 3; n <= NMAX && n <= 5; n += 2) { ll dmin = 1LL << 60, dmax = -(1LL << 60), zero = 0;
    for (size_t i = 0; i < per[n].n; i++) { ll *m = per[n].m[i]; ll d = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
      if (!d) zero++; if (d < dmin) dmin = d; if (d > dmax) dmax = d; if (n == 3) printf("  order-3 matrix [%lld %lld %lld; %lld %lld %lld; %lld %lld %lld] det %lld\n", m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], d); }
    printf("order %d: %zu matrices, det 0 for %lld, det range [%lld, %lld]\n", n, per[n].n, zero, dmin, dmax); }
  return 0;
}
