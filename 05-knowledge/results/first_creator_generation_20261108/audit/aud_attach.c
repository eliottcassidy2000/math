/* aud_attach.c - independent auditor's check of PROOF.md section 1 (the attachment form).
 * Written from scratch; shares no code with ext.c.
 *
 * For a tournament T on V = {0..n-1} and a pattern x (bit a of x = [a -> v]) the claim is
 *   H(T + v_x) = H(T) + sum_a x_a L_a - sum_{a<b} x_a x_b M_ab,
 *   W[a][b] = sum_{S: a in S, b notin S} f[S][a] g[V\S][b],  L_a = f[V][a] - g[V][a] + sum_{b!=a} W[a][b],
 *   M_ab = W[a][b] + W[b][a]   (f = paths ending at, g = paths starting at).
 * Checks, for every tested T and EVERY pattern x:
 *   (1) the form equals a direct subset-DP count of Hamiltonian paths of the (n+1)-vertex tournament;
 *   (2) L_a = H(T + v_{e_a}) - H(T) >= 0;
 *   (3) H-neutral patterns (H(T+v_x) = H(T)) are exactly the x whose support {a : a -> v} is a union of INITIAL strong
 *       components (prefix of the condensation order), and there are #SCC + 1 of them;
 *   (4) for x with x_a = 1, x_b = 0: W[a][b] = #Ham paths of T + v_x with a, v, b consecutive (direct DP on a pruned
 *       digraph, independent of f/g);
 *   (5) the coefficient bounds W, L, M <= max H(n+1) (A038375), and the sharper M_ab <= H(T) + ... (recorded maxima).
 * usage: aud_attach SEED                                  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
typedef long long i64;
#define MX 12

static u64 rs = 88172645463325252ULL;
static u64 rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static const i64 maxH[13] = {1, 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095, 531757};   /* A038375, index = order */

static u64 dpbuf[1 << MX][MX];

/* Hamiltonian paths of a digraph on m vertices, B[i][j] = arc i->j.  mode 0: all; 1: starting at vertex s;
   2: ending at vertex s. */
static u64 hp_digraph(int m, unsigned char B[MX][MX], int mode, int s) {
  int F = 1 << m;
  for (int S = 0; S < F; S++) for (int v = 0; v < m; v++) dpbuf[S][v] = 0;
  for (int v = 0; v < m; v++) if (mode != 1 || v == s) dpbuf[1 << v][v] = 1;
  for (int S = 1; S < F; S++) for (int v = 0; v < m; v++) {
    if (!((S >> v) & 1)) continue;
    u64 x = dpbuf[S][v]; if (!x) continue;
    for (int w = 0; w < m; w++) if (!((S >> w) & 1) && B[v][w]) dpbuf[S | (1 << w)][w] += x;
  }
  u64 h = 0;
  for (int v = 0; v < m; v++) if (mode != 2 || v == s) h += dpbuf[F - 1][v];
  return h;
}

static i64 fe[1 << MX][MX], gs[1 << MX][MX];
static i64 W[MX][MX], L[MX], M[MX][MX];
static long nfail = 0, nchecks = 0, ntour = 0;
static i64 maxW_ratio_num = 0, maxM_seen[MX + 2], maxL_seen[MX + 2], maxW_seen[MX + 2];
static long neutral_hist[MX + 2];

static void form(int n, unsigned char A[MX][MX], i64 *H) {
  int F = 1 << n;
  for (int S = 0; S < F; S++) for (int a = 0; a < n; a++) { fe[S][a] = 0; gs[S][a] = 0; }
  for (int a = 0; a < n; a++) { fe[1 << a][a] = 1; gs[1 << a][a] = 1; }
  for (int S = 1; S < F; S++) for (int a = 0; a < n; a++) if ((S >> a) & 1) {
    if (fe[S][a]) for (int b = 0; b < n; b++) if (!((S >> b) & 1) && A[a][b]) fe[S | (1 << b)][b] += fe[S][a];
    if (gs[S][a]) for (int b = 0; b < n; b++) if (!((S >> b) & 1) && A[b][a]) gs[S | (1 << b)][b] += gs[S][a];
  }
  *H = 0; for (int a = 0; a < n; a++) *H += fe[F - 1][a];
  i64 Hg = 0; for (int a = 0; a < n; a++) Hg += gs[F - 1][a];
  if (Hg != *H) { printf("FAIL: sum f != sum g\n"); nfail++; }
  for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) W[a][b] = 0;
  for (int S = 1; S < F - 1; S++) for (int a = 0; a < n; a++) if ((S >> a) & 1) {
    if (!fe[S][a]) continue;
    int C = (F - 1) & ~S;
    for (int b = 0; b < n; b++) if ((C >> b) & 1) W[a][b] += fe[S][a] * gs[C][b];
  }
  for (int a = 0; a < n; a++) {
    L[a] = fe[F - 1][a] - gs[F - 1][a];
    for (int b = 0; b < n; b++) if (b != a) L[a] += W[a][b];
  }
  for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) M[a][b] = (a == b) ? 0 : W[a][b] + W[b][a];
}

static void build_ext(int n, unsigned char A[MX][MX], unsigned x, unsigned char B[MX][MX]) {
  for (int i = 0; i <= n; i++) for (int j = 0; j <= n; j++) B[i][j] = 0;
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) B[i][j] = A[i][j];
  for (int a = 0; a < n; a++) { if ((x >> a) & 1) B[a][n] = 1; else B[n][a] = 1; }
}

/* strong components in condensation order: comp[v] = index, 0 = top (beats all later) */
static int scc_order(int n, unsigned char A[MX][MX], int *comp) {
  unsigned char R[MX][MX];
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) R[i][j] = (i == j) || A[i][j];
  for (int k = 0; k < n; k++) for (int i = 0; i < n; i++) if (R[i][k]) for (int j = 0; j < n; j++) if (R[k][j]) R[i][j] = 1;
  /* rank of v = number of vertices u that reach v but are not reached by v (strictly above) */
  int rank[MX];
  for (int v = 0; v < n; v++) { rank[v] = 0; for (int u = 0; u < n; u++) if (R[u][v] && !R[v][u]) rank[v]++; }
  /* distinct ranks -> component index */
  int vals[MX], nv = 0;
  for (int v = 0; v < n; v++) { int seen = 0; for (int i = 0; i < nv; i++) if (vals[i] == rank[v]) seen = 1; if (!seen) vals[nv++] = rank[v]; }
  for (int i = 1; i < nv; i++) { int t = vals[i], j = i - 1; while (j >= 0 && vals[j] > t) { vals[j + 1] = vals[j]; j--; } vals[j + 1] = t; }
  for (int v = 0; v < n; v++) for (int i = 0; i < nv; i++) if (vals[i] == rank[v]) comp[v] = i;
  /* sanity: same component iff mutual reachability */
  for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if ((comp[u] == comp[v]) != (R[u][v] && R[v][u])) { printf("FAIL scc\n"); nfail++; }
  return nv;
}

static void test_tournament(int n, unsigned char A[MX][MX], int check_w) {
  i64 H; form(n, A, &H);
  i64 Hd = (i64)hp_digraph(n, A, 0, 0);
  if (H != Hd) { printf("FAIL H(T) n=%d %lld %lld\n", n, H, Hd); nfail++; }
  int comp[MX]; int k = scc_order(n, A, comp);
  unsigned prefix_ok[1 << MX]; memset(prefix_ok, 0, sizeof(unsigned) * (1u << n));
  for (int j = 0; j <= k; j++) { unsigned m = 0; for (int v = 0; v < n; v++) if (comp[v] < j) m |= 1u << v; prefix_ok[m] = 1; }
  unsigned char B[MX][MX];
  int nneutral = 0;
  for (unsigned x = 0; x < (1u << n); x++) {
    i64 hf = H;
    for (int a = 0; a < n; a++) if ((x >> a) & 1) hf += L[a];
    for (int a = 0; a < n; a++) if ((x >> a) & 1) for (int b = a + 1; b < n; b++) if ((x >> b) & 1) hf -= M[a][b];
    build_ext(n, A, x, B);
    i64 hd = (i64)hp_digraph(n + 1, B, 0, 0);
    nchecks++;
    if (hf != hd) { if (nfail < 20) printf("FAIL form n=%d x=%u form=%lld direct=%lld\n", n, x, hf, hd); nfail++; }
    int neutral = (hd == H);
    if (neutral) nneutral++;
    if (neutral != (int)prefix_ok[x]) { if (nfail < 20) printf("FAIL neutral n=%d x=%u\n", n, x); nfail++; }
    if (hd < H) { printf("FAIL insertion lemma\n"); nfail++; }
  }
  if (nneutral != k + 1) { printf("FAIL neutral count %d vs %d\n", nneutral, k + 1); nfail++; }
  if (n + 1 < 14) neutral_hist[k]++;
  for (int a = 0; a < n; a++) {
    if (L[a] < 0) { printf("FAIL L<0\n"); nfail++; }
    build_ext(n, A, 1u << a, B);
    i64 he = (i64)hp_digraph(n + 1, B, 0, 0);
    if (L[a] != he - H) { printf("FAIL L_a\n"); nfail++; }
    if (L[a] > maxL_seen[n]) maxL_seen[n] = L[a];
    for (int b = 0; b < n; b++) if (b != a) {
      if (W[a][b] > maxW_seen[n]) maxW_seen[n] = W[a][b];
      if (b > a && M[a][b] > maxM_seen[n]) maxM_seen[n] = M[a][b];
      if (W[a][b] > maxH[n + 1] || M[a][b] > maxH[n + 1] || L[a] > maxH[n + 1]) { printf("FAIL coefficient bound n=%d\n", n); nfail++; }
    }
  }
  if (check_w) {
    /* (4): W[a][b] = paths of T+v_x with a,v,b consecutive, for one random x with x_a=1, x_b=0, all ordered pairs */
    for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if (a != b) {
      unsigned x = (unsigned)(rnd() & ((1u << n) - 1)); x |= 1u << a; x &= ~(1u << b);
      build_ext(n, A, x, B);
      /* prune: only arc into v is a->v, only arc out of v is v->b */
      for (int c = 0; c < n; c++) { if (c != a) B[c][n] = 0; if (c != b) B[n][c] = 0; }
      i64 tot = (i64)hp_digraph(n + 1, B, 0, 0), st = (i64)hp_digraph(n + 1, B, 1, n), en = (i64)hp_digraph(n + 1, B, 2, n);
      if (tot - st - en != W[a][b]) { if (nfail < 20) printf("FAIL W interpretation n=%d a=%d b=%d %lld vs %lld\n", n, a, b, tot - st - en, W[a][b]); nfail++; }
    }
  }
  ntour++;
}

static void random_tournament(int n, unsigned char A[MX][MX], int bias) {
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] = 0;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
    int fwd;
    if (bias) fwd = (rnd() % 100) < (u64)bias; else fwd = rnd() & 1;
    if (fwd) A[i][j] = 1; else A[j][i] = 1;
  }
  /* random relabelling */
  int p[MX]; for (int i = 0; i < n; i++) p[i] = i;
  for (int i = n - 1; i > 0; i--) { int j = (int)(rnd() % (u64)(i + 1)); int t = p[i]; p[i] = p[j]; p[j] = t; }
  unsigned char C[MX][MX];
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) C[p[i]][p[j]] = A[i][j];
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] = C[i][j];
}

int main(int argc, char **argv) {
  if (argc > 1) rs ^= (u64)strtoull(argv[1], 0, 10) * 0x9E3779B97F4A7C15ULL;
  unsigned char A[MX][MX];
  /* exhaustive labelled n = 1..6 */
  for (int n = 1; n <= 6; n++) {
    int m = n * (n - 1) / 2;
    for (u64 c = 0; c < (1ULL << m); c++) {
      int b = 0;
      for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] = 0;
      for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if ((c >> b) & 1) A[i][j] = 1; else A[j][i] = 1; b++; }
      test_tournament(n, A, n <= 4);
    }
    printf("exhaustive labelled n=%d done: tournaments %ld, pattern checks %ld, failures %ld\n", n, ntour, nchecks, nfail);
    fflush(stdout);
  }
  /* random n = 7..10: uniform, and biased toward transitive (to get many strong components) */
  int counts[11] = {0, 0, 0, 0, 0, 0, 0, 600, 300, 120, 50};
  for (int n = 7; n <= 10; n++) {
    for (int t = 0; t < counts[n]; t++) {
      int bias = (t % 3 == 0) ? 0 : (t % 3 == 1 ? 80 : 93);
      random_tournament(n, A, bias);
      test_tournament(n, A, (t % 10) == 0);
    }
    printf("random n=%d done: tournaments %ld, pattern checks %ld, failures %ld\n", n, ntour, nchecks, nfail);
    fflush(stdout);
  }
  printf("max coefficients seen (n: W, L, M; max H(n+1)):\n");
  for (int n = 2; n <= 10; n++) printf("  n=%d: W %lld, L %lld, M %lld; maxH(n+1) %lld\n", n, maxW_seen[n], maxL_seen[n], maxM_seen[n], maxH[n + 1]);
  printf("strong-component histogram of tested T (k: count):"); for (int k = 0; k <= 11; k++) if (neutral_hist[k]) printf(" %d:%ld", k, neutral_hist[k]); printf("\n");
  printf("TOTAL: tournaments %ld, pattern checks %ld, failures %ld -> %s\n", ntour, nchecks, nfail, nfail ? "FAIL" : "PASS");
  return nfail ? 1 : 0;
}
