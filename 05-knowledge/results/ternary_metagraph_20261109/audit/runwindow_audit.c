/* runwindow_audit.c - follow-up audit of THM-4951 (3f) (run-window theorem) and of the n = 3k remark in (3f).
 * Written by the auditor, independently of run_window_theorem_check.py / run_window_check.py.
 *
 * Frame: path 0 -> 1 -> ... -> n-1; tile (i, j), j >= i + 2, reversed (j -> i) when marked; t + 1 = complement;
 * positions a = 0..L-1 (L = n - 2), triple word x_a = [tile (a, a + 2) reversed]; sigma(i, j) = (n-1-j, n-1-i).
 * t_w marks exactly the tiles inside the windows [p, q + 2] of the 1-runs [p, q] of w (union of staircases).
 *
 * Mode "theorem": for every colouring w of [0, L-1] whose maximal runs all have length >= 2 (exhaustive for
 * n <= NEXH, uniformly sampled run compositions for larger n up to 32) check
 *   (i) the triple word of t_w is w;  (ii) nu(t_w) = sum over 1-runs of ceil(l/3) = P(A);
 *   (iii) nu(t_w + 1) = sum over 0-runs of ceil(l/3) = P(B);
 *   (iv) F = U_{0-runs [p0,q0]} { min(p0 + 3i + 2, q0 + 1) : i >= 0, p0 + 3i <= q0 } has |F| = P(B), each point is
 *        outside every 1-window or the first vertex of a 1-window, and (t_w + 1) - F has no 3-cycle (acyclic);
 *   (v) t_w is sigma-fixed iff w is a palindrome.
 * Mode "n3k": at n = 3k (and other n), for every TIGHT colouring w (any run lengths) test whether t_w is a tight line
 * with triple word w, and classify the failures and the words with an isolated 0 between two 1-runs.
 * usage: runwindow_audit theorem NMAX NEXH SAMPLES SEED | runwindow_audit n3k N
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
#define MAXN 33

static int n, L;
static uint32_t ot[MAXN], ou[MAXN];

static void build(const int *w) {
  /* windows of the 1-runs */
  static int rev[MAXN][MAXN];
  memset(rev, 0, sizeof rev);
  int a = 0;
  while (a < L) {
    if (!w[a]) { a++; continue; }
    int p = a;
    while (a < L && w[a]) a++;
    int q = a - 1;                              /* 1-run [p, q], window [p, q + 2] */
    for (int i = p; i <= q + 2; i++) for (int j = i + 2; j <= q + 2; j++) rev[i][j] = 1;
  }
  for (int v = 0; v < n; v++) ot[v] = ou[v] = 0;
  for (int v = 0; v + 1 < n; v++) { ot[v] |= 1u << (v + 1); ou[v] |= 1u << (v + 1); }
  for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) {
    if (rev[i][j]) { ot[j] |= 1u << i; ou[i] |= 1u << j; }
    else { ot[i] |= 1u << j; ou[j] |= 1u << i; }
  }
}

/* nu: triangles grouped by highest vertex, branch on the highest free vertex */
static uint32_t tri[MAXN][600]; static int ntri[MAXN];
static int bestp;
static void pack(uint32_t fr, int cur) {
  if (cur > bestp) bestp = cur;
  if (!fr) return;
  if (cur + __builtin_popcount(fr) / 3 <= bestp) return;
  int h = 31 - __builtin_clz(fr);
  for (int k = 0; k < ntri[h]; k++) if ((tri[h][k] & fr) == tri[h][k]) pack(fr & ~tri[h][k], cur + 1);
  pack(fr & ~(1u << h), cur);
}
static void triangles(const uint32_t *out) {
  uint32_t in[MAXN];
  for (int v = 0; v < n; v++) in[v] = 0;
  for (int v = 0; v < n; v++) { uint32_t o = out[v]; while (o) { int u = __builtin_ctz(o); o &= o - 1; in[u] |= 1u << v; } }
  for (int h = 0; h < n; h++) {
    ntri[h] = 0;
    uint32_t low = (h == 0) ? 0 : ((1u << h) - 1), A = out[h] & low;
    while (A) { int a = __builtin_ctz(A); A &= A - 1; uint32_t B = out[a] & in[h] & low;
      while (B) { int b = __builtin_ctz(B); B &= B - 1; if (ntri[h] < 600) tri[h][ntri[h]++] = (1u << h) | (1u << a) | (1u << b); else { fprintf(stderr, "tri overflow\n"); exit(3); } } }
  }
}
static int nu(const uint32_t *out) {
  triangles(out);
  bestp = 0;
  uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
  pack(full, 0);
  return bestp;
}
static int P_of(const int *w, int colour) {
  int c = 0, last = -100;
  for (int a = 0; a < L; a++) if (w[a] == colour && a - last >= 3) { c++; last = a; }
  return c;
}

static long long checked = 0, fails = 0;
static int check_theorem(const int *w) {
  build(w);
  int ok = 1;
  /* (i) */
  for (int a = 0; a < L; a++) { int bit = (ot[a + 2] >> a) & 1u; if (bit != w[a]) ok = 0; }
  /* run sums */
  int sumA = 0, sumB = 0;
  for (int a = 0; a < L; ) { int c = w[a], p = a; while (a < L && w[a] == c) a++; int l = a - p; if (c) sumA += (l + 2) / 3; else sumB += (l + 2) / 3; }
  if (P_of(w, 1) != sumA || P_of(w, 0) != sumB) ok = 0;
  /* (ii), (iii) */
  if (nu(ot) != sumA) ok = 0;
  int nuu = nu(ou);
  if (nuu != sumB) ok = 0;
  /* (iv) F; triangles of u are still in tri[] from the last nu call */
  uint32_t F = 0; int nF = 0;
  int inwin[MAXN] = {0}, firstwin[MAXN] = {0};
  for (int a = 0; a < L; ) {
    int c = w[a], p = a; while (a < L && w[a] == c) a++; int q = a - 1;
    if (c) { for (int v = p; v <= q + 2; v++) inwin[v] = 1; firstwin[p] = 1; }
    else for (int i = 0; p + 3 * i <= q; i++) { int x = p + 3 * i + 2; if (x > q + 1) x = q + 1; if (!((F >> x) & 1u)) { F |= 1u << x; nF++; } }
  }
  if (nF != sumB) ok = 0;
  for (int v = 0; v < n; v++) if (((F >> v) & 1u) && inwin[v] && !firstwin[v]) ok = 0;
  for (int h = 0; h < n; h++) for (int k = 0; k < ntri[h]; k++) if (!(tri[h][k] & F)) ok = 0;
  /* (v) colour */
  int blue = 1;
  for (int i = 0; i < n && blue; i++) for (int j = i + 2; j < n; j++) {
    int b1 = (ot[j] >> i) & 1u, b2 = (ot[n - 1 - i] >> (n - 1 - j)) & 1u;
    if (b1 != b2) { blue = 0; break; }
  }
  int pal = 1;
  for (int a = 0; a < L; a++) if (w[a] != w[L - 1 - a]) pal = 0;
  if (blue != pal) ok = 0;
  checked++;
  if (!ok) { fails++; if (fails <= 5) { printf("FAIL n=%d w=", n); for (int a = 0; a < L; a++) putchar('0' + w[a]); printf("\n"); } }
  return ok;
}

/* exhaustive enumeration of colourings with all runs >= 2 */
static int W[MAXN];
static void rec(int pos, int colour) {
  if (pos == L) { check_theorem(W); return; }
  for (int l = 2; pos + l <= L; l++) {
    for (int a = pos; a < pos + l; a++) W[a] = colour;
    if (pos + l == L || L - (pos + l) >= 2) rec(pos + l, 1 - colour);
  }
}
static u64 rs;
static u64 rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

int main(int argc, char **argv) {
  if (!strcmp(argv[1], "theorem")) {
    int nmax = atoi(argv[2]), nexh = atoi(argv[3]); long long S = atoll(argv[4]); rs = strtoull(argv[5], 0, 10) | 1;
    for (n = 4; n <= nmax; n++) {
      L = n - 2;
      long long c0 = checked, f0 = fails;
      if (n <= nexh) { rec(0, 0); rec(0, 1); }
      else {
        for (long long s = 0; s < S; s++) {
          /* random composition of L into parts >= 2, random first colour */
          int pos = 0, colour = (int)(rnd() & 1);
          while (pos < L) {
            int rem = L - pos, l;
            if (rem <= 3) l = rem;                  /* 2 or 3 left: one run */
            else { do { l = 2 + (int)(rnd() % (rem - 1)); } while (rem - l == 1); if (l > rem) l = rem; }
            for (int a = pos; a < pos + l; a++) W[a] = colour;
            pos += l; colour = 1 - colour;
          }
          check_theorem(W);
        }
      }
      printf("n=%d L=%d: %s %lld colourings with all runs >= 2, failures %lld\n", n, L, n <= nexh ? "all" : "sampled",
             checked - c0, fails - f0);
      fflush(stdout);
    }
    printf("total checked %lld, failures %lld\n", checked, fails);
    return 0;
  }
  /* mode n3k */
  n = atoi(argv[2]); L = n - 2; int K = n / 3;
  int ntight = 0, nreal = 0, niso = 0, niso_real = 0;
  for (u64 x = 0; x < (1ULL << L); x++) {
    if (!(x & 1ULL)) continue;
    int w[MAXN];
    for (int a = 0; a < L; a++) w[a] = (int)((x >> a) & 1ULL);
    if (P_of(w, 1) + P_of(w, 0) != K) continue;
    ntight++;
    build(w);
    int okw = 1;
    for (int a = 0; a < L; a++) if ((int)((ot[a + 2] >> a) & 1u) != w[a]) okw = 0;
    int real = okw && (nu(ot) + nu(ou) == K);
    nreal += real;
    int iso = 0;   /* an isolated 0 between two 1-runs: 1 0 1 at positions a-1, a, a+1 */
    for (int a = 1; a + 1 < L; a++) if (!w[a] && w[a - 1] && w[a + 1]) iso = 1;
    niso += iso; niso_real += iso && real;
    if (!real || iso) {
      printf("  w=");
      for (int a = 0; a < L; a++) putchar('0' + w[a]);
      printf(" isolated-0-between-1-runs=%d run-window realises=%d\n", iso, real);
    }
  }
  printf("n=%d: tight colourings (x_0 = 1) %d, realised by the run-window tiling %d; with an isolated 0 between two "
         "1-runs %d, of which realised %d\n", n, ntight, nreal, niso, niso_real);
  return 0;
}
