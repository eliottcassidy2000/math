/* AUDIT (independent): ingredients and censuses of THM-5071 (4b), (4c).
 * (i)   f(n) = min H over strong tournaments of order n, n = 3..7 (all labelled tournaments), and Moon's deletion
 *       property (every strong T of order 4..7 has v with T - v strong).
 * (ii)  Strong-resolution lemma: every strong semicomplete digraph of order 3..7 with exactly one reciprocal pair
 *       (placed at {0,1}; relabelling) has a strong tournament resolution; and, for ANY number of reciprocal pairs,
 *       every strong semicomplete digraph of order 3..5 has one (Camion / Session234 DG4 form).
 * (iii) One reciprocal pair, no missing pair (pair {0,1} doubled), orders 2..7: the set of H values, per order and
 *       cumulatively; strong ones at orders <= 4.
 * (iv)  One missing pair, no reciprocal pair (pair {0,1} empty), orders 2..7: the set of H values.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;
static int n, A[8][8];
static ll dp[1 << 8][8];

static ll ham(void) {
  int N = 1 << n;
  for (int S = 0; S < N; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
  for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) if (dp[S][v])
    for (int t = 0; t < n; t++) if (!(S >> t & 1) && A[v][t]) dp[S | 1 << t][t] += dp[S][v];
  ll h = 0;
  for (int v = 0; v < n; v++) h += dp[N - 1][v];
  return h;
}
static int strong_on(int mask) { /* induced subgraph on mask strongly connected? */
  if (!mask) return 1;
  int s = __builtin_ctz(mask), f = 1 << s, b = 1 << s, ch = 1;
  while (ch) { ch = 0;
    for (int v = 0; v < n; v++) if (mask >> v & 1) {
      if (f >> v & 1) for (int t = 0; t < n; t++) if ((mask >> t & 1) && A[v][t] && !(f >> t & 1)) { f |= 1 << t; ch = 1; }
      if (b >> v & 1) for (int t = 0; t < n; t++) if ((mask >> t & 1) && A[t][v] && !(b >> t & 1)) { b |= 1 << t; ch = 1; }
    } }
  return f == mask && b == mask;
}
static int strong(void) { return strong_on((1 << n) - 1); }
static void set_pairs(ll code, int skip01, int state01) {
  memset(A, 0, sizeof A); int k = 0;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
    if (skip01 && i == 0 && j == 1) { A[0][1] = state01 & 1; A[1][0] = state01 >> 1; continue; }
    if (code >> k & 1) A[i][j] = 1; else A[j][i] = 1;
    k++;
  }
}

int main(void) {
  /* (i) */
  for (n = 3; n <= 7; n++) {
    int m = n * (n - 1) / 2; ll fmin = -1, nstrong = 0, moonfail = 0;
    for (ll code = 0; code < (1LL << m); code++) {
      set_pairs(code, 0, 0);
      if (!strong()) continue;
      nstrong++;
      ll h = ham(); if (fmin < 0 || h < fmin) fmin = h;
      if (n >= 4) { int ok = 0; for (int v = 0; v < n && !ok; v++) if (strong_on(((1 << n) - 1) & ~(1 << v))) ok = 1; if (!ok) moonfail++; }
    }
    printf("(i) n=%d: strong labelled tournaments %lld, f(n) = %lld, strong ones without a strong one-vertex deletion: %lld\n",
           n, nstrong, fmin, moonfail);
    fflush(stdout);
  }
  /* (ii) general Camion form, orders 3..5, any number of reciprocal pairs */
  for (n = 3; n <= 5; n++) {
    int m = n * (n - 1) / 2; ll pw = 1; for (int i = 0; i < m; i++) pw *= 3;
    ll nstr = 0, bad = 0;
    for (ll code = 0; code < pw; code++) {
      memset(A, 0, sizeof A); ll c = code; int bi[10], bj[10], nb = 0;
      for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { int s2 = 1 + c % 3; c /= 3;
        A[i][j] = s2 & 1; A[j][i] = s2 >> 1; if (s2 == 3) { bi[nb] = i; bj[nb] = j; nb++; } }
      if (!strong()) continue;
      nstr++;
      int found = 0;
      for (int r = 0; r < (1 << nb) && !found; r++) {
        for (int t = 0; t < nb; t++) { if (r >> t & 1) A[bj[t]][bi[t]] = 0; else A[bi[t]][bj[t]] = 0; }
        if (strong()) found = 1;
        for (int t = 0; t < nb; t++) { A[bi[t]][bj[t]] = 1; A[bj[t]][bi[t]] = 1; }
      }
      if (!found) bad++;
    }
    printf("(ii) n=%d: strong semicomplete digraphs (any b) %lld, without a strong tournament resolution: %lld\n", n, nstr, bad);
  }
  /* (iii) one reciprocal pair; (ii) one-pair lemma up to order 7 */
  static char seenB[8][100000], seenZ[8][100000];
  for (n = 2; n <= 7; n++) {
    int m = n * (n - 1) / 2 - 1; ll nstr = 0, bad = 0;
    for (ll code = 0; code < (1LL << m); code++) {
      set_pairs(code, 1, 3);
      ll h = ham(); if (h < 100000) seenB[n][h] = 1;
      if (n >= 3 && strong()) {
        nstr++;
        A[1][0] = 0; int s1 = strong(); A[1][0] = 1; A[0][1] = 0; int s2 = strong(); A[0][1] = 1;
        if (!s1 && !s2) bad++;
        if (n <= 4 && (h == 4 || h == 7 || h == 1)) printf("   strong one-B-pair digraph with H=%lld at n=%d\n", h, n);
      }
      set_pairs(code, 1, 0);
      h = ham(); if (h < 100000) seenZ[n][h] = 1;
    }
    if (n >= 3) printf("(ii) n=%d: strong semicomplete digraphs with exactly one reciprocal pair (at {0,1}) %lld, with no strong resolution: %lld\n", n, nstr, bad);
  }
  for (n = 2; n <= 7; n++) {
    printf("(iii) one reciprocal pair, order %d: absent values in 1..40:", n);
    for (int h = 1; h <= 40; h++) if (!seenB[n][h]) printf(" %d", h);
    printf("\n");
  }
  char cumB[100000] = {0}, cumZ[100000] = {0}; int firstB[100000], firstZ[100000];
  for (int h = 0; h < 100000; h++) { firstB[h] = firstZ[h] = 0; }
  for (n = 2; n <= 7; n++) for (int h = 0; h < 100000; h++) {
    if (seenB[n][h] && !cumB[h]) { cumB[h] = 1; firstB[h] = n; }
    if (seenZ[n][h] && !cumZ[h]) { cumZ[h] = 1; firstZ[h] = n; }
  }
  printf("(iii) one reciprocal pair, orders <= 7: absent values in 1..250:");
  for (int h = 1; h <= 250; h++) if (!cumB[h]) printf(" %d", h);
  int lastB = 0; for (int h = 1; h < 100000; h++) if (cumB[h]) lastB = h;
  printf("\n      largest value %d; first order of 21: %d; values 4, 7 seen: %d %d\n", lastB, firstB[21], cumB[4], cumB[7]);
  printf("(iv) one missing pair, orders <= 7: absent values in 0..200:");
  for (int h = 0; h <= 200; h++) if (!cumZ[h]) printf(" %d", h);
  printf("\n      first order of 7: %d; first order of 21: %d\n", firstZ[7], firstZ[21]);
  return 0;
}
