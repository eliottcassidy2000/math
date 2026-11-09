/* HYP-9590 (B) probe 2 (auditor): a structured blue family.  For a palindrome w of length L = n - 2 and a set M
 * of mirror tiles (a, n-1-a) (sigma-fixed), let t = R(w) XOR M (flip the tiles of M).  t is sigma-fixed, so the
 * line {t, t+1} is blue.  Records which residual pairs (blue (k,k) at n = 0,1 mod 3; blue non-tight both-odd
 * pairs at even n) occur, with one witness each (re-verified).  |M| <= MMAX.
 * usage: aud_mirror NMIN NMAX MMAX
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef uint64_t u64;
#define MAXN 40
static int n, L;
static unsigned char rv[MAXN][MAXN];
static u64 outm[MAXN];
static u64 tri_by[MAXN][800];
static int ntri_by[MAXN], best;
static void rec(u64 avail, int cur) {
  if (cur > best) best = cur;
  if (!avail) return;
  if (cur + __builtin_popcountll(avail) / 3 <= best) return;
  int v = __builtin_ctzll(avail);
  for (int t = 0; t < ntri_by[v]; t++) {
    u64 T = tri_by[v][t];
    if ((T & avail) == T) rec(avail & ~T, cur + 1);
  }
  rec(avail & ~(1ULL << v), cur);
}
static int nu_cur(void) {
  for (int v = 0; v < n; v++) ntri_by[v] = 0;
  for (int a = 0; a < n; a++)
    for (int b = a + 1; b < n; b++)
      for (int c = b + 1; c < n; c++) {
        int ab = (outm[a] >> b) & 1, bc = (outm[b] >> c) & 1, ca = (outm[c] >> a) & 1;
        if ((ab && bc && ca) || (!ab && !bc && !ca)) tri_by[a][ntri_by[a]++] = (1ULL << a) | (1ULL << b) | (1ULL << c);
      }
  best = 0;
  rec((1ULL << n) - 1, 0);
  return best;
}
static void build(int comp) {
  for (int i = 0; i < n; i++) outm[i] = 0;
  for (int i = 0; i + 1 < n; i++) outm[i] |= 1ULL << (i + 1);
  for (int i = 0; i < n; i++)
    for (int j = i + 2; j < n; j++) {
      int r = rv[i][j] ^ comp;
      if (r) outm[j] |= 1ULL << i; else outm[i] |= 1ULL << j;
    }
}
static int w[MAXN];
static void set_runwindow(void) {
  memset(rv, 0, sizeof rv);
  int a = 0;
  while (a < L) {
    if (w[a] != 1) { a++; continue; }
    int p = a;
    while (a < L && w[a] == 1) a++;
    int q = a - 1;
    for (int i = p; i <= q + 2; i++)
      for (int j = i + 2; j <= q + 2; j++) rv[i][j] = 1;
  }
}
static int sigma_fixed(void) {
  for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) if (rv[i][j] != rv[n - 1 - j][n - 1 - i]) return 0;
  return 1;
}
static int found[MAXN][MAXN];
int main(int argc, char **argv) {
  int nmin = atoi(argv[1]), nmax = atoi(argv[2]), mmax = atoi(argv[3]);
  for (n = nmin; n <= nmax; n++) {
    L = n - 2;
    int k = n / 3;
    memset(found, 0, sizeof found);
    int mir[MAXN], nm = 0;
    for (int a = 0; n - 1 - a >= a + 2; a++) mir[nm++] = a;
    for (u64 half = 0; half < (1ULL << ((L + 1) / 2)); half++) {
      for (int i = 0; i < (L + 1) / 2; i++) { w[i] = (half >> i) & 1; w[L - 1 - i] = w[i]; }
      for (u64 M = 0; M < (1ULL << nm); M++) {
        if (__builtin_popcountll(M) > mmax) continue;
        set_runwindow();
        for (int e = 0; e < nm; e++) if ((M >> e) & 1) rv[mir[e]][n - 1 - mir[e]] ^= 1;
        build(0);
        int a = nu_cur();
        build(1);
        int b = nu_cur();
        int lo = a < b ? a : b, hi = a < b ? b : a;
        int fam = 0;
        if (lo == k && hi == k && (n % 3 != 2)) fam = 1;
        if (n % 2 == 0 && lo % 2 == 1 && hi % 2 == 1 && lo + hi > k) fam = 1;
        if (fam && !found[lo][hi]) {
          if (!sigma_fixed()) { printf("BUG: not sigma-fixed\n"); return 1; }
          found[lo][hi] = 1;
          printf("n=%d k=%d blue (%d, %d): w = ", n, k, lo, hi);
          for (int i = 0; i < L; i++) putchar('0' + w[i]);
          printf(", flipped mirror tiles:");
          for (int e = 0; e < nm; e++) if ((M >> e) & 1) printf(" (%d,%d)", mir[e], n - 1 - mir[e]);
          printf("; nu(t) = %d, nu(t+1) = %d\n", a, b);
          fflush(stdout);
        }
      }
    }
    /* list family members not found */
    printf("n=%d k=%d: family members NOT found in this structured family:", n, k);
    for (int i = 1; i <= k; i++)
      for (int j = i; j <= k; j++) {
        int fam = 0;
        if (i == k && j == k && (n % 3 != 2)) fam = 1;
        if (n % 2 == 0 && i % 2 == 1 && j % 2 == 1 && i + j > k) fam = 1;
        if (fam && !found[i][j]) printf(" (%d, %d)", i, j);
      }
    printf("\n");
    fflush(stdout);
  }
  return 0;
}
