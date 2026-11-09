/* Auditor: explicit lines for the residual families of HYP-9590 (B).
 * Nested mirror tiles: t_m reverses exactly the sigma-fixed tiles (a, n-1-a), 0 <= a < m (0 <= m <= k).
 *   Claim: t_m is blue with nu(t_m) = m and nu(t_m + 1) = k.
 *   - t_m: its backward arcs are n-1-a -> a (a < m), pairwise vertex-disjoint; every cycle uses one, so
 *     nu <= m; the triangles {a, x_a, n-1-a} with distinct x_a in [m, n-1-m] (n - 2m >= m) are disjoint 3-cycles.
 *   - t_m + 1: no flipped tile has the form (x, x+2) when m <= k < (n-1)/2 (n >= 6), so every consecutive window
 *     is a 3-cycle and nu = P([0, n-3]) = k.
 *   m = k gives blue (k, k) at every n >= 6; at even n with k odd, odd m gives the both-odd pairs (m, k).
 * Black (k, k) at n = 3k: t_k with the tile (0, 3) also reversed.  The packings {a, x_a, n-1-a} (x_0 != 3) and
 *   the consecutive windows {3b, 3b+1, 3b+2} avoid the arc 0-3, so both survive; (0, 3) != sigma(0, 3).
 * Checked with exact nu.
 */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64;
static int n;
static unsigned char rv[64][64];
static u64 outm[64];
static u64 tri_by[64][2000];
static int ntri_by[64], best;
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
static int sigma_fixed(void) {
  for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) if (rv[i][j] != rv[n - 1 - j][n - 1 - i]) return 0;
  return 1;
}
int main(int argc, char **argv) {
  int nmax = 40, bad = 0, checks = 0;
  for (n = 6; n <= nmax; n++) {
    int k = n / 3;
    for (int m = 0; m <= k; m++) {
      memset(rv, 0, sizeof rv);
      for (int a = 0; a < m; a++) rv[a][n - 1 - a] = 1;
      build(0);
      int x = nu_cur();
      build(1);
      int y = nu_cur();
      int blue = sigma_fixed();
      checks++;
      if (!(x == m && y == k && blue)) { bad++; printf("FAIL nested n=%d m=%d: nu = (%d, %d), blue %d\n", n, m, x, y, blue); }
    }
    printf("n=%d k=%d: nested mirror lines t_m, m = 0..k: nu-pairs (m, k), all blue%s\n", n, k,
           (n % 2 == 0 && k % 2 == 1) ? "  [both-odd members (m, k), m odd, realised]" : "");
    if (n % 3 == 0) {
      memset(rv, 0, sizeof rv);
      for (int a = 0; a < k; a++) rv[a][n - 1 - a] = 1;
      rv[0][3] ^= 1;
      build(0);
      int x = nu_cur();
      build(1);
      int y = nu_cur();
      int blue = sigma_fixed();
      checks++;
      if (!(x == k && y == k && !blue)) { bad++; printf("FAIL black n=%d: nu = (%d, %d), blue %d\n", n, x, y, blue); }
      else printf("n=%d k=%d: t_k + tile (0,3): black (k, k) = (%d, %d)\n", n, k, x, y);
    }
    fflush(stdout);
  }
  printf("checks %d, failures %d\n", checks, bad);
  return 0;
}
