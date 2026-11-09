/* Auditor: the single sigma-fixed tile (0, n-1) reversed gives a BLUE line with nu-pair (1, k) at every n >= 4.
 * (Every cycle of t uses the only backward arc n-1 -> 0, so nu(t) = 1; t + 1 has all consecutive 3-cycles, so
 * nu(t + 1) = k.)  At even n with k odd (n = 4 mod 6) this is a member of the "both-odd" family of HYP-9590.
 * Checked here with exact nu for n = 4..40.
 */
#include <stdio.h>
#include <stdint.h>
typedef uint64_t u64;
static int n;
static u64 outm[64];
static u64 tri_by[64][1200];
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
static void build(int comp) { /* t: only tile (0, n-1) reversed */
  for (int i = 0; i < n; i++) outm[i] = 0;
  for (int i = 0; i + 1 < n; i++) outm[i] |= 1ULL << (i + 1);
  for (int i = 0; i < n; i++)
    for (int j = i + 2; j < n; j++) {
      int r = (i == 0 && j == n - 1) ^ comp;
      if (r) outm[j] |= 1ULL << i; else outm[i] |= 1ULL << j;
    }
}
int main(void) {
  int bad = 0;
  for (n = 4; n <= 40; n++) {
    build(0);
    int a = nu_cur();
    build(1);
    int b = nu_cur();
    int ok = (a == 1 && b == n / 3);
    if (!ok) bad++;
    printf("n=%d k=%d: nu(t) = %d, nu(t+1) = %d %s%s\n", n, n / 3, a, b, ok ? "(1, k)" : "MISMATCH",
           (n % 2 == 0 && (n / 3) % 2 == 1) ? "  <- both-odd family member (even n, k odd)" : "");
  }
  printf("tile (0, n-1) is sigma-fixed: sigma(0, n-1) = (n-1-(n-1), n-1-0) = (0, n-1). mismatches: %d\n", bad);
  return 0;
}
