/* HYP-9590 (B) probe (auditor): search for lines {t, t+1} of the tiling metagraph whose unordered nu-pair is one of
 * the three residual families (black (k,k) at n = 3k; blue (k,k) at n = 0,1 mod 3; blue non-tight both-odd pairs at
 * even n), by random restarts plus a greedy walk on tiles (black) or sigma-orbits of tiles (blue).
 * Every witness is re-verified (nu of both ends recomputed from scratch, colour recomputed) and printed as its list
 * of reversed tiles, so it can be checked independently.
 * usage: aud_residual NMIN NMAX SEED MAXSTEPS
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef uint64_t u64;
#define MAXN 32
static int n, m;
static int ti[600], tj[600], orb[600], norb, orbrep[600];
static unsigned char bit[600]; /* reversed in t */
static u64 outm[MAXN];
static u64 rng = 88172645463325252ULL;
static u64 rnd(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

static void build(int comp) {
  for (int i = 0; i < n; i++) outm[i] = 0;
  for (int i = 0; i + 1 < n; i++) outm[i] |= 1ULL << (i + 1);
  for (int e = 0; e < m; e++) {
    int r = bit[e] ^ comp;
    if (r) outm[tj[e]] |= 1ULL << ti[e]; else outm[ti[e]] |= 1ULL << tj[e];
  }
}
static u64 tri_by[MAXN][600];
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
static int tile_index(int i, int j) {
  for (int e = 0; e < m; e++) if (ti[e] == i && tj[e] == j) return e;
  return -1;
}
static int sigma_fixed(void) {
  for (int e = 0; e < m; e++) if (bit[e] != bit[tile_index(n - 1 - tj[e], n - 1 - ti[e])]) return 0;
  return 1;
}
static void pair_now(int *a, int *b) {
  build(0); *a = nu_cur();
  build(1); *b = nu_cur();
}
static int dist(int a, int b, int i, int j) {
  int d1 = abs(a - i) + abs(b - j), d2 = abs(a - j) + abs(b - i);
  return d1 < d2 ? d1 : d2;
}
static void flip(int unit, int blue) {
  if (!blue) { bit[unit] ^= 1; return; }
  for (int e = 0; e < m; e++) if (orb[e] == unit) bit[e] ^= 1;
}

static int search(int i, int j, int blue, long maxsteps) {
  int units = blue ? norb : m;
  long steps = 0;
  while (steps < maxsteps) {
    /* restart: random density */
    double dens = (rnd() % 1000) / 1000.0;
    if (blue) {
      for (int o = 0; o < norb; o++) {
        int val = (rnd() % 1000) < dens * 1000;
        for (int e = 0; e < m; e++) if (orb[e] == o) bit[e] = val;
      }
    } else
      for (int e = 0; e < m; e++) bit[e] = (rnd() % 1000) < dens * 1000;
    int a, b;
    pair_now(&a, &b);
    int d = dist(a, b, i, j);
    for (int it = 0; it < 400 && d > 0; it++, steps++) {
      int u = rnd() % units;
      flip(u, blue);
      int a2, b2;
      pair_now(&a2, &b2);
      int d2 = dist(a2, b2, i, j);
      if (d2 <= d) { d = d2; a = a2; b = b2; } else flip(u, blue);
    }
    if (d == 0) {
      /* re-verify from scratch */
      int a3, b3;
      pair_now(&a3, &b3);
      int lo = a3 < b3 ? a3 : b3, hi = a3 < b3 ? b3 : a3;
      int col = sigma_fixed();
      if (lo == i && hi == j && col == blue) {
        printf("n=%d k=%d %s (%d, %d): FOUND after %ld steps; nu(t) = %d, nu(t+1) = %d; reversed tiles of t:", n, n / 3,
               blue ? "blue" : "black", i, j, steps, a3, b3);
        for (int e = 0; e < m; e++) if (bit[e]) printf(" %d-%d", ti[e], tj[e]);
        printf("\n");
        fflush(stdout);
        return 1;
      }
    }
  }
  printf("n=%d k=%d %s (%d, %d): not found in %ld steps\n", n, n / 3, blue ? "blue" : "black", i, j, maxsteps);
  fflush(stdout);
  return 0;
}

int main(int argc, char **argv) {
  int nmin = atoi(argv[1]), nmax = atoi(argv[2]);
  rng ^= (u64)atoll(argv[3]) * 0x9E3779B97F4A7C15ULL;
  long maxsteps = atol(argv[4]);
  for (n = nmin; n <= nmax; n++) {
    m = 0;
    for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { ti[m] = i; tj[m] = j; m++; }
    norb = 0;
    for (int e = 0; e < m; e++) orb[e] = -1;
    for (int e = 0; e < m; e++) {
      if (orb[e] >= 0) continue;
      int f = tile_index(n - 1 - tj[e], n - 1 - ti[e]);
      orb[e] = norb; orb[f] = norb; orbrep[norb] = e; norb++;
    }
    int k = n / 3;
    if (n % 3 == 0) search(k, k, 0, maxsteps);
    if (n % 3 == 0 || n % 3 == 1) search(k, k, 1, maxsteps);
    if (n % 2 == 0)
      for (int i = 1; i <= k; i += 2)
        for (int j = i; j <= k; j += 2)
          if (i + j > k && !(i == k && j == k && (n % 3 != 2))) search(i, j, 1, maxsteps);
  }
  return 0;
}
