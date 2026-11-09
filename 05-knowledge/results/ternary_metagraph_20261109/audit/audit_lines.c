/* audit_lines.c - independent audit of THM-4951 (4), (5a), (5c) and of the census rows (written from scratch).
 *
 * Frame (as stated in THM-4951): base path 0 -> 1 -> ... -> n-1; tiles (i, j) with j >= i + 2; a tiling t reverses
 * the marked tiles (j -> i); the line through t is {t, t + 1} with t + 1 = complement of t; sigma(i, j) =
 * (n-1-j, n-1-i); a line is blue iff t = sigma(t).
 * nu(T) = maximum number of vertex-disjoint 3-cycles (computed by the auditor's own branch on the HIGHEST free vertex).
 *
 * Mode "all":  every line at order n (each line once, via its end with tile (0, 2) reversed), counted by unordered
 *              nu-pair and colour.  For blue lines also the pair (H mod 8) and (H mod 4) by an exact Hamiltonian-path DP.
 * Mode "blue": only the sigma-fixed tilings (enumerated orbit by orbit), with nu-pairs and H mod 4 / mod 8 pairs.
 * usage: audit_lines n all|blue
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
typedef unsigned long long u64;
#define MAXN 12

static int n, m;
static int TI[80], TJ[80], IDX[MAXN][MAXN], SIG[80];

static void build(u64 t, uint32_t *out) {
  for (int v = 0; v < n; v++) out[v] = 0;
  for (int v = 0; v + 1 < n; v++) out[v] |= 1u << (v + 1);
  for (int b = 0; b < m; b++) {
    if ((t >> b) & 1ULL) out[TJ[b]] |= 1u << TI[b]; else out[TI[b]] |= 1u << TJ[b];
  }
}

/* maximum number of disjoint 3-cycles: triangles are grouped by their highest vertex; branch on the highest free
 * vertex h: either h is left uncovered, or it is covered by a triangle whose other two vertices are free and lower. */
typedef struct { uint32_t tri[MAXN][64]; int nt[MAXN]; } trilist;
static void triangles(const uint32_t *out, trilist *L) {
  uint32_t in[MAXN];
  for (int v = 0; v < n; v++) in[v] = 0;
  for (int v = 0; v < n; v++) { uint32_t o = out[v]; while (o) { int u = __builtin_ctz(o); o &= o - 1; in[u] |= 1u << v; } }
  for (int h = 0; h < n; h++) {
    L->nt[h] = 0;
    uint32_t low = (1u << h) - 1;
    /* triangles h -> a -> b -> h with a, b < h */
    uint32_t A = out[h] & low;
    while (A) {
      int a = __builtin_ctz(A); A &= A - 1;
      uint32_t B = out[a] & in[h] & low;
      while (B) { int b = __builtin_ctz(B); B &= B - 1; L->tri[h][L->nt[h]++] = (1u << h) | (1u << a) | (1u << b); }
    }
  }
}
static int best;
#pragma omp threadprivate(best)
static void pack(const trilist *L, uint32_t freev, int cur) {
  if (cur > best) best = cur;
  if (!freev) return;
  if (cur + __builtin_popcount(freev) / 3 <= best) return;
  int h = 31 - __builtin_clz(freev);
  for (int k = 0; k < L->nt[h]; k++) if ((L->tri[h][k] & freev) == L->tri[h][k]) pack(L, freev & ~L->tri[h][k], cur + 1);
  pack(L, freev & ~(1u << h), cur);
}
static int nu_of(const uint32_t *out) {
  trilist L; triangles(out, &L);
  best = 0; pack(&L, (1u << n) - 1, 0);
  return best;
}

static long long hpaths(const uint32_t *out) {
  static long long h[1 << MAXN][MAXN];
#pragma omp threadprivate(h)
  int F = 1 << n;
  for (int S = 0; S < F; S++) for (int v = 0; v < n; v++) h[S][v] = 0;
  for (int v = 0; v < n; v++) h[1 << v][v] = 1;
  for (int S = 1; S < F; S++)
    for (int v = 0; v < n; v++) {
      long long c = h[S][v];
      if (!c) continue;
      uint32_t nb = out[v] & ~(uint32_t)S;
      while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; h[S | (1 << u)][u] += c; }
    }
  long long t = 0;
  for (int v = 0; v < n; v++) t += h[F - 1][v];
  return t;
}

int main(int argc, char **argv) {
  n = atoi(argv[1]);
  int blue_only = argc > 2 && (!strcmp(argv[2], "blue") || !strcmp(argv[2], "bluenu"));
  int skipH = argc > 2 && !strcmp(argv[2], "bluenu");     /* blue lines, nu-pairs only (used at n = 12) */
  m = 0;
  for (int j = 2; j < n; j++) for (int i = 0; i + 2 <= j; i++) { TI[m] = i; TJ[m] = j; IDX[i][j] = m; m++; }
  for (int b = 0; b < m; b++) SIG[b] = IDX[n - 1 - TJ[b]][n - 1 - TI[b]];
  u64 all = (m == 64) ? ~0ULL : ((1ULL << m) - 1);
  int b02 = IDX[0][2];
  static long long cnt[5][5][2];
  static long long h4[2][2], h8[8][8];
  long long nlines = 0, nblue = 0;
  if (!blue_only) {
    /* enumerate tilings with tile (0,2) reversed: the other m - 1 bits are free */
    u64 half = 1ULL << (m - 1);
#pragma omp parallel
    {
      long long lc[5][5][2]; memset(lc, 0, sizeof lc);
      long long l4[2][2], l8[8][8]; memset(l4, 0, sizeof l4); memset(l8, 0, sizeof l8);
      long long ln = 0, lb = 0;
#pragma omp for schedule(dynamic, 4096)
      for (long long r = 0; r < (long long)half; r++) {
        /* insert the forced bit b02 = 1 into r */
        u64 lowmask = (1ULL << b02) - 1;
        u64 t = ((u64)r & lowmask) | (1ULL << b02) | (((u64)r & ~lowmask) << 1);
        u64 u = t ^ all;
        uint32_t ot[MAXN], ou[MAXN];
        build(t, ot); build(u, ou);
        int a = nu_of(ot), c = nu_of(ou);
        u64 st = 0;
        for (int b = 0; b < m; b++) if ((t >> b) & 1ULL) st |= 1ULL << SIG[b];
        int blue = (st == t);
        int lo = a < c ? a : c, hi = a < c ? c : a;
        lc[lo][hi][blue]++;
        ln++;
        if (blue) {
          lb++;
          long long H1 = hpaths(ot), H2 = hpaths(ou);
          int x = (int)(H1 & 7), y = (int)(H2 & 7);
          if (x > y) { int z = x; x = y; y = z; }
          l8[x][y]++;
          int p = (int)((H1 >> 1) & 1), q = (int)((H2 >> 1) & 1);
          if (p > q) { int z = p; p = q; q = z; }
          l4[p][q]++;
        }
      }
#pragma omp critical
      {
        for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) for (int k = 0; k < 2; k++) cnt[i][j][k] += lc[i][j][k];
        for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) h4[i][j] += l4[i][j];
        for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) h8[i][j] += l8[i][j];
        nlines += ln; nblue += lb;
      }
    }
  } else {
    /* sigma-orbits of tiles */
    int orb[80], no = 0, seenb[80] = {0};
    u64 omask[80];
    for (int b = 0; b < m; b++) if (!seenb[b]) { omask[no] = (1ULL << b) | (1ULL << SIG[b]); seenb[b] = seenb[SIG[b]] = 1; orb[no++] = b; }
    long long half = 1LL << (no - 1);
    int o02 = -1;
    for (int o = 0; o < no; o++) if ((omask[o] >> b02) & 1ULL) o02 = o;
#pragma omp parallel
    {
      long long lc[5][5][2]; memset(lc, 0, sizeof lc);
      long long l4[2][2], l8[8][8]; memset(l4, 0, sizeof l4); memset(l8, 0, sizeof l8);
      long long ln = 0;
#pragma omp for schedule(dynamic, 1024)
      for (long long r = 0; r < half; r++) {
        u64 t = omask[o02];
        int k = 0;
        for (int o = 0; o < no; o++) { if (o == o02) continue; if ((r >> k) & 1) t |= omask[o]; k++; }
        u64 u = t ^ all;
        uint32_t ot[MAXN], ou[MAXN];
        build(t, ot); build(u, ou);
        int a = nu_of(ot), c = nu_of(ou);
        int lo = a < c ? a : c, hi = a < c ? c : a;
        lc[lo][hi][1]++;
        ln++;
        if (skipH) continue;
        long long H1 = hpaths(ot), H2 = hpaths(ou);
        int x = (int)(H1 & 7), y = (int)(H2 & 7);
        if (x > y) { int z = x; x = y; y = z; }
        l8[x][y]++;
        int p = (int)((H1 >> 1) & 1), q = (int)((H2 >> 1) & 1);
        if (p > q) { int z = p; p = q; q = z; }
        l4[p][q]++;
      }
#pragma omp critical
      {
        for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) for (int kk = 0; kk < 2; kk++) cnt[i][j][kk] += lc[i][j][kk];
        for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) h4[i][j] += l4[i][j];
        for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) h8[i][j] += l8[i][j];
        nlines += ln; nblue += ln;
      }
    }
  }
  printf("n=%d tiles m=%d mode=%s: lines %lld (expected %llu), blue lines %lld\n", n, m, blue_only ? "blue" : "all",
         nlines, blue_only ? 0ULL : (1ULL << (m - 1)), nblue);
  for (int i = 0; i < 5; i++) for (int j = i; j < 5; j++) for (int k = 0; k < 2; k++)
    if (cnt[i][j][k]) printf("  nu-pair (%d,%d) %s %lld\n", i, j, k ? "blue" : "black", cnt[i][j][k]);
  printf("  blue lines by H mod 4 pair: (1,1) %lld (1,3) %lld (3,3) %lld\n", h4[0][0], h4[0][1], h4[1][1]);
  printf("  blue lines by H mod 8 pair:");
  for (int i = 0; i < 8; i++) for (int j = i; j < 8; j++) if (h8[i][j]) printf(" (%d,%d):%lld", i, j, h8[i][j]);
  printf("\n");
  return 0;
}
