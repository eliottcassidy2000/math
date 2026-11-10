// d_single_realiser.c -- AUDIT check of the sentence "THM-4876 shows the witnesses must change with n" (section 6):
// for every class T read (order m <= 10, H(T) = h < 1024) count the odd n in [3h, 5h) that are extension values
// H(T + z_X) of THIS SINGLE base T, over all 2^m out-sets.  Output: for each h the best single-base coverage.
// usage: gentourng -q m [res/mod] | ./d_single_realiser m > file   (lines: h  window_size  best_covered  classes)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
#define MAXM 10
#define HM 1024
static int m, full, inm[MAXM], outm[MAXM];
static i64 E[1 << MAXM][MAXM], B[1 << MAXM][MAXM], Mx[MAXM][MAXM], a[MAXM], mu[MAXM][MAXM], val[1 << MAXM];
static unsigned char seen[5 * HM + 8];
int main(int argc, char **argv) {
  m = atoi(argv[1]); full = (1 << m) - 1; char line[256];
  static int best[HM], ncls[HM], distinctmax[HM]; memset(best, -1, sizeof best);
  while (fgets(line, sizeof line, stdin)) {
    if (line[0] != '0' && line[0] != '1') continue; int k = 0; for (int i = 0; i < m; i++) inm[i] = outm[i] = 0;
    for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if (line[k++] == '1') { outm[i] |= 1 << j; inm[j] |= 1 << i; } else { outm[j] |= 1 << i; inm[i] |= 1 << j; } }
    for (int S = 1; S <= full; S++) { int vv = S; while (vv) { int v = __builtin_ctz(vv); vv &= vv - 1; int R = S ^ (1 << v);
        if (!R) { E[S][v] = 1; continue; } i64 e = 0; int pp = R & inm[v]; while (pp) { int u = __builtin_ctz(pp); pp &= pp - 1; e += E[R][u]; } E[S][v] = e; } }
    i64 H = 0; for (int v = 0; v < m; v++) H += E[full][v]; if (H >= HM) continue;
    for (int S = 1; S <= full; S++) { int vv = S; while (vv) { int v = __builtin_ctz(vv); vv &= vv - 1; int R = S ^ (1 << v);
        if (!R) { B[S][v] = 1; continue; } i64 b = 0; int pp = R & outm[v]; while (pp) { int u = __builtin_ctz(pp); pp &= pp - 1; b += B[R][u]; } B[S][v] = b; } }
    for (int j = 0; j < m; j++) for (int i = 0; i < m; i++) { Mx[j][i] = 0; if (i == j) continue; int rest = full & ~(1 << i) & ~(1 << j); i64 s = 0;
      for (int S0 = rest;; S0 = (S0 - 1) & rest) { int S = S0 | (1 << j); s += E[S][j] * B[full ^ S][i]; if (!S0) break; } Mx[j][i] = s; }
    for (int i = 0; i < m; i++) { a[i] = B[full][i] - E[full][i]; for (int j = 0; j < m; j++) if (j != i) a[i] += Mx[j][i]; }
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) mu[i][j] = (i == j) ? 0 : Mx[i][j] + Mx[j][i];
    val[0] = H; int lo = 3 * (int)H, hi = 5 * (int)H; memset(seen, 0, (size_t)hi + 1); int cov = 0;
    for (int X = 1; X <= full; X++) { int t = 31 - __builtin_clz(X); int Y = X ^ (1 << t); i64 v = val[Y] + a[t]; int yy = Y; while (yy) { int i = __builtin_ctz(yy); yy &= yy - 1; v -= mu[i][t]; } val[X] = v;
      if (v >= lo && v < hi && !seen[v]) { seen[v] = 1; cov++; } }
    ncls[H]++; if (cov > best[H]) best[H] = cov;
  }
  for (int h = 1; h < HM; h += 2) if (ncls[h]) printf("%d %d %d %d\n", h, h, best[h], ncls[h]);   /* window [3h,5h) has exactly h odd numbers */
  return 0;
}
