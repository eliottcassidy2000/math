// d_window.c -- AUDIT checks (c) and (d): H-histogram of tournament classes and the deletion-pair census
// (H(T), H(T+z_X)) over all out-sets X, for bases T read in gentourng format.  Written from scratch by the auditor.
//
// usage: gentourng -q m [res/mod] | ./d_window m prefix
// writes prefix.hist  ("H count" lines, all classes read),
//        prefix.A.bin (bitset [h/2][n/2], h < 1024 odd, n < 16384 odd: bit set iff some class T with H(T)=h and some
//                      out-set X has H(T+z_X) = n),
//        prefix.Z.bin (same, out-sets of size one only),
//        and prints a summary line with the number of classes, of bases with H < 1024, and of failed direct controls.
// Method: E[S][v] / B[S][v] = Hamiltonian paths of T[S] ending / starting at v;  M[j][i] = sum_{S: j in S, i not in S}
// E[S][j] B[V-S][i];  H(T+z_X) = H + sum_{i in X} a_i - sum_{pairs in X} mu_ij, a_i = g(i) - f(i) + sum_{j != i} M[j][i],
// mu_ij = M[i][j] + M[j][i].  CONTROL: for every base with H < 1024 one out-set (pseudo-random) is recounted by a direct
// DP on the (m+1)-vertex tournament that does not use the formula.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
#define MAXM 10
#define HM 1024
#define NM 16384
static int m, full;
static int inm[MAXM], outm[MAXM];
static i64 E[1 << MAXM][MAXM], B[1 << MAXM][MAXM], Dv[1 << MAXM], Dw[1 << MAXM][MAXM];
static i64 Mx[MAXM][MAXM], a[MAXM], mu[MAXM][MAXM], val[1 << MAXM];
static unsigned char *tA, *tZ;
#define ROW (NM / 16)
static void setbit(unsigned char *t, i64 h, i64 n) { size_t k = (size_t)(n >> 1); t[(size_t)(h >> 1) * ROW + (k >> 3)] |= (unsigned char)(1u << (k & 7)); }
static i64 direct_ext(int X) {
  Dv[0] = 1;
  for (int S = 1; S <= full; S++) {
    i64 t = 0; int uu = S & ~X;
    while (uu) { int u = __builtin_ctz(uu); uu &= uu - 1; t += E[S][u]; }
    Dv[S] = t;
    for (int w = 0; w < m; w++) {
      if (!(S >> w & 1)) { Dw[S][w] = 0; continue; }
      int R = S ^ (1 << w); i64 s = 0;
      if (X >> w & 1) s += Dv[R];
      int pp = R & inm[w];
      while (pp) { int u = __builtin_ctz(pp); pp &= pp - 1; s += Dw[R][u]; }
      Dw[S][w] = s;
    }
  }
  i64 tot = Dv[full]; for (int w = 0; w < m; w++) tot += Dw[full][w];
  return tot;
}
int main(int argc, char **argv) {
  if (argc < 3) return 2;
  m = atoi(argv[1]); full = (1 << m) - 1; const char *prefix = argv[2];
  int HH = 200000; u64 *hist = calloc(HH, sizeof(u64));
  tA = calloc((size_t)(HM / 2) * ROW, 1); tZ = calloc((size_t)(HM / 2) * ROW, 1);
  char line[256], fn[1024]; u64 nt = 0, nbase = 0, fail = 0, overflow = 0;
  while (fgets(line, sizeof line, stdin)) {
    if (line[0] != '0' && line[0] != '1') continue;
    int k = 0; for (int i = 0; i < m; i++) inm[i] = outm[i] = 0;
    for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if (line[k++] == '1') { outm[i] |= 1 << j; inm[j] |= 1 << i; } else { outm[j] |= 1 << i; inm[i] |= 1 << j; } }
    nt++;
    for (int S = 1; S <= full; S++) { int vv = S; while (vv) { int v = __builtin_ctz(vv); vv &= vv - 1; int R = S ^ (1 << v);
        if (!R) { E[S][v] = 1; continue; } i64 e = 0; int pp = R & inm[v]; while (pp) { int u = __builtin_ctz(pp); pp &= pp - 1; e += E[R][u]; } E[S][v] = e; } }
    i64 H = 0; for (int v = 0; v < m; v++) H += E[full][v];
    if (H >= HH) { fprintf(stderr, "H too large\n"); return 1; }
    hist[H]++;
    if (H >= HM) continue;
    nbase++;
    for (int S = 1; S <= full; S++) { int vv = S; while (vv) { int v = __builtin_ctz(vv); vv &= vv - 1; int R = S ^ (1 << v);
        if (!R) { B[S][v] = 1; continue; } i64 b = 0; int pp = R & outm[v]; while (pp) { int u = __builtin_ctz(pp); pp &= pp - 1; b += B[R][u]; } B[S][v] = b; } }
    for (int j = 0; j < m; j++) for (int i = 0; i < m; i++) { Mx[j][i] = 0; if (i == j) continue;
      int rest = full & ~(1 << i) & ~(1 << j); i64 s = 0;
      for (int S0 = rest;; S0 = (S0 - 1) & rest) { int S = S0 | (1 << j); if (S != full) s += E[S][j] * B[full ^ S][i]; if (!S0) break; }
      Mx[j][i] = s; }
    for (int i = 0; i < m; i++) { a[i] = B[full][i] - E[full][i]; for (int j = 0; j < m; j++) if (j != i) a[i] += Mx[j][i]; }
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) mu[i][j] = (i == j) ? 0 : Mx[i][j] + Mx[j][i];
    val[0] = H;
    for (int X = 1; X <= full; X++) { int t = 31 - __builtin_clz(X); int Y = X ^ (1 << t); i64 v = val[Y] + a[t];
      int yy = Y; while (yy) { int i = __builtin_ctz(yy); yy &= yy - 1; v -= mu[i][t]; } val[X] = v; }
    int Xc = (int)(((nt + 12345) * 2654435761ULL >> 9) & (u64)full);
    if (direct_ext(Xc) != val[Xc]) fail++;
    for (int X = 0; X <= full; X++) { i64 n = val[X]; if (n < 1 || !(n & 1)) { fail++; continue; } if (n >= NM) { overflow++; continue; }
      setbit(tA, H, n); if (__builtin_popcount(X) == 1) setbit(tZ, H, n); }
  }
  FILE *fp;
  snprintf(fn, sizeof fn, "%s.hist", prefix); fp = fopen(fn, "w"); for (int h = 0; h < HH; h++) if (hist[h]) fprintf(fp, "%d %llu\n", h, hist[h]); fclose(fp);
  snprintf(fn, sizeof fn, "%s.A.bin", prefix); fp = fopen(fn, "wb"); fwrite(tA, 1, (size_t)(HM / 2) * ROW, fp); fclose(fp);
  snprintf(fn, sizeof fn, "%s.Z.bin", prefix); fp = fopen(fn, "wb"); fwrite(tZ, 1, (size_t)(HM / 2) * ROW, fp); fclose(fp);
  printf("order %d: classes read %llu, bases with H < %d: %llu, failed controls %llu, extension values >= %d (not recorded): %llu\n", m, nt, HM, nbase, fail, NM, overflow);
  return 0;
}
