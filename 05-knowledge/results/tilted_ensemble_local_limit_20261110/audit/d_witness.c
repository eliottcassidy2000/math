// d_witness.c -- AUDIT check (d): explicit witnesses for deletion pairs (h, n) = (H(T), H(T + z_X)).
// usage: gentourng -q m | ./d_witness m targets.txt      (targets.txt: lines "h n")
// For each target not yet found prints:  h n m <gentourng line of T> <X as a list of vertices 0..m-1 beaten by z>
// The value H(T+z_X) is found with the cut quadratic and then RECOUNTED by a direct DP on the (m+1)-vertex tournament;
// only pairs confirmed by the direct count are printed.  Written from scratch by the auditor.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
#define MAXM 10
static int m, full, inm[MAXM], outm[MAXM];
static i64 E[1 << MAXM][MAXM], B[1 << MAXM][MAXM], Dv[1 << MAXM], Dw[1 << MAXM][MAXM], Mx[MAXM][MAXM], a[MAXM], mu[MAXM][MAXM], val[1 << MAXM];
static i64 direct_ext(int X) {
  Dv[0] = 1;
  for (int S = 1; S <= full; S++) {
    i64 t = 0; int uu = S & ~X; while (uu) { int u = __builtin_ctz(uu); uu &= uu - 1; t += E[S][u]; } Dv[S] = t;
    for (int w = 0; w < m; w++) { if (!(S >> w & 1)) { Dw[S][w] = 0; continue; } int R = S ^ (1 << w); i64 s = 0;
      if (X >> w & 1) s += Dv[R]; int pp = R & inm[w]; while (pp) { int u = __builtin_ctz(pp); pp &= pp - 1; s += Dw[R][u]; } Dw[S][w] = s; } }
  i64 tot = Dv[full]; for (int w = 0; w < m; w++) tot += Dw[full][w]; return tot;
}
int main(int argc, char **argv) {
  m = atoi(argv[1]); full = (1 << m) - 1; FILE *tf = fopen(argv[2], "r"); if (!tf) return 2;
  i64 th[256], tn[256]; int found[256], nt = 0; while (nt < 256 && fscanf(tf, "%lld %lld", &th[nt], &tn[nt]) == 2) { found[nt] = 0; nt++; } fclose(tf);
  char line[256]; int left = nt;
  while (left > 0 && fgets(line, sizeof line, stdin)) {
    if (line[0] != '0' && line[0] != '1') continue; line[strcspn(line, "\r\n")] = 0;
    int k = 0; for (int i = 0; i < m; i++) inm[i] = outm[i] = 0;
    for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if (line[k++] == '1') { outm[i] |= 1 << j; inm[j] |= 1 << i; } else { outm[j] |= 1 << i; inm[i] |= 1 << j; } }
    for (int S = 1; S <= full; S++) { int vv = S; while (vv) { int v = __builtin_ctz(vv); vv &= vv - 1; int R = S ^ (1 << v);
        if (!R) { E[S][v] = 1; continue; } i64 e = 0; int pp = R & inm[v]; while (pp) { int u = __builtin_ctz(pp); pp &= pp - 1; e += E[R][u]; } E[S][v] = e; } }
    i64 H = 0; for (int v = 0; v < m; v++) H += E[full][v];
    int want = 0; for (int t = 0; t < nt; t++) if (!found[t] && th[t] == H) want = 1; if (!want) continue;
    for (int S = 1; S <= full; S++) { int vv = S; while (vv) { int v = __builtin_ctz(vv); vv &= vv - 1; int R = S ^ (1 << v);
        if (!R) { B[S][v] = 1; continue; } i64 b = 0; int pp = R & outm[v]; while (pp) { int u = __builtin_ctz(pp); pp &= pp - 1; b += B[R][u]; } B[S][v] = b; } }
    for (int j = 0; j < m; j++) for (int i = 0; i < m; i++) { Mx[j][i] = 0; if (i == j) continue; i64 s = 0;
      for (int S = 1; S < full; S++) if ((S >> j & 1) && !(S >> i & 1)) s += E[S][j] * B[full ^ S][i]; Mx[j][i] = s; }
    for (int i = 0; i < m; i++) { a[i] = B[full][i] - E[full][i]; for (int j = 0; j < m; j++) if (j != i) a[i] += Mx[j][i]; }
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) mu[i][j] = (i == j) ? 0 : Mx[i][j] + Mx[j][i];
    val[0] = H;
    for (int X = 1; X <= full; X++) { int t = 31 - __builtin_clz(X); int Y = X ^ (1 << t); i64 v = val[Y] + a[t]; int yy = Y; while (yy) { int i = __builtin_ctz(yy); yy &= yy - 1; v -= mu[i][t]; } val[X] = v; }
    for (int t = 0; t < nt; t++) if (!found[t] && th[t] == H) for (int X = 0; X <= full; X++) if (val[X] == tn[t] && direct_ext(X) == tn[t]) {
        printf("%lld %lld %d %s", th[t], tn[t], m, line); for (int i = 0; i < m; i++) if (X >> i & 1) printf(" %d", i); printf("\n"); found[t] = 1; left--; break; }
  }
  for (int t = 0; t < nt; t++) if (!found[t]) fprintf(stderr, "NOT FOUND at order %d: %lld %lld\n", m, th[t], tn[t]);
  return 0;
}
