// Base-4 trees: two new vertices per step, H -> 4H + eps, eps in {-3,-1,1,3}; branches whose target is 7 or 21 are skipped.
#include "armed_lib.h"
typedef struct { int X, Y; } pr;
static int collect(int n, const int *out, i64 H, pr lists[4][4096], int cnt[4]) {
  i64 Pi[MAXN*MAXN], H0; build(n, out, Pi, &H0);
  for (int k = 0; k < 4; k++) cnt[k] = 0;
  int out1[MAXN];
  for (int X = 0; X < (1 << n); X++) {
    i64 W = 0; for (int u = 0; u < n; u++) if ((X >> u) & 1) for (int v = 0; v < n; v++) if (!((X >> v) & 1)) W += Pi[u*n+v];
    i64 H1 = H + 2*W;
    // quick bound: H2 >= H1 and H2 <= ... skip if H1 > 4H+3
    if (H1 > 4*H + 3) continue;
    add_vertex(n, out, X, out1);
    i64 Pi1[MAXN*MAXN], Hc1; build(n+1, out1, Pi1, &Hc1);
    int plus[1], minus[1]; (void)plus; (void)minus;
    // enumerate Y by Gray code
    int m = n + 1; i64 rowX[MAXN], colX[MAXN], rowAll[MAXN];
    for (int u = 0; u < m; u++) { rowX[u]=0; colX[u]=0; rowAll[u]=0; for (int v = 0; v < m; v++) if (v!=u) rowAll[u] += Pi1[u*m+v]; }
    int Y = 0; i64 W1 = 0;
    for (long long i = 0; i < (1LL << m); i++) {
      if (i > 0) { int u = __builtin_ctzll((unsigned long long)i);
        if (!((Y >> u) & 1)) { W1 += (rowAll[u]-rowX[u]) - colX[u]; Y |= 1<<u; for (int w=0; w<m; w++){ rowX[w]+=Pi1[w*m+u]; colX[w]+=Pi1[u*m+w]; } }
        else { Y &= ~(1<<u); for (int w=0; w<m; w++){ rowX[w]-=Pi1[w*m+u]; colX[w]-=Pi1[u*m+w]; } W1 -= (rowAll[u]-rowX[u]) - colX[u]; } }
      i64 H2 = H1 + 2*W1; i64 eps = H2 - 4*H;
      int k = eps == -3 ? 0 : eps == -1 ? 1 : eps == 1 ? 2 : eps == 3 ? 3 : -1;
      if (k >= 0 && cnt[k] < 4096) { lists[k][cnt[k]].X = X; lists[k][cnt[k]].Y = Y; cnt[k]++; }
    }
  }
  return 0;
}
static long long nodes;
static int armed4(int n, const int *out, int d) {
  if (d == 0) return 1;
  if (n + 2 > CAPN) return 0;
  nodes++;
  i64 Pi[MAXN*MAXN], H; build(n, out, Pi, &H);
  static pr L[24][4][4096]; int cnt[4];
  pr (*lists)[4096] = L[d];
  collect(n, out, H, lists, cnt);
  int eps[4] = {-3,-1,1,3};
  for (int k = 0; k < 4; k++) {
    i64 tgt = 4*H + eps[k]; if (tgt == 7 || tgt == 21 || tgt < 1) continue;
    int ok = 0; int out1[MAXN], out2[MAXN];
    int c = cnt[k];
    for (int i = 0; i < c && !ok; i++) { add_vertex(n, out, lists[k][i].X, out1); add_vertex(n+1, out1, lists[k][i].Y, out2);
      if (armed4(n+2, out2, d-1)) ok = 1; }
    if (!ok) return 0;
  }
  return 1;
}
int main(int argc, char **argv) {
  int n = atoi(argv[1]), D = atoi(argv[2]); CAPN = n + 2*D + 1; if (CAPN > 15) CAPN = 15;
  size_t NN = (size_t)1 << CAPN; he = malloc(sizeof(i64)*NN*CAPN); g = malloc(sizeof(i64)*NN*CAPN); Hc = malloc(sizeof(i64)*NN);
  char line[4096]; int idx = 0;
  while (fgets(line, sizeof line, stdin)) {
    if (line[0] == '>') continue; idx++;
    int out[MAXN]; parse(line, n, out); i64 H = direct_H(n, out); if (H == 1) continue;
    int best = 0;
    for (int d = 1; d <= D; d++) { nodes = 0; if (armed4(n, out, d)) best = d; else break; }
    line[strcspn(line, "\n")] = 0; printf("class %d H=%lld base4depth=%d code=%s\n", idx, H, best, line); fflush(stdout);
  }
  return 0;
}
