// a_labjoint.c -- AUDIT: exact joint distribution of (number of back arcs, H) over ALL labelled tournaments on n <= 7
// vertices. Written from scratch by the auditor. Arc {i<j} is a back arc when j -> i.  Output lines: "b H count".
// H is computed with a start-vertex DP: S[mask][v] = number of Hamiltonian paths of the sub-tournament on mask that start at v.
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64;
#define NMAX 7
int main(int argc, char **argv) {
  int n = atoi(argv[1]); if (n < 2 || n > NMAX) return 2;
  int np = n * (n - 1) / 2, full = (1 << n) - 1, HM = 4096;
  u64 *cnt = calloc((size_t)(np + 1) * HM, sizeof(u64));
  static u64 S[1 << NMAX][NMAX];
  for (u64 code = 0; code < (1ULL << np); code++) {
    int outm[NMAX] = {0}, k = 0, back = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (code >> k & 1) { outm[j] |= 1 << i; back++; } else outm[i] |= 1 << j; k++; }
    for (int mask = 1; mask <= full; mask++) for (int v = 0; v < n; v++) {
      if (!(mask >> v & 1)) continue; int rest = mask ^ (1 << v);
      if (!rest) { S[mask][v] = 1; continue; }
      u64 t = 0; int ww = rest & outm[v]; while (ww) { int w = __builtin_ctz(ww); ww &= ww - 1; t += S[rest][w]; } S[mask][v] = t; }
    u64 H = 0; for (int v = 0; v < n; v++) H += S[full][v];
    if (H >= (u64)HM) { fprintf(stderr, "H overflow\n"); return 1; }
    cnt[(size_t)back * HM + H]++;
  }
  for (int b = 0; b <= np; b++) for (int h = 0; h < HM; h++) if (cnt[(size_t)b * HM + h]) printf("%d %d %llu\n", b, h, cnt[(size_t)b * HM + h]);
  return 0;
}
