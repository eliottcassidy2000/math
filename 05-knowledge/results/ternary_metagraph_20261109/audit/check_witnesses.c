/* check_witnesses.c - auditor's re-check of witnesses_order10.txt (lines 'H_3 code') with the auditor's own code:
 * decode (bit b over pairs i < j in lexicographic order, set iff i -> j), order 10, strong, H_3 = I(Omega, 3) by the
 * subset recursion, and OCF: I(Omega, 2) equals the Hamiltonian-path count.
 * usage: check_witnesses witnesses_order10.txt
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
typedef long long i64;
#define NN 10

static void decode(u64 code, int n, uint32_t *out) {
  int b = 0;
  for (int i = 0; i < n; i++) out[i] = 0;
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++, b++) {
      if ((code >> b) & 1ULL) out[i] |= 1u << j; else out[j] |= 1u << i;
    }
}
static int strong(const uint32_t *out, int n) {
  uint32_t full = (1u << n) - 1, in[NN] = {0};
  for (int v = 0; v < n; v++) for (int u = 0; u < n; u++) if ((out[v] >> u) & 1u) in[u] |= 1u << v;
  uint32_t r = 1, l = 0;
  while (r != l) { l = r; for (int v = 0; v < n; v++) if ((r >> v) & 1u) r |= out[v]; }
  if (r != full) return 0;
  r = 1; l = 0;
  while (r != l) { l = r; for (int v = 0; v < n; v++) if ((r >> v) & 1u) r |= in[v]; }
  return r == full;
}
static i64 g[1 << NN][NN], cyc[1 << NN], z3[1 << NN], z2[1 << NN], hh[1 << NN][NN];
static void values(const uint32_t *out, int n, i64 *H3, i64 *H2, i64 *HP) {
  int F = 1 << n;
  for (int S = 0; S < F; S++) cyc[S] = 0;
  for (int m = 0; m < n; m++) {
    for (int S = 0; S < F; S++) for (int v = 0; v < n; v++) g[S][v] = 0;
    g[1 << m][m] = 1;
    for (int S = 1 << m; S < F; S++) {
      if (!((S >> m) & 1) || (S & ((1 << m) - 1))) continue;
      for (int v = m; v < n; v++) {
        i64 c = g[S][v];
        if (!c) continue;
        uint32_t nb = out[v] & ~(uint32_t)S;
        while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; g[S | (1 << u)][u] += c; }
      }
      int sz = __builtin_popcount(S);
      if (sz >= 3 && (sz & 1)) for (int v = m + 1; v < n; v++) if ((out[v] >> m) & 1u) cyc[S] += g[S][v];
    }
  }
  z3[0] = z2[0] = 1;
  for (int U = 1; U < F; U++) {
    int m = __builtin_ctz(U), rest = U & ~(1 << m);
    z3[U] = z3[rest]; z2[U] = z2[rest];
    for (int T = rest; ; T = (T - 1) & rest) {
      int S = T | (1 << m);
      if (cyc[S]) { z3[U] += 3 * cyc[S] * z3[U & ~S]; z2[U] += 2 * cyc[S] * z2[U & ~S]; }
      if (!T) break;
    }
  }
  *H3 = z3[F - 1]; *H2 = z2[F - 1];
  for (int S = 0; S < F; S++) for (int v = 0; v < n; v++) hh[S][v] = 0;
  for (int v = 0; v < n; v++) hh[1 << v][v] = 1;
  for (int S = 1; S < F; S++) for (int v = 0; v < n; v++) {
    i64 c = hh[S][v];
    if (!c) continue;
    uint32_t nb = out[v] & ~(uint32_t)S;
    while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; hh[S | (1 << u)][u] += c; }
  }
  i64 t = 0;
  for (int v = 0; v < n; v++) t += hh[F - 1][v];
  *HP = t;
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r");
  long long v; u64 code; long tot = 0, bad = 0, prev = -1, unsorted = 0;
  while (fscanf(f, "%lld %llu", &v, &code) == 2) {
    uint32_t out[NN]; decode(code, NN, out);
    i64 H3, H2, HP;
    values(out, NN, &H3, &H2, &HP);
    int ok = strong(out, NN) && H3 == v && H2 == HP && (code >> 45) == 0;
    if (!ok) { bad++; if (bad < 10) printf("FAIL value %lld code %llu: H3 %lld H2 %lld HP %lld\n", v, code, H3, H2, HP); }
    if (v <= prev) unsorted++;
    prev = v; tot++;
  }
  printf("witnesses checked %ld, failures %ld, non-increasing value lines %ld\n", tot, bad, unsorted);
  return 0;
}
