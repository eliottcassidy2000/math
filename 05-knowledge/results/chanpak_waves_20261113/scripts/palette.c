// Wave 26: passive-palette census. For every class R of order m (gentourng) and every mark set V
// (V = vertices beaten by the apex), print the block data of the exact cut transfer
//   a = H(R), q = P_R(2), (s,e) = total / odd-minus-even start weights from V,
//   (t,f) = total / odd-minus-even end weights into the complement of V, n = H(R + apex).
// Start weight of a nonempty path P = H(R - P). Mode "check" verifies the cut transfer on random
// two-block spines against direct Hamiltonian-path DP.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
#define MAXM 12

static i64 Hsub[1 << MAXM];
static i64 dpS[(1 << MAXM) * MAXM];

static i64 hamcount(int n, const int *out) {
  int N = 1 << n; i64 *dp = calloc((size_t)N * n, sizeof(i64));
  for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) { i64 x = dp[(size_t)S * n + v]; if (!x) continue;
    int nb = out[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | 1 << u) * n + u] += x; } }
  i64 s = 0; for (int v = 0; v < n; v++) s += dp[(size_t)(N - 1) * n + v]; free(dp); return s;
}

// per-vertex start and end weights
static void vertex_data(int m, const int *out, i64 *Sv, i64 *ev, i64 *Tv, i64 *fv, i64 *a, i64 *q) {
  int N = 1 << m;
  // Hsub[S] = H(R[S])
  memset(dpS, 0, sizeof(i64) * N * m);
  for (int v = 0; v < m; v++) dpS[(1 << v) * m + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < m; v++) { i64 x = dpS[S * m + v]; if (!x) continue;
    int nb = out[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dpS[(S | 1 << u) * m + u] += x; } }
  Hsub[0] = 1; for (int S = 1; S < N; S++) { i64 s = 0; for (int v = 0; v < m; v++) s += dpS[S * m + v]; Hsub[S] = s; }
  *a = Hsub[N - 1];
  i64 qq = 0; for (int S = 0; S < N; S++) qq += Hsub[S] * Hsub[(N - 1) ^ S]; *q = qq;
  for (int v = 0; v < m; v++) { Sv[v] = ev[v] = Tv[v] = fv[v] = 0; }
  // ends: dpS[S*m+v] = # Ham paths of R[S] ending at v
  for (int S = 1; S < N; S++) { int pc = __builtin_popcount(S); int sg = (pc & 1) ? 1 : -1;
    for (int v = 0; v < m; v++) { i64 x = dpS[S * m + v]; if (!x) continue; i64 w = x * Hsub[(N - 1) ^ S];
      Tv[v] += w; fv[v] += sg * w; } }
  // starts: paths starting at v covering S = Ham paths of reversed R[S] ending at v
  static i64 dpR[(1 << MAXM) * MAXM];
  memset(dpR, 0, sizeof(i64) * N * m);
  int in[MAXM]; for (int v = 0; v < m; v++) in[v] = 0;
  for (int v = 0; v < m; v++) for (int u = 0; u < m; u++) if ((out[v] >> u) & 1) in[u] |= 1 << v;
  for (int v = 0; v < m; v++) dpR[(1 << v) * m + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < m; v++) { i64 x = dpR[S * m + v]; if (!x) continue;
    int nb = in[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dpR[(S | 1 << u) * m + u] += x; } }
  for (int S = 1; S < N; S++) { int pc = __builtin_popcount(S); int sg = (pc & 1) ? 1 : -1;
    for (int v = 0; v < m; v++) { i64 x = dpR[S * m + v]; if (!x) continue; i64 w = x * Hsub[(N - 1) ^ S];
      Sv[v] += w; ev[v] += sg * w; } }
}

static void parse(const char *line, int n, int *out) {
  int k = 0; for (int i = 0; i < n; i++) out[i] = 0;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') out[i] |= 1 << j; else out[j] |= 1 << i; k++; }
}

// apex beats V: apex index m
static i64 apexcount(int m, const int *out, int V) {
  int o2[MAXM + 1]; for (int v = 0; v < m; v++) o2[v] = out[v] | (((V >> v) & 1) ? 0 : (1 << m)); o2[m] = V;
  return hamcount(m + 1, o2);
}

int main(int argc, char **argv) {
  if (!strcmp(argv[1], "check")) {
    srand(7);
    for (int trial = 0; trial < 300; trial++) {
      int m1 = 1 + rand() % 5, m2 = 1 + rand() % 5, n = m1 + m2;
      int out[MAXM + 1]; for (int i = 0; i <= n; i++) out[i] = 0;
      for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        int bi = i < m1, bj = j < m1;
        if (bi != bj) { out[i] |= 1 << j; continue; }           // block 1 => block 2
        if (rand() & 1) out[i] |= 1 << j; else out[j] |= 1 << i; }
      int V = rand() % (1 << n);
      // direct
      int o2[MAXM + 1]; for (int v = 0; v < n; v++) o2[v] = out[v] | (((V >> v) & 1) ? 0 : (1 << n)); o2[n] = V;
      i64 Nd = hamcount(n + 1, o2);
      // transfer: block 1 then block 2
      int out1[MAXM], out2[MAXM];
      for (int i = 0; i < m1; i++) out1[i] = out[i] & ((1 << m1) - 1);
      for (int i = 0; i < m2; i++) out2[i] = (out[m1 + i] >> m1) & ((1 << m2) - 1);
      int V1 = V & ((1 << m1) - 1), V2 = V >> m1;
      i64 S1[MAXM], e1[MAXM], T1[MAXM], f1[MAXM], a1, q1, S2[MAXM], e2[MAXM], T2[MAXM], f2[MAXM], a2, q2;
      vertex_data(m1, out1, S1, e1, T1, f1, &a1, &q1);
      i64 s1 = 0, ee1 = 0; for (int v = 0; v < m1; v++) if ((V1 >> v) & 1) { s1 += S1[v]; ee1 += e1[v]; }
      i64 n1 = apexcount(m1, out1, V1);
      vertex_data(m2, out2, S2, e2, T2, f2, &a2, &q2);
      i64 t2 = 0, ff2 = 0; for (int v = 0; v < m2; v++) if (!((V2 >> v) & 1)) { t2 += T2[v]; ff2 += f2[v]; }
      i64 n2 = apexcount(m2, out2, V2);
      // state after block 1: (h,N,S,d) = (a1, n1, s1, ee1)
      i64 Nt = a2 * n1 + (n2 - a2) * a1 + t2 * s1 + ff2 * ee1;
      if (Nt != Nd) { printf("MISMATCH trial %d: direct %lld transfer %lld\n", trial, Nd, Nt); return 1; }
    }
    printf("cut transfer verified on 300 random two-block spines\n");
    return 0;
  }
  int m = atoi(argv[2]); char line[4096];
  while (fgets(line, sizeof line, stdin)) {
    if (line[0] == '>') continue;
    int out[MAXM]; parse(line, m, out);
    i64 Sv[MAXM], ev[MAXM], Tv[MAXM], fv[MAXM], a, q;
    vertex_data(m, out, Sv, ev, Tv, fv, &a, &q);
    line[strcspn(line, "\n")] = 0;
    for (int V = 0; V < (1 << m); V++) {
      i64 s = 0, e = 0, t = 0, f = 0;
      for (int v = 0; v < m; v++) { if ((V >> v) & 1) { s += Sv[v]; e += ev[v]; } else { t += Tv[v]; f += fv[v]; } }
      i64 n = apexcount(m, out, V);
      printf("%lld %lld %lld %lld %lld %lld %lld %s %d\n", a, q, s, e, t, f, n, line, V);
    }
  }
  return 0;
}
