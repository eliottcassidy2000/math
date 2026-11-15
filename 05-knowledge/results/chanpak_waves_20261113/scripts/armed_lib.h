// Wave 26: unit extensions by one new vertex with an arbitrary out-set.
// H(T + x_X) = H + 2 W(X), W(X) = sum_{u in X, v notin X} Pi(u,v),
// Pi(u,v) = sum over paths u -> v of even order >= 2 of H(T - path).
// Modes:
//   check  : verify the cut formula against direct DP on random tournaments
//   census : read gentourng lines (order n) and report, per class, H and whether 2H+1 / 2H-1 are
//            one-vertex extension counts (and how many out-sets give each)
//   tree D : read gentourng lines, and for each class with H in the allowed set, search an adaptive
//            full binary tree of unit extensions to depth D (out-set re-chosen at every node)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long i64;
#define MAXN 20
static i64 *he, *Hc, *g;
static int CAPN;
static long long nodes_visited = 0;
static long long node_budget = 0;

static void build(int n, const int *out, i64 *Pi, i64 *Hout) {
  int N = 1 << n;
  memset(he, 0, sizeof(i64) * (size_t)N * n);
  for (int v = 0; v < n; v++) he[(size_t)(1 << v) * n + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) {
    i64 x = he[(size_t)S * n + v]; if (!x) continue;
    int nb = out[v] & ~S;
    while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; he[(size_t)(S | 1 << u) * n + u] += x; }
  }
  Hc[0] = 1;
  for (int S = 1; S < N; S++) { i64 s = 0; for (int v = 0; v < n; v++) s += he[(size_t)S * n + v]; Hc[S] = s; }
  *Hout = Hc[N - 1];
  memset(Pi, 0, sizeof(i64) * n * n);
  for (int u = 0; u < n; u++) {
    memset(g, 0, sizeof(i64) * (size_t)N * n);
    g[(size_t)(1 << u) * n + u] = 1;
    for (int S = 1; S < N; S++) {
      if (!((S >> u) & 1)) continue;
      int pc = __builtin_popcount(S);
      for (int v = 0; v < n; v++) {
        i64 x = g[(size_t)S * n + v]; if (!x) continue;
        if (pc >= 2 && !(pc & 1)) Pi[u * n + v] += x * Hc[(N - 1) ^ S];
        int nb = out[v] & ~S;
        while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; g[(size_t)(S | 1 << w) * n + w] += x; }
      }
    }
  }
}

static i64 direct_H(int n, const int *out) {
  int N = 1 << n;
  i64 *dp = calloc((size_t)N * n, sizeof(i64));
  for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) {
    i64 x = dp[(size_t)S * n + v]; if (!x) continue;
    int nb = out[v] & ~S;
    while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | 1 << u) * n + u] += x; }
  }
  i64 s = 0; for (int v = 0; v < n; v++) s += dp[(size_t)(N - 1) * n + v];
  free(dp); return s;
}

// enumerate all out-sets X by Gray code; record X with 2W == H+1 (plus) and 2W == H-1 (minus)
static int enum_units(int n, const i64 *Pi, i64 H, int *plus, int *np, int *minus, int *nm, int cap) {
  i64 rowX[MAXN], colX[MAXN], rowAll[MAXN];
  for (int u = 0; u < n; u++) { rowX[u] = 0; colX[u] = 0; rowAll[u] = 0; for (int v = 0; v < n; v++) if (v != u) rowAll[u] += Pi[u * n + v]; }
  int X = 0; i64 W = 0; *np = 0; *nm = 0;
  long long total = 1LL << n;
  for (long long i = 0; i < total; i++) {
    if (i > 0) {
      int u = __builtin_ctzll((unsigned long long)i);
      if (!((X >> u) & 1)) {  // add u
        W += (rowAll[u] - rowX[u]) - colX[u];
        X |= 1 << u;
        for (int w = 0; w < n; w++) { rowX[w] += Pi[w * n + u]; colX[w] += Pi[u * n + w]; }
      } else {                // remove u
        X &= ~(1 << u);
        for (int w = 0; w < n; w++) { rowX[w] -= Pi[w * n + u]; colX[w] -= Pi[u * n + w]; }
        W -= (rowAll[u] - rowX[u]) - colX[u];
      }
    }
    if (2 * W == H + 1 && *np < cap) plus[(*np)++] = X;
    if (2 * W == H - 1 && *nm < cap) minus[(*nm)++] = X;
  }
  return 0;
}

static void add_vertex(int n, const int *out, int X, int *out2) {
  for (int v = 0; v < n; v++) out2[v] = out[v] | (((X >> v) & 1) ? 0 : (1 << n));
  out2[n] = X;
}

// stats per depth
static long long st_nodes[64], st_plus[64], st_minus[64], st_armed[64];

static int armed(int n, const int *out, int d, int depth) {
  if (d == 0) return 1;
  if (n + 1 > CAPN) return 0;
  if (node_budget && nodes_visited >= node_budget) return 0;
  nodes_visited++;
  i64 Pi[MAXN * MAXN], H;
  build(n, out, Pi, &H);
  int *plus = malloc(sizeof(int) * 4096), *minus = malloc(sizeof(int) * 4096), np, nm;
  enum_units(n, Pi, H, plus, &np, minus, &nm, 4096);
  st_nodes[depth]++; st_plus[depth] += np; st_minus[depth] += nm;
  int ok = 0;
  int out2[MAXN];
  for (int i = 0; i < np && !ok; i++) { add_vertex(n, out, plus[i], out2); if (armed(n + 1, out2, d - 1, depth + 1)) ok = 1; }
  if (ok) {
    ok = 0;
    for (int i = 0; i < nm && !ok; i++) { add_vertex(n, out, minus[i], out2); if (armed(n + 1, out2, d - 1, depth + 1)) ok = 1; }
  }
  if (ok) st_armed[depth]++;
  free(plus); free(minus);
  return ok;
}

static int parse(const char *line, int n, int *out) {
  int k = 0; for (int i = 0; i < n; i++) out[i] = 0;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
    if (line[k] == '1') out[i] |= 1 << j; else out[j] |= 1 << i; k++;
  }
  return 0;
}

