// Wave 26: search matched slides C z == w C (equal a, q, s, t) for every marked strong core of order m,
// tile widths d <= D. Slide step Delta = N(p C z q) - N(p w C q), p = last bit of w, q = first bit of z.
// Block data and spine counts use the exact cut transfer (h,N,S,d) -> (ah, aN+(n-a)h+tS+fd, qS+hs, he).
// Prints every solution with its step; a final check recomputes the two guarded spines by direct DP.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
#define MAXM 10
typedef struct { i64 a, q, s, e, t, f, n; } blk;

static i64 Hsub[1 << MAXM], dpS[(1 << MAXM) * MAXM], dpR[(1 << MAXM) * MAXM];
static i64 hamcount(int n, const int *out) {
  int N = 1 << n; i64 *dp = calloc((size_t)N * n, sizeof(i64));
  for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) { i64 x = dp[(size_t)S * n + v]; if (!x) continue;
    int nb = out[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | 1 << u) * n + u] += x; } }
  i64 s = 0; for (int v = 0; v < n; v++) s += dp[(size_t)(N - 1) * n + v]; free(dp); return s;
}
static void vertex_data(int m, const int *out, i64 *Sv, i64 *ev, i64 *Tv, i64 *fv, i64 *a, i64 *q) {
  int N = 1 << m;
  memset(dpS, 0, sizeof(i64) * N * m);
  for (int v = 0; v < m; v++) dpS[(1 << v) * m + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < m; v++) { i64 x = dpS[S * m + v]; if (!x) continue;
    int nb = out[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dpS[(S | 1 << u) * m + u] += x; } }
  Hsub[0] = 1; for (int S = 1; S < N; S++) { i64 s = 0; for (int v = 0; v < m; v++) s += dpS[S * m + v]; Hsub[S] = s; }
  *a = Hsub[N - 1]; i64 qq = 0; for (int S = 0; S < N; S++) qq += Hsub[S] * Hsub[(N - 1) ^ S]; *q = qq;
  for (int v = 0; v < m; v++) Sv[v] = ev[v] = Tv[v] = fv[v] = 0;
  for (int S = 1; S < N; S++) { int sg = (__builtin_popcount(S) & 1) ? 1 : -1;
    for (int v = 0; v < m; v++) { i64 x = dpS[S * m + v]; if (!x) continue; i64 w = x * Hsub[(N - 1) ^ S]; Tv[v] += w; fv[v] += sg * w; } }
  int in[MAXM]; for (int v = 0; v < m; v++) in[v] = 0;
  for (int v = 0; v < m; v++) for (int u = 0; u < m; u++) if ((out[v] >> u) & 1) in[u] |= 1 << v;
  memset(dpR, 0, sizeof(i64) * N * m);
  for (int v = 0; v < m; v++) dpR[(1 << v) * m + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < m; v++) { i64 x = dpR[S * m + v]; if (!x) continue;
    int nb = in[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dpR[(S | 1 << u) * m + u] += x; } }
  for (int S = 1; S < N; S++) { int sg = (__builtin_popcount(S) & 1) ? 1 : -1;
    for (int v = 0; v < m; v++) { i64 x = dpR[S * m + v]; if (!x) continue; i64 w = x * Hsub[(N - 1) ^ S]; Sv[v] += w; ev[v] += sg * w; } }
}
static i64 apexcount(int m, const int *out, int V) {
  int o2[MAXM + 1]; for (int v = 0; v < m; v++) o2[v] = out[v] | (((V >> v) & 1) ? 0 : (1 << m)); o2[m] = V;
  return hamcount(m + 1, o2);
}
static blk single(int mu) { blk b = {1, 2, mu, mu, 1 - mu, 1 - mu, 1}; return b; }
// state transfer
typedef struct { i64 h, N, S, d; } st;
static st start(void) { st x = {1, 1, 0, 0}; return x; }
static st step(st x, blk b) { st y; y.h = b.a * x.h; y.N = b.a * x.N + (b.n - b.a) * x.h + b.t * x.S + b.f * x.d; y.S = b.q * x.S + x.h * b.s; y.d = x.h * b.e; return y; }
static void parse(const char *line, int n, int *out) {
  int k = 0; for (int i = 0; i < n; i++) out[i] = 0;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') out[i] |= 1 << j; else out[j] |= 1 << i; k++; }
}
static int strong(int n, const int *out) {
  int reach = 1, frontier = 1;
  while (frontier) { int nf = 0; for (int v = 0; v < n; v++) if ((frontier >> v) & 1) nf |= out[v]; nf &= ~reach; reach |= nf; frontier = nf; }
  if (reach != (1 << n) - 1) return 0;
  int in[16]; for (int v = 0; v < n; v++) in[v] = 0;
  for (int v = 0; v < n; v++) for (int u = 0; u < n; u++) if ((out[v] >> u) & 1) in[u] |= 1 << v;
  reach = 1; frontier = 1;
  while (frontier) { int nf = 0; for (int v = 0; v < n; v++) if ((frontier >> v) & 1) nf |= in[v]; nf &= ~reach; reach |= nf; frontier = nf; }
  return reach == (1 << n) - 1;
}
// direct check of a guarded spine: p, (C or word), (word or C), q with apex; build tournament explicitly
static i64 direct_spine(int m, const int *outC, int VC, int pre_len, i64 pre_bits, int post_len, i64 post_bits, int g1, int g2) {
  // order: guard1, pre-word (pre_len bits, MSB first), C (m), post-word, guard2; ordinal between pieces
  int n = 1 + pre_len + m + post_len + 1; if (n + 1 > 24) return -1;
  int *out = calloc(n + 1, sizeof(int)); int marks = 0, idx = 0;
  int pos_C = 1 + pre_len;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
    int ci = (i >= pos_C && i < pos_C + m), cj = (j >= pos_C && j < pos_C + m);
    if (ci && cj) { if ((outC[i - pos_C] >> (j - pos_C)) & 1) out[i] |= 1 << j; else out[j] |= 1 << i; }
    else out[i] |= 1 << j;  // ordinal: earlier beats later
  }
  // marks
  marks |= g1 << 0; idx = 1;
  for (int b = 0; b < pre_len; b++) { if ((pre_bits >> (pre_len - 1 - b)) & 1) marks |= 1 << idx; idx++; }
  for (int v = 0; v < m; v++) { if ((VC >> v) & 1) marks |= 1 << idx; idx++; }
  for (int b = 0; b < post_len; b++) { if ((post_bits >> (post_len - 1 - b)) & 1) marks |= 1 << idx; idx++; }
  marks |= g2 << idx;
  for (int v = 0; v < n; v++) if (!((marks >> v) & 1)) out[v] |= 1 << n;
  out[n] = marks;
  i64 r = hamcount(n + 1, out); free(out); return r;
}
int main(int argc, char **argv) {
  int m = atoi(argv[1]), D = atoi(argv[2]);
  char line[4096]; int cls = 0; long long nsol = 0;
  while (fgets(line, sizeof line, stdin)) {
    if (line[0] == '>') continue; cls++;
    int out[MAXM]; parse(line, m, out);
    if (!strong(m, out)) continue;
    i64 Sv[MAXM], ev[MAXM], Tv[MAXM], fv[MAXM], a, q;
    vertex_data(m, out, Sv, ev, Tv, fv, &a, &q);
    for (int V = 0; V < (1 << m); V++) {
      blk C; C.a = a; C.q = q; C.s = C.e = C.t = C.f = 0;
      for (int v = 0; v < m; v++) { if ((V >> v) & 1) { C.s += Sv[v]; C.e += ev[v]; } else { C.t += Tv[v]; C.f += fv[v]; } }
      C.n = apexcount(m, out, V);
      for (int d = 1; d <= D; d++) {
        i64 M = (1LL << d) - 1;
        for (i64 z = 0; z <= M; z++) {
          i64 num = M * C.s + C.a * z; if (num % C.q) continue; i64 w = num / C.q; if (w > M) continue;
          // t-values of words: t_word = sum over unmarked positions i (1-based from first) of 2^(i-1)
          i64 tz = 0, tw = 0;
          for (int i = 1; i <= d; i++) { if (!((z >> (d - i)) & 1)) tz += 1LL << (i - 1); if (!((w >> (d - i)) & 1)) tw += 1LL << (i - 1); }
          if (C.q * tz != M * C.t + C.a * tw) continue;
          int p = (int)(w & 1), qg = (int)((z >> (d - 1)) & 1);
          // spine p C z q
          st x = start(); x = step(x, single(p)); x = step(x, C);
          for (int i = 1; i <= d; i++) x = step(x, single((z >> (d - i)) & 1));
          x = step(x, single(qg)); i64 N1 = x.N;
          st y = start(); y = step(y, single(p));
          for (int i = 1; i <= d; i++) y = step(y, single((w >> (d - i)) & 1));
          y = step(y, C); y = step(y, single(qg)); i64 N2 = y.N;
          i64 Delta = N1 - N2;
          if (Delta == 0) continue;
          nsol++;
          int ok = 1;
          if (1 + d + m + 1 + 1 <= 20) {
            i64 D1 = direct_spine(m, out, V, 0, 0, d, z, p, qg), D2 = direct_spine(m, out, V, d, w, 0, 0, p, qg);
            ok = (D1 == N1 && D2 == N2);
          }
          line[strcspn(line, "\n")] = 0;
          printf("m=%d class=%d code=%s mask=%d h=%lld Q=%lld d=%d z=%lld w=%lld guards=(%d,%d) N1=%lld N2=%lld Delta=%lld direct=%s\n",
                 m, cls, line, V, C.a, C.q, d, z, w, p, qg, N1, N2, Delta, ok ? "ok" : "MISMATCH");
        }
      }
    }
  }
  fprintf(stderr, "order %d: %lld nonzero slides\n", m, nsol);
  return 0;
}
