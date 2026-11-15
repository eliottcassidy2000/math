// Union of two-tape carriers: spine [guard 1][x-tape][C][y-tape][guard 0] + one apex, over all strong marked cores C of
// orders 3..MAXO (gentourng classes, all masks). Counts via the exact cut transfer. Reports coverage of odd N <= X
// and, optionally, multiplicities in a window.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
#define MAXM 10
typedef struct { i64 a, q, s, e, t, f, n; } blk;
typedef struct { i64 h, N, S, d; } st;
static st stepb(st x, blk b) { st y; y.h = b.a * x.h; y.N = b.a * x.N + (b.n - b.a) * x.h + b.t * x.S + b.f * x.d; y.S = b.q * x.S + x.h * b.s; y.d = x.h * b.e; return y; }
static st step1(st x, int mu) { st y; y.h = x.h; y.N = x.N + (1 - mu) * (x.S + x.d); y.S = 2 * x.S + mu * x.h; y.d = mu * x.h; return y; }
static i64 Hsub[1 << MAXM], dpS[(1 << MAXM) * MAXM], dpR[(1 << MAXM) * MAXM];
static i64 hamcount(int n, const int *out) { int N = 1 << n; i64 *dp = calloc((size_t)N * n, sizeof(i64));
  for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) { i64 x = dp[(size_t)S * n + v]; if (!x) continue; int nb = out[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | 1 << u) * n + u] += x; } }
  i64 s = 0; for (int v = 0; v < n; v++) s += dp[(size_t)(N - 1) * n + v]; free(dp); return s; }
static void vdata(int m, const int *out, i64 *Sv, i64 *ev, i64 *Tv, i64 *fv, i64 *a, i64 *q) {
  int N = 1 << m; memset(dpS, 0, sizeof(i64) * N * m); for (int v = 0; v < m; v++) dpS[(1 << v) * m + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < m; v++) { i64 x = dpS[S * m + v]; if (!x) continue; int nb = out[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dpS[(S | 1 << u) * m + u] += x; } }
  Hsub[0] = 1; for (int S = 1; S < N; S++) { i64 s = 0; for (int v = 0; v < m; v++) s += dpS[S * m + v]; Hsub[S] = s; }
  *a = Hsub[N - 1]; i64 qq = 0; for (int S = 0; S < N; S++) qq += Hsub[S] * Hsub[(N - 1) ^ S]; *q = qq;
  for (int v = 0; v < m; v++) Sv[v] = ev[v] = Tv[v] = fv[v] = 0;
  for (int S = 1; S < N; S++) { int sg = (__builtin_popcount(S) & 1) ? 1 : -1; for (int v = 0; v < m; v++) { i64 x = dpS[S * m + v]; if (!x) continue; i64 w = x * Hsub[(N - 1) ^ S]; Tv[v] += w; fv[v] += sg * w; } }
  int in[MAXM]; for (int v = 0; v < m; v++) in[v] = 0; for (int v = 0; v < m; v++) for (int u = 0; u < m; u++) if ((out[v] >> u) & 1) in[u] |= 1 << v;
  memset(dpR, 0, sizeof(i64) * N * m); for (int v = 0; v < m; v++) dpR[(1 << v) * m + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < m; v++) { i64 x = dpR[S * m + v]; if (!x) continue; int nb = in[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dpR[(S | 1 << u) * m + u] += x; } }
  for (int S = 1; S < N; S++) { int sg = (__builtin_popcount(S) & 1) ? 1 : -1; for (int v = 0; v < m; v++) { i64 x = dpR[S * m + v]; if (!x) continue; i64 w = x * Hsub[(N - 1) ^ S]; Sv[v] += w; ev[v] += sg * w; } } }
static i64 apexcount(int m, const int *out, int V) { int o2[MAXM + 1]; for (int v = 0; v < m; v++) o2[v] = out[v] | (((V >> v) & 1) ? 0 : (1 << m)); o2[m] = V; return hamcount(m + 1, o2); }
static void parse(const char *line, int n, int *out) { int k = 0; for (int i = 0; i < n; i++) out[i] = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') out[i] |= 1 << j; else out[j] |= 1 << i; k++; } }
static int strong(int n, const int *out) { int reach = 1, fr = 1; while (fr) { int nf = 0; for (int v = 0; v < n; v++) if ((fr >> v) & 1) nf |= out[v]; nf &= ~reach; reach |= nf; fr = nf; }
  if (reach != (1 << n) - 1) return 0; int in[16]; for (int v = 0; v < n; v++) in[v] = 0; for (int v = 0; v < n; v++) for (int u = 0; u < n; u++) if ((out[v] >> u) & 1) in[u] |= 1 << v;
  reach = 1; fr = 1; while (fr) { int nf = 0; for (int v = 0; v < n; v++) if ((fr >> v) & 1) nf |= in[v]; nf &= ~reach; reach |= nf; fr = nf; } return reach == (1 << n) - 1; }
static unsigned char *hit; static int *mult; static i64 X, W0, W1; static long long cfg;
static void ydfs(st x) {  // x: state after C and some y-bits
  st fin = step1(x, 0);
  if (fin.N > X) return;
  cfg++; if (fin.N <= X) { i64 k = (fin.N - 1) / 2; hit[k >> 3] |= 1 << (k & 7); if (mult && fin.N >= W0 && fin.N <= W1) mult[(fin.N - W0) / 2]++; }
  ydfs(step1(x, 0)); ydfs(step1(x, 1));
}
static void xdfs(st x, blk C, int depth) {
  st c = stepb(x, C);
  if (c.S > X && depth > 0) return;
  ydfs(c);
  if (depth < 40) { xdfs(step1(x, 0), C, depth + 1); xdfs(step1(x, 1), C, depth + 1); }
}
int main(int argc, char **argv) {
  X = atoll(argv[1]); hit = calloc(X / 16 + 2, 1);
  int want_m[4096], want_c[4096], want_V[4096], nw = 0; FILE *fb = fopen("greedy_basis.txt", "r"); char ln[512];
  while (fgets(ln, sizeof ln, fb)) { if (sscanf(ln, "order %d class %d mask %d", &want_m[nw], &want_c[nw], &want_V[nw]) == 3) nw++; } fclose(fb);
  long long ncores = 0;
  for (int m = 3; m <= 5; m++) {
    char cmd[64]; sprintf(cmd, "gentourng -c %d 2>/dev/null", m); FILE *fp = popen(cmd, "r"); char line[512]; int cls = 0;
    while (fgets(line, sizeof line, fp)) { if (line[0] == '>') continue; cls++; int out[MAXM]; parse(line, m, out); if (!strong(m, out)) continue;
      i64 Sv[MAXM], ev[MAXM], Tv[MAXM], fv[MAXM], a, q; vdata(m, out, Sv, ev, Tv, fv, &a, &q);
      for (int V = 0; V < (1 << m); V++) { int ok = 0; for (int w = 0; w < nw; w++) if (want_m[w] == m && want_c[w] == cls && want_V[w] == V) ok = 1; if (!ok) continue;
        blk C; C.a = a; C.q = q; C.s = C.e = C.t = C.f = 0;
        for (int v = 0; v < m; v++) { if ((V >> v) & 1) { C.s += Sv[v]; C.e += ev[v]; } else { C.t += Tv[v]; C.f += fv[v]; } }
        C.n = apexcount(m, out, V); ncores++;
        st s0 = {1, 1, 0, 0}; s0 = step1(s0, 1); xdfs(s0, C, 0); } }
    pclose(fp); }
  i64 miss = 0, last = 0; for (i64 N = 333; N <= X; N += 2) { i64 k = (N - 1) / 2; if (!(hit[k >> 3] & (1 << (k & 7)))) { miss++; last = N; } }
  printf("greedy basis: %lld cores, X=%lld, odd misses in [333, X]: %lld, largest %lld\n", ncores, X, miss, last); return 0; }
