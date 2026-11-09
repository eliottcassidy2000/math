/* gen_classes.c - independent audit of THM-4950 (opus-2026-11-09-S87 packet), written from scratch by the auditor.
 *
 * Generates the isomorphism classes of tournaments of order 1..NMAX (default 9) by one-vertex extension, with an
 * individualisation-refinement canonical labelling of the auditor's own.  Any key of the form "min code over some
 * relabellings" separates non-isomorphic tournaments, so the generated lists are complete by induction whatever the
 * quality of the canonical form; the class counts (A000568: 1, 1, 2, 4, 12, 56, 456, 6880, 191536) then show that
 * there are no duplicates.
 *
 * For every class it computes: strong or not; the odd-cycle packing vector alpha_1..alpha_4 by a subset recursion on
 * the independence polynomial Z(U) = sum_k alpha_k(T[U]) x^k; H_3 = Z(V)(3); H = Z(V)(2) and, independently, the number
 * of Hamiltonian paths (OCF check).
 *
 * Code convention (shared with the packet so that witness codes can be decoded): bit b runs over the pairs i < j in
 * lexicographic order and is set iff i -> j.
 * Output: OUTDIR/acls_n.txt with lines "code strong a1 a2 a3 a4 H H3".
 * usage: gen_classes OUTDIR [NMAX]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
typedef long long i64;
#define MAXN 10

static int N;
static uint32_t OUTM[MAXN];

static void decode(u64 code, int n, uint32_t *out) {
  int b = 0;
  for (int i = 0; i < n; i++) out[i] = 0;
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++, b++) {
      if ((code >> b) & 1ULL) out[i] |= 1u << j; else out[j] |= 1u << i;
    }
}

static u64 encode_perm(const int *pos) { /* pos[v] = new label of vertex v */
  int who[MAXN];
  for (int v = 0; v < N; v++) who[pos[v]] = v;
  u64 c = 0; int b = 0;
  for (int p = 0; p < N; p++)
    for (int q = p + 1; q < N; q++, b++)
      if ((OUTM[who[p]] >> who[q]) & 1u) c |= 1ULL << b;
  return c;
}

/* ordered partition: col[v] in 0..k-1.  Refine by (old colour, out-degree into each colour class) until stable. */
static int lexcmp(const int *a, const int *b, int len) {
  for (int i = 0; i < len; i++) if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
  return 0;
}
static int refine(int *col) {
  int k = 0;
  for (int v = 0; v < N; v++) if (col[v] + 1 > k) k = col[v] + 1;
  for (;;) {
    int key[MAXN][MAXN + 1];
    for (int v = 0; v < N; v++) {
      key[v][0] = col[v];
      for (int c = 1; c <= k; c++) key[v][c] = 0;
      uint32_t o = OUTM[v];
      while (o) { int u = __builtin_ctz(o); o &= o - 1; key[v][col[u] + 1]++; }
    }
    int idx[MAXN];
    for (int v = 0; v < N; v++) idx[v] = v;
    for (int a = 1; a < N; a++) {           /* insertion sort by key */
      int x = idx[a], b = a - 1;
      while (b >= 0 && lexcmp(key[idx[b]], key[x], k + 1) > 0) { idx[b + 1] = idx[b]; b--; }
      idx[b + 1] = x;
    }
    int nk = 0, nc[MAXN];
    nc[idx[0]] = 0;
    for (int a = 1; a < N; a++) { if (lexcmp(key[idx[a]], key[idx[a - 1]], k + 1)) nk++; nc[idx[a]] = nk; }
    nk++;
    for (int v = 0; v < N; v++) col[v] = nc[v];
    if (nk == k) return k;
    k = nk;
  }
}

static u64 bestc; static int found;
static void canon_rec(int *col) {
  int k = refine(col);
  if (k == N) {
    u64 c = encode_perm(col);
    if (!found || c < bestc) { bestc = c; found = 1; }
    return;
  }
  int sz[MAXN] = {0};
  for (int v = 0; v < N; v++) sz[col[v]]++;
  int t = 0;
  while (sz[t] < 2) t++;                       /* first non-singleton colour class */
  for (int v = 0; v < N; v++) {
    if (col[v] != t) continue;
    int c2[MAXN];
    for (int u = 0; u < N; u++) c2[u] = (col[u] > t) ? col[u] + 1 : col[u];
    for (int u = 0; u < N; u++) if (col[u] == t && u != v) c2[u] = t + 1;
    canon_rec(c2);
  }
}
static u64 canonical(const uint32_t *out, int n) {
  N = n;
  for (int v = 0; v < n; v++) OUTM[v] = out[v];
  int col[MAXN];
  for (int v = 0; v < n; v++) col[v] = 0;
  found = 0;
  canon_rec(col);
  return bestc;
}

/* ---------- invariants ---------- */
static int is_strong(const uint32_t *out, int n) {
  if (n == 1) return 1;
  uint32_t full = (1u << n) - 1, in[MAXN];
  for (int v = 0; v < n; v++) in[v] = 0;
  for (int v = 0; v < n; v++) for (int u = 0; u < n; u++) if ((out[v] >> u) & 1u) in[u] |= 1u << v;
  uint32_t r = 1, last = 0;
  while (r != last) { last = r; for (int v = 0; v < n; v++) if ((r >> v) & 1u) r |= out[v]; }
  if (r != full) return 0;
  r = 1; last = 0;
  while (r != last) { last = r; for (int v = 0; v < n; v++) if ((r >> v) & 1u) r |= in[v]; }
  return r == full;
}

static i64 g[1 << MAXN][MAXN];
static i64 cyc[1 << MAXN];
static i64 Z[1 << MAXN][5];
/* cyc[S] = number of directed cycles with vertex set exactly S (|S| >= 3) */
static void cycles(const uint32_t *out, int n) {
  int F = 1 << n;
  for (int S = 0; S < F; S++) cyc[S] = 0;
  for (int m = 0; m < n; m++) {
    /* paths starting at m using vertices >= m only; index subsets of [m, n-1] directly in F */
    uint32_t upper = (uint32_t)(F - 1) & ~((1u << m) - 1);
    /* iterate subsets S of upper containing m in increasing numeric order */
    for (uint32_t S = upper; ; S = (S - 1) & upper) { for (int v = 0; v < n; v++) g[S][v] = 0; if (!S) break; }
    g[1u << m][m] = 1;
    for (uint32_t S = 1u << m; S <= upper; S++) {
      if ((S & upper) != S || !((S >> m) & 1u)) continue;
      for (int v = m; v < n; v++) {
        i64 c = g[S][v];
        if (!c) continue;
        uint32_t nb = out[v] & upper & ~S;
        while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; g[S | (1u << u)][u] += c; }
      }
      if (__builtin_popcount(S) >= 3) {
        i64 t = 0;
        for (int v = m + 1; v < n; v++) if (g[S][v] && ((out[v] >> m) & 1u)) t += g[S][v];
        cyc[S] = t;
      }
    }
  }
}
static void packing_poly(int n, i64 *alpha) {
  int F = 1 << n;
  for (int k = 0; k < 5; k++) Z[0][k] = 0;
  Z[0][0] = 1;
  for (int U = 1; U < F; U++) {
    int m = __builtin_ctz(U);
    int rest = U & ~(1 << m);
    for (int k = 0; k < 5; k++) Z[U][k] = Z[rest][k];
    for (int T = rest; ; T = (T - 1) & rest) {
      int S = T | (1 << m);
      int sz = __builtin_popcount(S);
      if (sz >= 3 && (sz & 1) && cyc[S]) {
        for (int k = 0; k < 4; k++) Z[U][k + 1] += cyc[S] * Z[U & ~S][k];
      }
      if (!T) break;
    }
  }
  for (int k = 0; k < 5; k++) alpha[k] = Z[F - 1][k];
}
static i64 ham_paths(const uint32_t *out, int n) {
  static i64 h[1 << MAXN][MAXN];
  int F = 1 << n;
  for (int S = 0; S < F; S++) for (int v = 0; v < n; v++) h[S][v] = 0;
  for (int v = 0; v < n; v++) h[1 << v][v] = 1;
  for (int S = 1; S < F; S++)
    for (int v = 0; v < n; v++) {
      if (!h[S][v]) continue;
      uint32_t nb = out[v] & ~(uint32_t)S;
      while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; h[S | (1 << u)][u] += h[S][v]; }
    }
  i64 t = 0;
  for (int v = 0; v < n; v++) t += h[F - 1][v];
  return t;
}

/* ---------- hash set of codes ---------- */
typedef struct { u64 *k; size_t cap, cnt; } hset;
static void hs_init(hset *h, size_t cap) { h->cap = cap; h->cnt = 0; h->k = malloc(cap * sizeof(u64)); memset(h->k, 0xff, cap * sizeof(u64)); }
static int hs_add(hset *h, u64 key) {
  size_t i = (size_t)((key * 0xD6E8FEB86659FD93ULL) >> 17) & (h->cap - 1);
  while (h->k[i] != ~0ULL) { if (h->k[i] == key) return 0; i = (i + 1) & (h->cap - 1); }
  h->k[i] = key; h->cnt++; return 1;
}
static int cmpu64(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }

int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: gen_classes OUTDIR [NMAX]\n"); return 1; }
  int nmax = argc > 2 ? atoi(argv[2]) : 9;
  u64 *cur = malloc(sizeof(u64)); size_t ncur = 1; cur[0] = 0;   /* order 1 */
  for (int n = 1; n <= nmax; n++) {
    if (n > 1) {
      hset hs; hs_init(&hs, 1u << 22);
      for (size_t r = 0; r < ncur; r++) {
        uint32_t base[MAXN];
        decode(cur[r], n - 1, base);
        for (uint32_t O = 0; O < (1u << (n - 1)); O++) {   /* O = out-set of the new vertex n-1 */
          uint32_t t[MAXN];
          for (int v = 0; v < n - 1; v++) t[v] = base[v];
          t[n - 1] = 0;
          for (int v = 0; v < n - 1; v++) { if ((O >> v) & 1u) t[n - 1] |= 1u << v; else t[v] |= 1u << (n - 1); }
          hs_add(&hs, canonical(t, n));
        }
      }
      free(cur);
      cur = malloc(hs.cnt * sizeof(u64)); ncur = 0;
      for (size_t i = 0; i < hs.cap; i++) if (hs.k[i] != ~0ULL) cur[ncur++] = hs.k[i];
      free(hs.k);
      qsort(cur, ncur, sizeof(u64), cmpu64);
    }
    char fn[1024]; snprintf(fn, sizeof fn, "%s/acls_%d.txt", argv[1], n);
    FILE *fo = fopen(fn, "w");
    long nstrong = 0, ocf_bad = 0;
    static unsigned char seen3[1 << 17];
    memset(seen3, 0, sizeof seen3);
    long distinct3 = 0; i64 maxh3 = 0;
    for (size_t r = 0; r < ncur; r++) {
      uint32_t t[MAXN]; decode(cur[r], n, t);
      int st = is_strong(t, n);
      cycles(t, n);
      i64 al[5]; packing_poly(n, al);
      i64 H2 = 0, H3 = 0, p2 = 1, p3 = 1;
      for (int k = 0; k < 5; k++) { H2 += al[k] * p2; H3 += al[k] * p3; p2 *= 2; p3 *= 3; }
      i64 hp = ham_paths(t, n);
      if (hp != H2) ocf_bad++;
      if (st) { nstrong++; if (!seen3[H3]) { seen3[H3] = 1; distinct3++; } if (H3 > maxh3) maxh3 = H3; }
      fprintf(fo, "%llu %d %lld %lld %lld %lld %lld %lld\n", cur[r], st, al[1], al[2], al[3], al[4], hp, H3);
    }
    fclose(fo);
    printf("n=%d classes %zu strong %ld distinct strong H_3 %ld max strong H_3 %lld OCF failures %ld\n",
           n, ncur, nstrong, distinct3, maxh3, ocf_bad);
    fflush(stdout);
  }
  return 0;
}
