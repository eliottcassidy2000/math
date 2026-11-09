/* intercyclic_gen.c - opus-2026-11-08-S86.
 * Exhaustive generation of the STRONG INTERCYCLIC tournaments (no two vertex-disjoint cycles; equivalently, by Moon,
 * no two vertex-disjoint cyclic triangles) of orders 3..NMAX, up to isomorphism.
 *
 * Completeness: intercyclicity is hereditary, and every strong tournament of order n >= 4 has a vertex whose deletion
 * leaves a strong tournament (Moon). Hence every strong intercyclic class of order n+1 is a one-vertex extension of a
 * strong intercyclic class of order n, and level-by-level extension with isomorphism rejection is exhaustive.
 * Canonical forms: individualisation-refinement with ordered colour refinement, adapted from kind-pasteur's classgen.c
 * (05-knowledge/results/first_creators_genealogy_20261107/classgen.c), widened to 128-bit codes and 32-bit masks.
 *
 * For every class it reports: n, alpha_1 (number of directed odd cycles; H = 1 + 2 alpha_1 since alpha_2 = 0),
 * H (Hamiltonian paths, computed independently by DP as a check), the cycle-transversal number tau (1 or 2; tau <= 2 is
 * Bang-Jensen-Bessy Theorem 1.8), and self-converse.
 * usage: intercyclic_gen NMAX [HCHECK] > out.txt  (line: n alpha1 H tau SC code_hex; H computed for n <= HCHECK)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned __int128 u128;
typedef unsigned long long u64;
#define MAXN 16      /* 128-bit canonical codes hold n(n-1)/2 <= 120 bits: orders <= 16 only */
static int nn;
static uint32_t A[MAXN];
static u128 best; static int have;
static u128 code_of(const int *pos){
  int inv[MAXN]; for (int v = 0; v < nn; v++) inv[pos[v]] = v;
  u128 c = 0; int b = 0;
  for (int i = 0; i < nn; i++) for (int j = i + 1; j < nn; j++){ if (A[inv[i]] >> inv[j] & 1) c |= ((u128)1) << b; b++; }
  return c;
}
static int ncells(const int *cell){ int m = -1; for (int v = 0; v < nn; v++) if (cell[v] > m) m = cell[v]; return m + 1; }
static void refine(int *cell){
  for (;;){
    int k = ncells(cell);
    static int sig[MAXN][MAXN + 1];
    for (int v = 0; v < nn; v++){ sig[v][0] = cell[v]; for (int c = 0; c < k; c++) sig[v][c + 1] = 0;
      for (int u = 0; u < nn; u++) if (A[v] >> u & 1) sig[v][cell[u] + 1]++; }
    int ord[MAXN]; for (int v = 0; v < nn; v++) ord[v] = v;
    for (int i = 1; i < nn; i++){ int x = ord[i], j = i - 1; while (j >= 0 && memcmp(sig[ord[j]], sig[x], (k + 1) * sizeof(int)) > 0){ ord[j + 1] = ord[j]; j--; } ord[j + 1] = x; }
    int nc[MAXN]; int c = 0; nc[ord[0]] = 0;
    for (int i = 1; i < nn; i++){ if (memcmp(sig[ord[i]], sig[ord[i - 1]], (k + 1) * sizeof(int))) c++; nc[ord[i]] = c; }
    int changed = (c + 1 != k);
    memcpy(cell, nc, sizeof(int) * nn);
    if (!changed) return;
  }
}
static void search(int *cell){
  refine(cell);
  int k = ncells(cell);
  if (k == nn){ u128 c = code_of(cell); if (!have || c < best){ best = c; have = 1; } return; }
  int cnt[MAXN] = {0}; for (int v = 0; v < nn; v++) cnt[cell[v]]++;
  int target = -1; for (int c = 0; c < k; c++) if (cnt[c] > 1){ target = c; break; }
  for (int v = 0; v < nn; v++) if (cell[v] == target){
    int c2[MAXN]; for (int u = 0; u < nn; u++){ c2[u] = cell[u] > target ? cell[u] + 1 : cell[u]; if (cell[u] == target && u != v) c2[u] = target + 1; }
    search(c2);
  }
}
static u128 canon(const uint32_t *adj, int n){
  nn = n; memcpy(A, adj, sizeof(uint32_t) * n); have = 0;
  int cell[MAXN]; for (int v = 0; v < n; v++) cell[v] = 0;
  search(cell); return best;
}
static void decode(u128 c, int n, uint32_t *adj){ memset(adj, 0, sizeof(uint32_t) * n); int b = 0;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++){ if (c >> b & 1) adj[i] |= 1u << j; else adj[j] |= 1u << i; b++; } }
static int strongq(const uint32_t *adj, int n){
  uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1), r = 1, prev = 0;
  while (r != prev){ prev = r; for (int v = 0; v < n; v++) if (r >> v & 1) r |= adj[v]; }
  if (r != full) return 0;
  r = 1; prev = 0;
  while (r != prev){ prev = r; for (int v = 0; v < n; v++) if (!(r >> v & 1)){ if (adj[v] & r) r |= 1u << v; } }
  return r == full;
}
/* list cyclic triangles as masks */
static int triangles(const uint32_t *adj, int n, uint32_t *tri){
  int t = 0;
  for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) for (int c = b + 1; c < n; c++){
    int ab = adj[a] >> b & 1, bc = adj[b] >> c & 1, ca = adj[c] >> a & 1;
    if (ab == bc && bc == ca) tri[t++] = (1u << a) | (1u << b) | (1u << c);
  }
  return t;
}
static int intercyclic(const uint32_t *adj, int n){
  static uint32_t tri[2048]; int t = triangles(adj, n, tri);
  for (int i = 0; i < t; i++) for (int j = i + 1; j < t; j++) if (!(tri[i] & tri[j])) return 0;
  return 1;
}
static int tau1(const uint32_t *adj, int n){        /* some vertex meets all cyclic triangles */
  static uint32_t tri[2048]; int t = triangles(adj, n, tri);
  uint32_t common = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
  for (int i = 0; i < t; i++) common &= tri[i];
  return t == 0 || common != 0;
}
/* odd cycles: DP over subsets, paths from the least vertex; n <= 18 */
static u64 odd_cycles(const uint32_t *adj, int n){
  u64 total = 0; static u64 dp[1 << 18];
  for (int s = 0; s < n; s++){
    /* vertices > s only; index subsets of {s+1..n-1} plus endpoint: use dp over (S, v) packed */
    int m = n - s - 1; if (m < 2) continue;
    /* dp2[S][v]: S subset of the m upper vertices, v in S, paths s -> ... -> v covering S */
    u64 *dp2 = calloc((size_t)(1u << m) * m, sizeof(u64));
    for (int v = 0; v < m; v++) if (adj[s] >> (s + 1 + v) & 1) dp2[(size_t)(1u << v) * m + v] = 1;
    for (uint32_t S = 1; S < (1u << m); S++) for (int v = 0; v < m; v++){
      u64 c = dp2[(size_t)S * m + v]; if (!c) continue;
      int sz = __builtin_popcount(S);
      if ((sz % 2 == 0) && (adj[s + 1 + v] >> s & 1)) total += c;     /* cycle length sz + 1 odd */
      uint32_t nb = (adj[s + 1 + v] >> (s + 1)) & ~S & ((1u << m) - 1);
      while (nb){ int u = __builtin_ctz(nb); nb &= nb - 1; dp2[(size_t)(S | (1u << u)) * m + u] += c; }
    }
    free(dp2);
  }
  (void)dp;
  return total;
}
static u64 hpaths(const uint32_t *adj, int n){
  size_t F = (size_t)1 << n; u64 *f = calloc(F * n, sizeof(u64));
  for (int v = 0; v < n; v++) f[((size_t)1 << v) * n + v] = 1;
  for (size_t S = 1; S < F; S++) for (int v = 0; v < n; v++){ u64 x = f[S * n + v]; if (!x) continue;
    uint32_t o = adj[v] & ~(uint32_t)S; while (o){ int u = __builtin_ctz(o); o &= o - 1; f[(S | ((size_t)1 << u)) * n + u] += x; } }
  u64 h = 0; for (int v = 0; v < n; v++) h += f[(F - 1) * n + v]; free(f); return h;
}
typedef struct { u128 *key; size_t cap, cnt; } hset;
static int hs_insert(hset *h, u128 k){
  u64 lo = (u64)k, hi = (u64)(k >> 64); size_t i = ((lo ^ (hi * 0x9E3779B97F4A7C15ULL)) * 0x9E3779B97F4A7C15ULL) >> 17; i &= h->cap - 1;
  while (h->key[i] != (u128)0 || ((u64)h->key[i] == 0 && 0)){ if (h->key[i] == k) return 0; i = (i + 1) & (h->cap - 1); }
  h->key[i] = k; h->cnt++; return 1;
}
static void print_hex(u128 c){ char buf[40]; int i = 0; if (!c){ printf("0"); return; } while (c){ int d = (int)(c & 15); buf[i++] = "0123456789abcdef"[d]; c >>= 4; } while (i) putchar(buf[--i]); }
int main(int argc, char **argv){
  int NMAX = argc > 1 ? atoi(argv[1]) : 12; int with_h = argc > 2 ? atoi(argv[2]) : 13;   /* check H = 1 + 2 alpha_1 up to this order */
  if (NMAX > MAXN) NMAX = MAXN;
  /* order 3: the cyclic triangle; codes never 0 for strong classes, so 0 can mark empty hash slots */
  uint32_t c3[3] = { 1u << 1, 1u << 2, 1u << 0 };
  u128 *reps = malloc(sizeof(u128)); size_t nreps = 1; reps[0] = canon(c3, 3);
  for (int n = 3; n <= NMAX; n++){
    if (n > 3){
      hset h; h.cap = 1u << 16; while (h.cap < 64 * nreps) h.cap <<= 1; h.key = calloc(h.cap, sizeof(u128)); h.cnt = 0;
      u128 *nr = NULL; size_t nnr = 0, cap = 0;
      for (size_t r = 0; r < nreps; r++){ uint32_t adj[MAXN]; decode(reps[r], n - 1, adj);
        /* cover[p]: vertices q such that {p,q} meets every cyclic triangle of the order-(n-1) class (which is strong
           and intercyclic). The extension by z with out-set P is intercyclic iff every new triangle z->p->q->z
           (p in P, q not in P, p->q) has {p,q} a cover; it is strong iff P is neither empty nor everything. */
        static uint32_t tri[4096]; int t = triangles(adj, n - 1, tri); uint32_t cover[MAXN];
        for (int p = 0; p < n - 1; p++){ cover[p] = 0; for (int q = 0; q < n - 1; q++){ if (q == p) continue; uint32_t pq = (1u << p) | (1u << q); int ok = 1; for (int i = 0; i < t; i++) if (!(tri[i] & pq)){ ok = 0; break; } if (ok) cover[p] |= 1u << q; } }
        uint32_t fullm = (1u << (n - 1)) - 1;
        for (uint32_t P = 1; P < fullm; P++){
          int ok = 1; uint32_t notP = fullm & ~P;
          for (uint32_t Q = P; Q; Q &= Q - 1){ int p = __builtin_ctz(Q); if ((adj[p] & notP) & ~cover[p]){ ok = 0; break; } }
          if (!ok) continue;
          uint32_t a2[MAXN]; for (int v = 0; v < n - 1; v++) a2[v] = adj[v]; a2[n - 1] = 0;
          for (int v = 0; v < n - 1; v++){ if (P >> v & 1) a2[n - 1] |= 1u << v; else a2[v] |= 1u << (n - 1); }
#ifdef CHECK
          if (!strongq(a2, n) || !intercyclic(a2, n)){ fprintf(stderr, "fast test disagrees\n"); exit(4); }
#endif
          u128 c = canon(a2, n);
          if (hs_insert(&h, c)){ if (nnr == cap){ cap = cap ? 2 * cap : 1024; nr = realloc(nr, cap * sizeof(u128)); } nr[nnr++] = c; }
        } }
      free(h.key); free(reps); reps = nr; nreps = nnr;
    }
    u64 amin = ~0ULL; int ntau2 = 0;
    for (size_t r = 0; r < nreps; r++){
      uint32_t adj[MAXN]; decode(reps[r], n, adj);
      u64 a1 = odd_cycles(adj, n);
      u64 H = (with_h && n <= with_h) ? hpaths(adj, n) : 0;
      if (with_h && n <= with_h && H != 1 + 2 * a1){ fprintf(stderr, "OCF MISMATCH at n=%d\n", n); return 3; }
      int t = tau1(adj, n) ? 1 : 2; if (t == 2) ntau2++;
      uint32_t op[MAXN]; for (int v = 0; v < n; v++) op[v] = ((1u << n) - 1) & ~adj[v] & ~(1u << v);
      int sc = canon(op, n) == reps[r];
      printf("%d %llu %llu %d %d ", n, a1, H, t, sc); print_hex(reps[r]); printf("\n");
      if (a1 < amin) amin = a1;
    }
    fflush(stdout);
    fprintf(stderr, "order %d: %zu strong intercyclic classes (tau=2: %d), min alpha_1 = %llu\n", n, nreps, ntau2, amin);
  }
  return 0;
}
