/* few_cycle_strong.c - opus-2026-11-08-S86.
 * All strong tournaments with at most M directed odd cycles (alpha_1 <= M), up to isomorphism, with their full
 * odd-cycle packing vectors (alpha_1, ..., alpha_5).
 * Completeness: alpha_1 cannot increase under vertex deletion (cycles of T - v are cycles of T), and every strong
 * tournament of order >= 4 has a vertex whose deletion leaves it strong (Moon). So the classes of order k+1 are
 * one-vertex extensions of the classes of order k, and level-by-level extension with isomorphism rejection is
 * exhaustive. Moon's sharp cycle count (at least k - l + 1 cycles of each length l) gives alpha_1 >= m(k) =
 * sum over odd l in [3,k] of (k - l + 1), so the levels stop at the first k with m(k) > M.
 * Pruning: an extension by z with out-set P adds at least t3 = #{(p,q): p in P, q not in P, p -> q} odd cycles.
 * Canonical forms: as in intercyclic_gen.c (kind-pasteur's classgen.c refinement), with 192-bit codes (n <= 20).
 * usage: few_cycle_strong M > out.txt     (line: n a1 a2 a3 a4 a5 H code_hex)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
#define MAXN 20
/* canonical codes need n(n-1)/2 bits: 190 at n = 20, so codes are 3 x 64-bit words (an earlier 128-bit version
   silently truncated at order 17) */
typedef struct { u64 w[3]; } code_t;
static int code_lt(code_t a, code_t b){ for (int i = 2; i >= 0; i--){ if (a.w[i] != b.w[i]) return a.w[i] < b.w[i]; } return 0; }
static int code_eq(code_t a, code_t b){ return a.w[0] == b.w[0] && a.w[1] == b.w[1] && a.w[2] == b.w[2]; }
static int code_zero(code_t a){ return !a.w[0] && !a.w[1] && !a.w[2]; }
static int code_bit(code_t a, int b){ return (int)(a.w[b >> 6] >> (b & 63) & 1); }
static int nn; static uint32_t A[MAXN]; static code_t best; static int have;
static code_t code_of(const int *pos){ int inv[MAXN]; for (int v = 0; v < nn; v++) inv[pos[v]] = v; code_t c = {{0, 0, 0}}; int b = 0;
  for (int i = 0; i < nn; i++) for (int j = i + 1; j < nn; j++){ if (A[inv[i]] >> inv[j] & 1) c.w[b >> 6] |= 1ULL << (b & 63); b++; } return c; }
static int ncells(const int *cell){ int m = -1; for (int v = 0; v < nn; v++) if (cell[v] > m) m = cell[v]; return m + 1; }
static void refine(int *cell){ for (;;){ int k = ncells(cell); static int sig[MAXN][MAXN + 1];
    for (int v = 0; v < nn; v++){ sig[v][0] = cell[v]; for (int c = 0; c < k; c++) sig[v][c + 1] = 0; for (int u = 0; u < nn; u++) if (A[v] >> u & 1) sig[v][cell[u] + 1]++; }
    int ord[MAXN]; for (int v = 0; v < nn; v++) ord[v] = v;
    for (int i = 1; i < nn; i++){ int x = ord[i], j = i - 1; while (j >= 0 && memcmp(sig[ord[j]], sig[x], (k + 1) * sizeof(int)) > 0){ ord[j + 1] = ord[j]; j--; } ord[j + 1] = x; }
    int nc[MAXN]; int c = 0; nc[ord[0]] = 0;
    for (int i = 1; i < nn; i++){ if (memcmp(sig[ord[i]], sig[ord[i - 1]], (k + 1) * sizeof(int))) c++; nc[ord[i]] = c; }
    int changed = (c + 1 != k); memcpy(cell, nc, sizeof(int) * nn); if (!changed) return; } }
static void search(int *cell){ refine(cell); int k = ncells(cell);
  if (k == nn){ code_t c = code_of(cell); if (!have || code_lt(c, best)){ best = c; have = 1; } return; }
  int cnt[MAXN] = {0}; for (int v = 0; v < nn; v++) cnt[cell[v]]++; int target = -1; for (int c = 0; c < k; c++) if (cnt[c] > 1){ target = c; break; }
  for (int v = 0; v < nn; v++) if (cell[v] == target){ int c2[MAXN]; for (int u = 0; u < nn; u++){ c2[u] = cell[u] > target ? cell[u] + 1 : cell[u]; if (cell[u] == target && u != v) c2[u] = target + 1; } search(c2); } }
static code_t canon(const uint32_t *adj, int n){ nn = n; memcpy(A, adj, sizeof(uint32_t) * n); have = 0; int cell[MAXN]; for (int v = 0; v < n; v++) cell[v] = 0; search(cell); return best; }
static void decode(code_t c, int n, uint32_t *adj){ memset(adj, 0, sizeof(uint32_t) * n); int b = 0;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++){ if (code_bit(c, b)) adj[i] |= 1u << j; else adj[j] |= 1u << i; b++; } }
/* odd cycles through vertex z (paths z -> ... -> z), DP over subsets of the other vertices */
static u64 odd_through(const uint32_t *adj, int n, int z, u64 cap){
  int m = n - 1; int idx[MAXN], k = 0; for (int v = 0; v < n; v++) if (v != z) idx[k++] = v;
  size_t F = (size_t)1 << m; u64 *dp = calloc(F * m, sizeof(u64)); u64 tot = 0;
  for (int i = 0; i < m; i++) if (adj[z] >> idx[i] & 1) dp[((size_t)1 << i) * m + i] = 1;
  for (size_t S = 1; S < F && tot <= cap; S++) for (int i = 0; i < m; i++){ u64 c = dp[S * m + i]; if (!c) continue;
    if ((__builtin_popcountll(S) % 2 == 0) && (adj[idx[i]] >> z & 1)) tot += c;
    for (int j = 0; j < m; j++) if (!(S >> j & 1) && (adj[idx[i]] >> idx[j] & 1)) dp[(S | ((size_t)1 << j)) * m + j] += c; }
  free(dp); return tot;
}
/* full packing vector by subset DP: cyc[S] = directed odd cycles on vertex set S */
static void packing(const uint32_t *adj, int n, u64 *alpha){
  size_t F = (size_t)1 << n; u64 *cyc = calloc(F, sizeof(u64)); u64 *dp = calloc(F * n, sizeof(u64));
  for (int s = 0; s < n; s++) dp[((size_t)1 << s) * n + s] = 1;
  for (size_t S = 1; S < F; S++){ int s = __builtin_ctzll(S);
    for (int v = 0; v < n; v++){ u64 c = dp[S * n + v]; if (!c) continue;
      uint32_t nb = adj[v] & ~(uint32_t)S & ~((1u << (s + 1)) - 1);
      while (nb){ int u = __builtin_ctz(nb); nb &= nb - 1; dp[(S | ((size_t)1 << u)) * n + u] += c; } }
    int k = __builtin_popcountll(S);
    if (k >= 3 && (k & 1)){ u64 t = 0; for (int v = 0; v < n; v++) if (dp[S * n + v] && (adj[v] >> s & 1)) t += dp[S * n + v]; cyc[S] = t; } }
  free(dp);
  u64 (*P)[7] = calloc(F, sizeof *P); P[0][0] = 1;
  for (size_t U = 1; U < F; U++){ int m = __builtin_ctzll(U); size_t rest = U & ~((size_t)1 << m);
    for (size_t T = rest;; T = (T - 1) & rest){ size_t S = T | ((size_t)1 << m);
      if (cyc[S]) for (int k = 0; k < 6; k++) if (P[U & ~S][k]) P[U][k + 1] += cyc[S] * P[U & ~S][k];
      if (!T) break; } }
  for (int k = 0; k < 7; k++){ alpha[k] = 0; for (size_t U = 0; U < F; U++) alpha[k] += P[U][k]; }
  free(P); free(cyc);
}
typedef struct { code_t *key; size_t cap, cnt; } hset;
static int hs_insert(hset *h, code_t k){ size_t i = ((k.w[0] ^ (k.w[1] * 0x9E3779B97F4A7C15ULL) ^ (k.w[2] * 0xC2B2AE3D27D4EB4FULL)) * 0x9E3779B97F4A7C15ULL) >> 17; i &= h->cap - 1;
  while (!code_zero(h->key[i])){ if (code_eq(h->key[i], k)) return 0; i = (i + 1) & (h->cap - 1); } h->key[i] = k; h->cnt++; return 1; }
static void print_hex(code_t c){ printf("%016llx%016llx%016llx", c.w[2], c.w[1], c.w[0]); }
static int moon_floor(int k){ int s = 0; for (int l = 3; l <= k; l += 2) s += k - l + 1; return s; }
int main(int argc, char **argv){
  u64 M = argc > 1 ? strtoull(argv[1], 0, 10) : 40;
  uint32_t c3[3] = { 1u << 1, 1u << 2, 1u << 0 };
  code_t *reps = malloc(sizeof(code_t)); u64 *a1s = malloc(sizeof(u64)); size_t nreps = 1; reps[0] = canon(c3, 3); a1s[0] = 1;
  for (int n = 3; n < MAXN; n++){
    if (moon_floor(n) > (int)M){ fprintf(stderr, "order %d: Moon floor %d > M, stop\n", n, moon_floor(n)); break; }
    if (n > 3){
      hset h; h.cap = 1u << 16; while (h.cap < 64 * nreps + 1024) h.cap <<= 1; h.key = calloc(h.cap, sizeof(code_t)); h.cnt = 0;
      code_t *nr = NULL; u64 *na = NULL; size_t nnr = 0, cap = 0;
      for (size_t r = 0; r < nreps; r++){ uint32_t adj[MAXN]; decode(reps[r], n - 1, adj); u64 base = a1s[r];
        uint32_t fullm = (1u << (n - 1)) - 1;
        for (uint32_t P = 1; P < fullm; P++){                    /* strong iff P neither empty nor everything */
          uint32_t notP = fullm & ~P; u64 t3 = 0;
          for (uint32_t Q = P; Q; Q &= Q - 1){ int p = __builtin_ctz(Q); t3 += __builtin_popcount(adj[p] & notP); }
          if (base + t3 > M) continue;
          uint32_t a2[MAXN]; for (int v = 0; v < n - 1; v++) a2[v] = adj[v]; a2[n - 1] = 0;
          for (int v = 0; v < n - 1; v++){ if (P >> v & 1) a2[n - 1] |= 1u << v; else a2[v] |= 1u << (n - 1); }
          u64 through = odd_through(a2, n, n - 1, M - base);
          if (base + through > M) continue;
          code_t c = canon(a2, n);
          if (hs_insert(&h, c)){ if (nnr == cap){ cap = cap ? 2 * cap : 1024; nr = realloc(nr, cap * sizeof(code_t)); na = realloc(na, cap * sizeof(u64)); } nr[nnr] = c; na[nnr] = base + through; nnr++; }
        } }
      free(h.key); free(reps); free(a1s); reps = nr; a1s = na; nreps = nnr;
    }
    for (size_t r = 0; r < nreps; r++){ uint32_t adj[MAXN]; decode(reps[r], n, adj); u64 al[7]; packing(adj, n, al);
      if (al[1] != a1s[r]){ fprintf(stderr, "alpha_1 mismatch\n"); return 3; }
      u64 H = 0; for (int k = 0; k < 7; k++) H += al[k] << k;
      printf("%d %llu %llu %llu %llu %llu %llu ", n, al[1], al[2], al[3], al[4], al[5], H); print_hex(reps[r]); printf("\n"); }
    fflush(stdout);
    fprintf(stderr, "order %d: %zu strong classes with alpha_1 <= %llu\n", n, nreps, M);
    if (!nreps) break;
  }
  return 0;
}
