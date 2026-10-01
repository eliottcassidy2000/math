/* Orchestrator's independent audit of the edim lane (2026-10-01): edge multiset dimension of Q_6.
 * Written from the definition (Allikvere, arXiv:2608.09983): for an edge uv and a vertex s,
 * d(uv,s) = min(d(u,s), d(v,s)); S is edge-multiset resolving iff the multisets {d(e,s) : s in S}
 * are pairwise distinct over all 192 edges of Q_6. The lane's code was not read.
 *
 * Reduction: split Q_6 by coordinate 5 (bit 5): S = A x {0} + B x {1}, A, B subsets of Q_5.
 * Flipping bit 5 swaps the layers, so WLOG |A| >= |B|; Aut(Q_5) (3840 elements: coordinate
 * permutations and translations of bits 0..4) acts on both layers at once, so WLOG A is the
 * minimum of its Aut(Q_5)-orbit. Every B of size k-|A| is enumerated (no further reduction).
 *
 * usage: ea KMIN KMAX FOUNDFILE   (prints per-k leaf counts; writes every resolving set found)
 * Driver: procgen_edim_20261001_orchestrator_check.py
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned long long u64;
static int EU[192], EV[192];
static uint32_t contrib[64][192];   /* 1 << (5*d(e,s)) */
static uint32_t PERM[3840][32];     /* Aut(Q5) as vertex permutations */
static uint32_t TAB[120][4][256];   /* coordinate-permutation tables for 32-bit masks */
static int NPERM5 = 0;

static int popc(unsigned x) { return __builtin_popcount(x); }

static void build_q6(void) {
  int ne = 0;
  for (int u = 0; u < 64; u++) for (int i = 0; i < 6; i++) { int v = u ^ (1 << i); if (u < v) { EU[ne] = u; EV[ne] = v; ne++; } }
  if (ne != 192) { fprintf(stderr, "edges %d\n", ne); exit(1); }
  for (int s = 0; s < 64; s++) for (int e = 0; e < 192; e++) {
    int du = popc(EU[e] ^ s), dv = popc(EV[e] ^ s); int d = du < dv ? du : dv;
    contrib[s][e] = 1u << (5 * d);
  }
}
/* Aut(Q5): x -> pi(x) ^ t, pi a coordinate permutation */
static int perms5[120][5];
static void build_aut5(void) {
  int p[5] = {0, 1, 2, 3, 4}; int n = 0;
  /* all permutations by Heap-free brute force */
  for (int a = 0; a < 5; a++) for (int b = 0; b < 5; b++) for (int c = 0; c < 5; c++) for (int d = 0; d < 5; d++) for (int e = 0; e < 5; e++) {
    int q[5] = {a, b, c, d, e}; int used = 0, ok = 1; for (int i = 0; i < 5; i++) { if (used >> q[i] & 1) ok = 0; used |= 1 << q[i]; }
    if (!ok) continue; memcpy(perms5[n], q, sizeof q); n++;
  }
  (void)p; NPERM5 = n; if (n != 120) exit(1);
  int g = 0;
  for (int pi = 0; pi < 120; pi++) for (int t = 0; t < 32; t++, g++) for (int x = 0; x < 32; x++) {
    int y = 0; for (int i = 0; i < 5; i++) if (x >> i & 1) y |= 1 << perms5[pi][i];
    PERM[g][x] = y ^ t;
  }
  for (int pi = 0; pi < 120; pi++) for (int bi = 0; bi < 4; bi++) for (int by = 0; by < 256; by++) {
    uint32_t m = 0; for (int j = 0; j < 8; j++) if (by >> j & 1) { int x = bi * 8 + j; int y = 0; for (int i = 0; i < 5; i++) if (x >> i & 1) y |= 1 << perms5[pi][i]; m |= 1u << y; }
    TAB[pi][bi][by] = m;
  }
}
static inline uint32_t xlate(uint32_t m, int t) { /* bit x -> bit x^t */
  if (t & 1) m = ((m & 0x55555555u) << 1) | ((m >> 1) & 0x55555555u);
  if (t & 2) m = ((m & 0x33333333u) << 2) | ((m >> 2) & 0x33333333u);
  if (t & 4) m = ((m & 0x0F0F0F0Fu) << 4) | ((m >> 4) & 0x0F0F0F0Fu);
  if (t & 8) m = ((m & 0x00FF00FFu) << 8) | ((m >> 8) & 0x00FF00FFu);
  if (t & 16) m = (m << 16) | (m >> 16);
  return m;
}
static uint32_t canon5(uint32_t m) {
  uint32_t best = 0xFFFFFFFFu;
  for (int pi = 0; pi < 120; pi++) {
    uint32_t pm = TAB[pi][0][m & 255] | TAB[pi][1][(m >> 8) & 255] | TAB[pi][2][(m >> 16) & 255] | TAB[pi][3][m >> 24];
    for (int t = 0; t < 32; t++) { uint32_t q = xlate(pm, t); if (q < best) best = q; }
  }
  return best;
}
/* simple open-addressing hash set of uint32 (nonzero keys stored +1 offset) */
typedef struct { uint32_t *k; size_t cap, n; } hset;
static void hs_init(hset *h, size_t cap) { h->cap = cap; h->n = 0; h->k = calloc(cap, 4); }
static int hs_add(hset *h, uint32_t key) {
  uint32_t kk = key + 1; /* key 0xFFFFFFFF never occurs for a < 32 */
  size_t i = (kk * 2654435761u) % h->cap;
  while (h->k[i]) { if (h->k[i] == kk) return 0; i = (i + 1) % h->cap; }
  h->k[i] = kk; h->n++; return 1;
}
static uint32_t *REPS[33]; static size_t NREPS[33];
static void gen_reps(int amax) {
  REPS[0] = malloc(4); REPS[0][0] = 0; NREPS[0] = 1;
  for (int a = 1; a <= amax; a++) {
    hset h; hs_init(&h, 4 * 1000000);
    size_t cap = 1024; uint32_t *out = malloc(cap * 4); size_t n = 0;
    for (size_t r = 0; r < NREPS[a - 1]; r++) {
      uint32_t R = REPS[a - 1][r];
      for (int x = 0; x < 32; x++) if (!(R >> x & 1)) {
        uint32_t c = canon5(R | 1u << x);
        if (hs_add(&h, c)) { if (n == cap) { cap *= 2; out = realloc(out, cap * 4); } out[n++] = c; }
      }
    }
    REPS[a] = out; NREPS[a] = n; free(h.k);
  }
}
/* Burnside count of a-subset orbits of Q5 under Aut(Q5), from the cycle types of the 3840 permutations */
static void burnside(int amax, u64 *orb) {
  /* polynomial coefficients as unsigned 128? counts small enough for u64: C(32,16)*3840 < 2^64 */
  u64 tot[33]; memset(tot, 0, sizeof tot);
  for (int g = 0; g < 3840; g++) {
    int seen[32] = {0}; u64 poly[33]; memset(poly, 0, sizeof poly); poly[0] = 1;
    for (int x = 0; x < 32; x++) if (!seen[x]) {
      int L = 0, y = x; while (!seen[y]) { seen[y] = 1; y = PERM[g][y]; L++; }
      for (int d = 32; d >= L; d--) poly[d] += poly[d - L];
    }
    for (int a = 0; a <= 32; a++) tot[a] += poly[a];
  }
  for (int a = 0; a <= amax; a++) orb[a] = tot[a] / 3840;
}

/* ------------------------------------------------------------------ search */
static uint32_t cur[16][192];   /* cur[level][e] */
static int hotE[512], hotF[512], nhot = 0;
static u64 leaves, fullchecks, found;
static int Bsel[16];
static uint32_t Acur;
static int kcur;
static FILE *foundf;

static int full_check(const uint32_t *v, int *ce, int *cf) {
  /* distinctness of 192 values by hashing */
  static uint32_t key[1024]; static int idx[1024]; static uint32_t stamp[1024]; static uint32_t gen = 0;
  gen++;
  for (int e = 0; e < 192; e++) {
    uint32_t x = v[e]; unsigned i = (x * 2654435761u) >> 22;
    while (stamp[i] == gen) { if (key[i] == x) { *ce = idx[i]; *cf = e; return 0; } i = (i + 1) & 1023; }
    stamp[i] = gen; key[i] = x; idx[i] = e;
  }
  return 1;
}
static void leaf(const uint32_t *prev, int last) {
  leaves++;
  const uint32_t *c = contrib[last];
  for (int h = 0; h < nhot; h++) {
    int e = hotE[h], f = hotF[h];
    if (prev[e] + c[e] == prev[f] + c[f]) {
      if (h > 0) { int te = hotE[h - 1], tf = hotF[h - 1]; hotE[h - 1] = e; hotF[h - 1] = f; hotE[h] = te; hotF[h] = tf; }
      return;
    }
  }
  fullchecks++;
  uint32_t v[192]; for (int e = 0; e < 192; e++) v[e] = prev[e] + c[e];
  int ce, cf;
  if (full_check(v, &ce, &cf)) {
    found++;
    uint64_t S = (uint64_t)Acur; for (int i = 0; i < kcur - popc(Acur) - 1; i++) S |= 1ull << (32 + Bsel[i]); S |= 1ull << last;
    fprintf(foundf, "k=%d %016llx\n", kcur, (unsigned long long)S);
  } else {
    if (nhot < 512) { hotE[nhot] = ce; hotF[nhot] = cf; nhot++; }
    else { hotE[511] = ce; hotF[511] = cf; }
  }
}
/* choose b more layer-1 vertices from x >= start; level = number already chosen */
static void dfsB(int level, int b, int start) {
  if (b == 0) return;
  if (b == 1) {
    for (int x = start; x < 32; x++) leaf(cur[level], 32 + x);
    return;
  }
  for (int x = start; x <= 32 - b; x++) {
    const uint32_t *c = contrib[32 + x];
    for (int e = 0; e < 192; e++) cur[level + 1][e] = cur[level][e] + c[e];
    Bsel[level] = x;
    dfsB(level + 1, b - 1, x + 1);
  }
}
static void search_k(int k) {
  u64 L0 = leaves, F0 = found, FC0 = fullchecks; kcur = k;
  for (int a = (k + 1) / 2; a <= k && a <= 32; a++) {
    int b = k - a;
    for (size_t r = 0; r < NREPS[a]; r++) {
      uint32_t A = REPS[a][r]; Acur = A;
      for (int e = 0; e < 192; e++) { uint32_t s = 0; for (int x = 0; x < 32; x++) if (A >> x & 1) s += contrib[x][e]; cur[0][e] = s; }
      if (b == 0) {
        /* the set is A alone: check directly */
        leaves++; int ce, cf; if (full_check(cur[0], &ce, &cf)) { found++; fprintf(foundf, "k=%d %016llx\n", k, (unsigned long long)A); }
      } else dfsB(0, b, 0);
    }
  }
  printf("k=%d leaves=%llu fullchecks=%llu resolving_found=%llu\n", k, leaves - L0, fullchecks - FC0, found - F0);
  fflush(stdout);
}

int main(int argc, char **argv) {
  int kmin = atoi(argv[1]), kmax = atoi(argv[2]);
  build_q6(); build_aut5();
  u64 orb[33]; burnside(32, orb);
  gen_reps(kmax);
  int ok = 1;
  for (int a = 0; a <= kmax; a++) { if (NREPS[a] != orb[a]) ok = 0; printf("a=%d reps=%zu burnside=%llu\n", a, NREPS[a], orb[a]); }
  printf("orbit representatives %s\n", ok ? "MATCH Burnside" : "MISMATCH");
  if (!ok) return 2;
  foundf = fopen(argc > 3 ? argv[3] : "ea_found.txt", "w");
  for (int k = kmin; k <= kmax; k++) {
    /* expected leaves = sum_a #reps(a) * C(32, k-a) */
    u64 exp = 0; for (int a = (k + 1) / 2; a <= k; a++) { u64 C = 1; int b = k - a; for (int i = 0; i < b; i++) C = C * (32 - i) / (i + 1); exp += NREPS[a] * C; }
    printf("k=%d expected leaves=%llu\n", k, exp);
    search_k(k);
  }
  fclose(foundf);
  return 0;
}
