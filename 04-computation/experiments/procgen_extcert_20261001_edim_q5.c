/* Orchestrator's independent audit of the Q_5 claims in Allikvere, "The edge multiset dimension of hypercubes"
 * (arXiv:2608.09983; certificate package Zenodo 10.5281/zenodo.21739363). External certificate audit, 2026-10-01.
 * Written from the definitions; the package's code was not read or run (only its result files are compared, by the
 * driver procgen_extcert_20261001_edim_paper_check.py).
 *
 * Definitions: for an edge uv of Q_5 and a vertex s, d(uv,s) = min(d(u,s), d(v,s)); S is edge-multiset resolving iff the
 * histograms h_e = (#{s in S : d(e,s) = r})_r are pairwise distinct over the 80 edges. Complements: h^{V-S}_e(r) =
 * h^V_e(r) - h^S_e(r) and h^V_e does not depend on e (Q_5 is edge-transitive), so S and V - S collide on the same pairs.
 *
 * mode "orbits" (default):
 *   every S with |S| <= 16 is enumerated up to Aut(Q_5) (3840 = 5! 2^5 elements); the representative system is checked
 *   two ways: (i) its size equals the Burnside count, (ii) the orbit sizes 3840/|Stab(S)| sum to C(32,a);
 *   - full test: no representative resolves the 80 edges  =>  edim_m(Q_5) = infinity;
 *   - directional test (the paper's filter): for each direction i the 16 edges {x, x+e_i} have pairwise distinct
 *     histograms; counts the directional survivors over all 2^32 subsets (complements double |S| < 16) and the survivor
 *     orbits under Aut(Q_5) x complement; prints each survivor orbit as "SURV a canon [canon of complement if a = 16]".
 * mode "r4":
 *   counts the weight functions w : V(Q_4) -> {0,1,2} whose vertex histograms (sum_s w(s) [d(x,s) = r])_r are pairwise
 *   distinct over the 16 vertices x (the paper's R4 table; directional survival of S is R4-membership of the 5
 *   projections of S, since d({x, x+e_i}, s) = d(pi_i x, pi_i s)).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned long long u64;
static uint32_t PERM[3840][32];     /* Aut(Q5) as vertex permutations x -> pi(x) ^ t */
static uint32_t TAB[120][4][256];   /* coordinate-permutation tables for 32-bit masks */
static int perms5[120][5];

static void build_aut5(void) {
  int n = 0;
  for (int a = 0; a < 5; a++) for (int b = 0; b < 5; b++) for (int c = 0; c < 5; c++) for (int d = 0; d < 5; d++) for (int e = 0; e < 5; e++) {
    int q[5] = {a, b, c, d, e}; int used = 0, ok = 1; for (int i = 0; i < 5; i++) { if (used >> q[i] & 1) ok = 0; used |= 1 << q[i]; }
    if (!ok) continue; memcpy(perms5[n], q, sizeof q); n++;
  }
  if (n != 120) exit(1);
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
static inline uint32_t pmap(uint32_t m, int pi) {
  return TAB[pi][0][m & 255] | TAB[pi][1][(m >> 8) & 255] | TAB[pi][2][(m >> 16) & 255] | TAB[pi][3][m >> 24];
}
static uint32_t canon5(uint32_t m) {
  uint32_t best = 0xFFFFFFFFu;
  for (int pi = 0; pi < 120; pi++) { uint32_t pm = pmap(m, pi); for (int t = 0; t < 32; t++) { uint32_t q = xlate(pm, t); if (q < best) best = q; } }
  return best;
}
static int stab5(uint32_t m) {
  int s = 0;
  for (int pi = 0; pi < 120; pi++) { uint32_t pm = pmap(m, pi); for (int t = 0; t < 32; t++) if (xlate(pm, t) == m) s++; }
  return s;
}
/* open-addressing hash set of uint32 (keys stored +1; key 0xFFFFFFFF never occurs for |S| <= 16) */
typedef struct { uint32_t *k; size_t cap, n; } hset;
static void hs_init(hset *h, size_t cap) { h->cap = cap; h->n = 0; h->k = calloc(cap, 4); }
static int hs_add(hset *h, uint32_t key) {
  uint32_t kk = key + 1;
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
static void burnside(u64 *orb) {
  u64 tot[33]; memset(tot, 0, sizeof tot);
  for (int g = 0; g < 3840; g++) {
    int seen[32] = {0}; u64 poly[33]; memset(poly, 0, sizeof poly); poly[0] = 1;
    for (int x = 0; x < 32; x++) if (!seen[x]) {
      int L = 0, y = x; while (!seen[y]) { seen[y] = 1; y = PERM[g][y]; L++; }
      for (int d = 32; d >= L; d--) poly[d] += poly[d - L];
    }
    for (int a = 0; a <= 32; a++) tot[a] += poly[a];
  }
  for (int a = 0; a <= 32; a++) orb[a] = tot[a] / 3840;
}
static u64 binom(int n, int k) { u64 r = 1; for (int i = 1; i <= k; i++) r = r * (n - k + i) / i; return r; }

static int EU5[80], EV5[80];
static inline uint32_t hist(int e, uint32_t S) {
  uint32_t h = 0;
  while (S) { int s = __builtin_ctz(S); S &= S - 1;
    int du = __builtin_popcount(EU5[e] ^ s), dv = __builtin_popcount(EV5[e] ^ s);
    h += 1u << (5 * (du < dv ? du : dv)); }
  return h;
}
static int resolves5(uint32_t S) {
  uint32_t key[80];
  for (int e = 0; e < 80; e++) key[e] = hist(e, S);
  for (int i = 0; i < 80; i++) for (int j = i + 1; j < 80; j++) if (key[i] == key[j]) return 0;
  return 1;
}
static int DIRE[5][16]; /* edge indices of direction i */
static int dirsurv5(uint32_t S) {
  for (int i = 0; i < 5; i++) {
    uint32_t key[16];
    for (int k = 0; k < 16; k++) key[k] = hist(DIRE[i][k], S);
    for (int a = 0; a < 16; a++) for (int b = a + 1; b < 16; b++) if (key[a] == key[b]) return 0;
  }
  return 1;
}

static int mode_orbits(void) {
  build_aut5();
  int ne = 0, nd[5] = {0};
  for (int u = 0; u < 32; u++) for (int i = 0; i < 5; i++) { int v = u ^ (1 << i); if (u < v) { EU5[ne] = u; EV5[ne] = v; DIRE[i][nd[i]++] = ne; ne++; } }
  u64 orb[33]; burnside(orb);
  gen_reps(16);
  int ok = 1; u64 total = 0, found = 0, surv_all = 0, surv_orbits_ac = 0, surv16 = 0, selfc16 = 0;
  for (int a = 0; a <= 16; a++) {
    u64 osum = 0, sorb = 0, ssets = 0;
    for (size_t r = 0; r < NREPS[a]; r++) {
      uint32_t S = REPS[a][r];
      u64 osz = 3840 / stab5(S);
      osum += osz; total++;
      if (resolves5(S)) { found++; printf("RESOLVING a=%d mask=%08x\n", a, S); }
      if (dirsurv5(S)) {
        sorb++; ssets += osz;
        if (a < 16) printf("SURV %d %08x\n", a, S);
        else { uint32_t cc = canon5(~S); printf("SURV %d %08x %08x\n", a, S, cc); surv16++; if (cc == S) selfc16++; }
      }
    }
    if (NREPS[a] != orb[a] || osum != binom(32, a)) ok = 0;
    surv_all += (a < 16 ? 2 : 1) * ssets;
    if (a < 16) surv_orbits_ac += sorb;
    printf("a=%d reps=%zu burnside=%llu orbit_sizes_sum=%llu C(32,a)=%llu survivor_orbits=%llu survivor_sets=%llu\n",
           a, NREPS[a], orb[a], osum, binom(32, a), sorb, ssets);
  }
  surv_orbits_ac += (surv16 + selfc16) / 2;
  printf("Q5: %llu orbit representatives with |S| <= 16 checked, resolving: %llu; reps match Burnside and orbit sizes sum to C(32,a): %s\n",
         total, found, ok ? "yes" : "NO");
  printf("Q5: directional survivors among all 2^32 subsets: %llu; survivor orbits under Aut(Q5) x complement: %llu (size-16 Aut-orbits %llu, self-complementary %llu)\n",
         surv_all, surv_orbits_ac, surv16, selfc16);
  return (ok && found == 0) ? 0 : 1;
}

static int mode_r4(void) {
  int D4[16][16];
  for (int x = 0; x < 16; x++) for (int s = 0; s < 16; s++) D4[x][s] = __builtin_popcount(x ^ s);
  int w[16] = {0}; uint32_t key[16] = {0};
  u64 good = 0, total = 0;
  for (;;) {
    total++;
    int distinct = 1;
    for (int a = 0; a < 16 && distinct; a++) for (int b = a + 1; b < 16; b++) if (key[a] == key[b]) { distinct = 0; break; }
    good += distinct;
    int s = 0;
    while (s < 16 && w[s] == 2) { w[s] = 0; for (int x = 0; x < 16; x++) key[x] -= 2u << (5 * D4[x][s]); s++; }
    if (s == 16) break;
    w[s]++; for (int x = 0; x < 16; x++) key[x] += 1u << (5 * D4[x][s]);
  }
  printf("R4: %llu resolving weight functions out of %llu\n", good, total);
  return 0;
}

int main(int argc, char **argv) {
  if (argc > 1 && !strcmp(argv[1], "r4")) return mode_r4();
  return mode_orbits();
}
