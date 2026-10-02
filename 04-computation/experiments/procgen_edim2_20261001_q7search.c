// procgen_edim2_20261001_q7search.c -- exhaustive search for edge-multiset resolving k-sets of Q_D (D = 6 or 7)
// with a cheap (non-canonical, but complete) symmetry normal form.
//
// Normal form (WLOG, every k-set S with k >= 1 is Aut(Q_D)-equivalent to one of these):
//   * coordinate D-1 has maximal imbalance |beta_{D-1}| = max_i |beta_i(S)|, beta_i(S) = #{s_i = 0} - #{s_i = 1},
//     and beta_{D-1} >= 0 (transpose coordinates, flip coordinate D-1); write S = A x {0} u B x {1},
//     |A| = a >= b = |B|, a - b = beta_{D-1} >= |beta_i(S)| for i < D-1;
//   * 0 in A (translate by an element of A on coordinates 0..D-2; |beta_i| unchanged);
//   * column sums of A non-increasing: c_0(A) >= c_1(A) >= ... >= c_{D-2}(A), c_i(A) = #{x in A : x_i = 1}
//     (permute coordinates 0..D-2; the constraints |beta_i(S)| <= a - b are permutation invariant).
// A ranges over all a-subsets of Q_{D-1} with these two properties (with redundancy), B over all b-subsets
// with |beta_i(A) + beta_i(B)| <= a - b (pruned DFS).  Leaf test: key(e) = sum_s 16^{d(e,s)} over levels
// 0..D-2 (level D-1 implied, counts <= 15 since k <= 15); exact duplicate detection.
// Optional strengthening (argv[5] = 1): A is translate-maximal, i.e. for every t in A the sorted column-sum
// vector of A xor t is lexicographically <= that of A (WLOG: translate by a maximising t, then sort).
// Mode tmax = 2 additionally keeps only the canonical representative under the residual group (one A per
// Aut(Q_{D-1})-orbit); with k < a the program only counts these A's (check against Burnside).
// usage: q7search D k a [maxprint] [tmax]   prints RESOLVING lines (vertex lists) and a SUMMARY line.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static int D, NV, H, NE, K, A_, B_, DELTA, MAXPRINT = 20, TMAX = 0;
static uint32_t inc[128][448];
static uint32_t lvlB[16][448];
static uint32_t tkey[2048]; static uint64_t tgen[2048]; static uint64_t gen = 1;
static long long nA = 0, leaves = 0, found = 0;
static int Aset[16], Bset[16];
static int betaA[8], betaB[8];
static int colA[8];
static int check(const uint32_t *key) {
  gen++;
  for (int e = 0; e < NE; e++) {
    uint32_t k = key[e], h = (k * 2654435761u) >> 21;
    while (tgen[h] == gen) { if (tkey[h] == k) return 0; h = (h + 1) & 2047; }
    tgen[h] = gen; tkey[h] = k;
  }
  return 1;
}
static void recB(int depth, int start) {
  int r = B_ - depth;
  for (int i = 0; i < D - 1; i++) {
    int bt = betaA[i] + betaB[i];
    if (bt - r > DELTA || bt + r < -DELTA) return;
  }
  if (r == 0) {
    leaves++;
    if (check(lvlB[depth])) {
      found++;
      if (found <= MAXPRINT) {
        printf("RESOLVING k=%d a=%d set:", K, A_);
        for (int j = 0; j < A_; j++) printf(" %d", Aset[j]);
        for (int j = 0; j < B_; j++) printf(" %d", Bset[j] + H);
        printf("\n"); fflush(stdout);
      }
    }
    return;
  }
  for (int v = start; v <= H - r; v++) {
    uint32_t *src = lvlB[depth], *dst = lvlB[depth + 1], *in = inc[H + v];
    for (int e = 0; e < NE; e++) dst[e] = src[e] + in[e];
    Bset[depth] = v;
    for (int i = 0; i < D - 1; i++) betaB[i] += (v >> i & 1) ? -1 : 1;
    recB(depth + 1, v + 1);
    for (int i = 0; i < D - 1; i++) betaB[i] -= (v >> i & 1) ? -1 : 1;
  }
}
static int cmpdesc(const void *x, const void *y) { return *(const int *)y - *(const int *)x; }
// --- mode 2: canonical form under the residual symmetry group of the filters ---------------------------
// The filters (0 in A, column sums non-increasing, translate-maximal) are met by at least one element of every
// Aut(Q_{D-1})-orbit.  An automorphism g = pi o (x -> x xor t) maps A to a set meeting the filters iff t in A,
// the sorted column-sum vector of A xor t equals that of A, and pi sorts the column sums of A xor t.  A is
// accepted iff its mask is <= the masks of all such images; then exactly one element per orbit is accepted.
static int c2g[8], usedq[8], permg[8], tg;
static uint64_t maskA;
static int rejectflag;
static void permrec(int p) {
  if (rejectflag) return;
  if (p == D - 1) {
    uint64_t m = 0;
    for (int j = 0; j < A_; j++) {
      int y = Aset[j] ^ tg, z = 0;
      for (int q = 0; q < D - 1; q++) z |= ((y >> permg[q]) & 1) << q;
      m |= 1ULL << z;
    }
    if (m < maskA) rejectflag = 1;
    return;
  }
  for (int q = 0; q < D - 1; q++) if (!usedq[q] && c2g[q] == colA[p]) {
    usedq[q] = 1; permg[p] = q; permrec(p + 1); usedq[q] = 0;
    if (rejectflag) return;
  }
}
static int canonical_residual(void) {
  maskA = 0; for (int j = 0; j < A_; j++) maskA |= 1ULL << Aset[j];
  rejectflag = 0;
  for (int j = 0; j < A_ && !rejectflag; j++) {
    tg = Aset[j]; int srt[8];
    for (int i = 0; i < D - 1; i++) srt[i] = c2g[i] = (tg >> i & 1) ? A_ - colA[i] : colA[i];
    qsort(srt, D - 1, sizeof(int), cmpdesc);
    int same = 1; for (int i = 0; i < D - 1; i++) if (srt[i] != colA[i]) { same = 0; break; }
    if (!same) continue;
    for (int q = 0; q < D - 1; q++) usedq[q] = 0;
    permrec(0);
  }
  return !rejectflag;
}
static void leafA(void) {
  // column sums sorted?
  for (int i = 0; i + 1 < D - 1; i++) if (colA[i] < colA[i + 1]) return;
  if (TMAX) {
    // translate-maximality: for every t in A, sorted column sums of A xor t are lexicographically <= those of A
    for (int j = 1; j < A_; j++) {
      int t = Aset[j], c2[8];
      for (int i = 0; i < D - 1; i++) c2[i] = (t >> i & 1) ? A_ - colA[i] : colA[i];
      qsort(c2, D - 1, sizeof(int), cmpdesc);
      for (int i = 0; i < D - 1; i++) { if (c2[i] < colA[i]) break; if (c2[i] > colA[i]) return; }
    }
  }
  if (TMAX >= 2 && !canonical_residual()) return;
  nA++;
  if (B_ < 0) return;
  memset(lvlB[0], 0, sizeof lvlB[0]);
  for (int j = 0; j < A_; j++) { int w = Aset[j]; for (int e = 0; e < NE; e++) lvlB[0][e] += inc[w][e]; }
  for (int i = 0; i < D - 1; i++) { betaA[i] = A_ - 2 * colA[i]; betaB[i] = 0; }
  recB(0, 0);
}
static void recA(int depth, int start) {
  int r = A_ - depth;
  // pruning: final c_i <= c_i + r, final c_{i+1} >= c_{i+1}; need c_i(final) >= c_{i+1}(final)
  for (int i = 0; i + 1 < D - 1; i++) if (colA[i + 1] > colA[i] + r) return;
  if (r == 0) { leafA(); return; }
  for (int v = start; v <= H - r; v++) {
    Aset[depth] = v;
    for (int i = 0; i < D - 1; i++) colA[i] += v >> i & 1;
    recA(depth + 1, v + 1);
    for (int i = 0; i < D - 1; i++) colA[i] -= v >> i & 1;
  }
}
int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: q7search D k a [maxprint]\n"); return 1; }
  D = atoi(argv[1]); K = atoi(argv[2]); A_ = atoi(argv[3]); B_ = K - A_; DELTA = A_ - B_;
  if (argc > 4) MAXPRINT = atoi(argv[4]);
  if (argc > 5) TMAX = atoi(argv[5]);
  if (D < 3 || D > 7 || K < 1 || K > 15 || A_ < 1 || (DELTA < 0 && B_ >= 0)) { fprintf(stderr, "bad args\n"); return 1; }
  if (B_ < 0 && TMAX < 2) { fprintf(stderr, "count-only mode (k < a) needs tmax = 2\n"); return 1; }
  NV = 1 << D; H = NV / 2; NE = D * H;
  int e = 0;
  for (int u = 0; u < NV; u++) for (int i = 0; i < D; i++) if (!(u >> i & 1)) {
    for (int w = 0; w < NV; w++) {
      int dd = __builtin_popcount((u ^ w) & ~(1 << i));
      inc[w][e] = (dd <= D - 2) ? (1u << (4 * dd)) : 0u;
    }
    e++;
  }
  // A must contain vertex 0
  Aset[0] = 0; for (int i = 0; i < D - 1; i++) colA[i] = 0;
  recA(1, 1);
  printf("SUMMARY D=%d k=%d a=%d b=%d tmax=%d Asets=%lld leaves=%lld found=%lld\n", D, K, A_, B_, TMAX, nA, leaves, found);
  return 0;
}
