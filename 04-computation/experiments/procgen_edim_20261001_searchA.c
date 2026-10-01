// procgen_edim_20261001_searchA.c -- METHOD A of the exhaustive search for edge-multiset resolving
// sets of Q6 (lane procgen_edim, 2026-10-01).
// Normal form A: coordinate 5 is a coordinate of MAXIMUM imbalance; after flipping it,
//   S = A x {0}  u  B x {1}  with |A| = a >= b = |B|, a - b = max_i |beta_i(S)|,
//   beta_i(S) = #{s in S: s_i = 0} - #{s in S: s_i = 1};
// A ranges over the Aut(Q5)-orbit representatives of a-subsets of Q5 (file), B over ALL b-subsets
// of Q5 with |beta_i(S)| <= a - b for i = 0..4 (pruned DFS, exact at the leaves).
// Leaf test: key(e) = sum_s 16^{d(e,s)} (levels 0..4, 4 bits each; level 5 implied since the
// levels sum to k <= 15); duplicate keys detected with a 1024-slot linear-probing table.
// usage: searchA k a repfile [maxdef]
//   prints "RESOLVING k= mask=" for every resolving leaf; with maxdef>0 also "NEAR k= def= mask="
//   for leaves with 1 <= defect <= maxdef (defect = 192 - #distinct histograms); final SUMMARY line.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define NE 192
static uint32_t inc[64][NE];
static int K, A_, B_, DELTA, MAXDEF = 0, lastdef = 0;
static uint32_t lvl[16][NE];
static uint32_t tkey[1024]; static uint64_t tgen[1024]; static uint64_t gen = 1;
static long long leaves = 0, found = 0, near = 0, pruned = 0;
static int chosen[16], imb[5];
static uint32_t curA;
static const uint32_t MASK20 = 0xFFFFF;
static int check(const uint32_t *key) {
  int def = 0; gen++;
  for (int e = 0; e < NE; e++) {
    uint32_t k = key[e] & MASK20, h = (k * 2654435761u) >> 22; int dup = 0;
    while (tgen[h] == gen) { if (tkey[h] == k) { dup = 1; break; } h = (h + 1) & 1023; }
    if (dup) { if (++def > MAXDEF) return 0; continue; }
    tgen[h] = gen; tkey[h] = k;
  }
  lastdef = def; return 1;
}
static void rec(int depth, int start) {
  int r = B_ - depth;
  for (int i = 0; i < 5; i++) if (imb[i] - r > DELTA || imb[i] + r < -DELTA) { pruned++; return; }
  if (r == 0) {
    leaves++;
    if (check(lvl[depth])) {
      uint64_t m = curA; for (int j = 0; j < B_; j++) m |= 1ULL << (32 + chosen[j]);
      if (lastdef) { near++; printf("NEAR k=%d def=%d mask=%016llx\n", K, lastdef, (unsigned long long)m); }
      else { found++; printf("RESOLVING k=%d mask=%016llx\n", K, (unsigned long long)m); }
      fflush(stdout);
    }
    return;
  }
  for (int v = start; v <= 32 - r; v++) {
    uint32_t *src = lvl[depth], *dst = lvl[depth + 1], *in = inc[32 + v];
    for (int e = 0; e < NE; e++) dst[e] = src[e] + in[e];
    chosen[depth] = v; for (int i = 0; i < 5; i++) imb[i] += (v >> i & 1) ? -1 : 1;
    rec(depth + 1, v + 1);
    for (int i = 0; i < 5; i++) imb[i] -= (v >> i & 1) ? -1 : 1;
  }
}
int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: searchA k a repfile [maxdef]\n"); return 1; }
  K = atoi(argv[1]); A_ = atoi(argv[2]); B_ = K - A_; DELTA = A_ - B_;
  if (DELTA < 0 || K > 15 || K < 1) { fprintf(stderr, "need 1<=k<=15, a>=k-a\n"); return 1; }
  if (argc > 4) MAXDEF = atoi(argv[4]);
  FILE *f = fopen(argv[3], "r"); if (!f) { fprintf(stderr, "no file\n"); return 1; }
  int e = 0;
  for (int u = 0; u < 64; u++) for (int i = 0; i < 6; i++) if (!(u >> i & 1)) {
    int v = u | 1 << i;
    for (int w = 0; w < 64; w++) { int a = __builtin_popcount(u ^ w), b = __builtin_popcount(v ^ w); inc[w][e] = 1u << (4 * (a < b ? a : b)); }
    e++;
  }
  char line[64]; long nrep = 0;
  while (fgets(line, sizeof line, f)) {
    uint32_t m = (uint32_t)strtoul(line, 0, 16); nrep++;
    if (__builtin_popcount(m) != A_) { fprintf(stderr, "bad rep size\n"); return 1; }
    curA = m; memset(lvl[0], 0, sizeof lvl[0]); for (int i = 0; i < 5; i++) imb[i] = 0;
    for (int v = 0; v < 32; v++) if (m >> v & 1) {
      for (int e2 = 0; e2 < NE; e2++) lvl[0][e2] += inc[v][e2];
      for (int i = 0; i < 5; i++) imb[i] += (v >> i & 1) ? -1 : 1;
    }
    rec(0, 0);
  }
  printf("SUMMARY-A k=%d a=%d b=%d reps=%ld leaves=%lld pruned=%lld found=%lld near=%lld\n", K, A_, B_, nrep, leaves, pruned, found, near);
  return 0;
}
