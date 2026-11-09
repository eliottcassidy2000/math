/* aud_twins_search.c - audit upgrade for HYP-9530 (C), second half (opus-2026-11-08-S86 audit).
 * For each target value h = 1 + 2a (a a carry-free singleton read from stdin, one per line), search for a STRONG
 * tournament with H = h that is NOT intercyclic (it has two vertex-disjoint cyclic triangles), i.e. a "carried twin"
 * of the unique strong intercyclic realiser.  Search: random restarts + greedy/annealed single-arc flips on order n,
 * H by subset DP.  Every reported witness is re-verified (H by DP, strong, two disjoint triangles) and printed as
 * out-neighbour masks.
 * usage: aud_twins_search n tries_per_target < values.txt
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
static int N;
static u64 rng = 0x9E3779B97F4A7C15ULL;
static u64 rnd(void){ rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static u64 hp(const uint32_t *adj){
  static u64 f[1 << 12][12]; int F = 1 << N;
  for (int S = 0; S < F; S++) for (int v = 0; v < N; v++) f[S][v] = 0;
  for (int v = 0; v < N; v++) f[1 << v][v] = 1;
  for (int S = 1; S < F; S++) for (int v = 0; v < N; v++){ u64 x = f[S][v]; if (!x) continue; uint32_t o = adj[v] & ~(uint32_t)S; while (o){ int u = __builtin_ctz(o); o &= o - 1; f[S | (1 << u)][u] += x; } }
  u64 h = 0; for (int v = 0; v < N; v++) h += f[F - 1][v]; return h;
}
static int strongq(const uint32_t *adj){
  uint32_t full = (1u << N) - 1, R = 1, fr = 1;
  while (fr){ uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= adj[__builtin_ctz(f)]; nx &= full & ~R; R |= nx; fr = nx; }
  if (R != full) return 0;
  uint32_t in[32] = {0}; for (int v = 0; v < N; v++) for (int u = 0; u < N; u++) if (adj[v] >> u & 1) in[u] |= 1u << v;
  R = 1; fr = 1;
  while (fr){ uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= in[__builtin_ctz(f)]; nx &= full & ~R; R |= nx; fr = nx; }
  return R == full;
}
static int disjoint_triangles(const uint32_t *adj){
  uint32_t tri[2000]; int t = 0;
  for (int a = 0; a < N; a++) for (int b = a + 1; b < N; b++) for (int c = b + 1; c < N; c++){
    int ab = adj[a] >> b & 1, bc = adj[b] >> c & 1, ca = adj[c] >> a & 1;
    if (ab == bc && bc == ca) tri[t++] = (1u << a) | (1u << b) | (1u << c); }
  for (int i = 0; i < t; i++) for (int j = i + 1; j < t; j++) if (!(tri[i] & tri[j])) return 1;
  return 0;
}
static void flip(uint32_t *adj, int i, int j){ adj[i] ^= 1u << j; adj[j] ^= 1u << i; }
int main(int argc, char **argv){
  N = argc > 1 ? atoi(argv[1]) : 10; long tries = argc > 2 ? atol(argv[2]) : 200;
  long a;
  while (scanf("%ld", &a) == 1){
    u64 target = 1 + 2 * (u64)a; int found = 0;
    for (long t = 0; t < tries && !found; t++){
      uint32_t adj[32] = {0};
      for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++){ if (rnd() & 1) adj[i] |= 1u << j; else adj[j] |= 1u << i; }
      u64 h = hp(adj);
      for (int step = 0; step < 4000; step++){
        if (h == target && strongq(adj) && disjoint_triangles(adj)){ found = 1; break; }
        int i = (int)(rnd() % N), j = (int)(rnd() % N); if (i == j) continue;
        flip(adj, i, j); u64 h2 = hp(adj);
        long long d1 = (long long)h - (long long)target, d2 = (long long)h2 - (long long)target;
        if (d1 < 0) d1 = -d1; if (d2 < 0) d2 = -d2;
        if (d2 <= d1 || (rnd() % 100) < 3) h = h2; else flip(adj, i, j);
      }
      if (found){
        /* re-verify from scratch */
        int ok = hp(adj) == target && strongq(adj) && disjoint_triangles(adj);
        printf("a = %ld H = %llu order %d strong non-intercyclic witness %s, out-masks:", a, target, N, ok ? "VERIFIED" : "FAILED");
        for (int v = 0; v < N; v++) printf(" %u", adj[v]); printf("\n"); fflush(stdout);
      }
    }
    if (!found){ printf("a = %ld H = %llu order %d: no witness found in %ld restarts\n", a, target, N, tries); fflush(stdout); }
  }
  return 0;
}
