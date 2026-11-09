/* sc_walk.c - mac-mini-2026-11-08-firstrealizers. Randomised witness search for STRONG self-converse tournaments of order
 * n (12..15) with small H. Every state is anti-invariant under the mirror mu(v) = n-1-v: start from the transitive TT_n
 * (i -> j for i < j), which mu reverses, and flip whole mirror orbits of arcs {i,j} <-> {n-1-j, n-1-i}; flipping an orbit
 * keeps mu an anti-automorphism, so every state is self-converse. A move is accepted iff the new H is <= HMAX (a random
 * walk on the anti-invariant tournaments with H <= HMAX). Every strong state with H in [HLO, HMAX] is recorded, with the
 * reversed orbits of the first witness. Only positive results (witnesses) are meaningful; nothing is exhaustive.
 * usage: sc_walk n HLO HMAX STEPS SEED > out   (output: "W h n 1 rev i-j ..." per new strong value; summary on stderr) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned long long u64;
#define MAXN 15
static int n;
static uint16_t adj[MAXN];
static u64 f[1 << MAXN][MAXN];
static u64 hpaths(void){
  int F = 1 << n; uint16_t inn[MAXN];
  for (int v = 0; v < n; v++) { inn[v] = 0; for (int u = 0; u < n; u++) if (adj[u] >> v & 1) inn[v] |= 1u << u; }
  for (int v = 0; v < n; v++) f[1 << v][v] = 1;
  for (int S = 1; S < F; S++) {
    if ((S & (S - 1)) == 0) continue;          /* singletons already set */
    uint16_t Sv = (uint16_t)S;
    while (Sv) { int v = __builtin_ctz(Sv); Sv &= Sv - 1;
      uint16_t S0 = (uint16_t)(S & ~(1u << v)), m = S0 & inn[v]; u64 x = 0;
      while (m) { int u = __builtin_ctz(m); m &= m - 1; x += f[S0][u]; }
      f[S][v] = x; } }
  u64 h = 0; for (int v = 0; v < n; v++) h += f[F - 1][v]; return h;
}
static int strongq(void){
  uint16_t all = (uint16_t)((1u << n) - 1), seen = 1, fr = 1;
  while (fr) { uint16_t nx = 0; for (int v = 0; v < n; v++) if (fr >> v & 1) nx |= adj[v]; nx &= ~seen; seen |= nx; fr = nx; }
  if (seen != all) return 0;
  seen = 1; fr = 1;
  while (fr) { uint16_t nx = 0; for (int v = 0; v < n; v++) if (!(seen >> v & 1)) for (int u = 0; u < n; u++) if ((fr >> u & 1) && (adj[v] >> u & 1)) { nx |= 1u << v; break; } seen |= nx; fr = nx; }
  return seen == all;
}
static int oi[128], oj[128], no;
static void flip(int r){                         /* reverse orbit r: arc {i,j} and its mirror */
  int i = oi[r], j = oj[r], mi = n - 1 - j, mj = n - 1 - i;
  adj[i] ^= 1u << j; adj[j] ^= 1u << i;
  if (!(mi == i && mj == j)) { adj[mi] ^= 1u << mj; adj[mj] ^= 1u << mi; }
}
static u64 rs;
static u64 rnd(void){ rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
int main(int argc, char **argv){
  n = atoi(argv[1]); long HLO = atol(argv[2]), HMAX = atol(argv[3]); long STEPS = atol(argv[4]); rs = strtoull(argv[5], NULL, 10) * 2654435761ULL + 88172645463325252ULL;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
    int mi = n - 1 - j, mj = n - 1 - i;
    if (mi < i || (mi == i && mj < j)) continue;
    oi[no] = i; oj[no] = j; no++; }
  for (int v = 0; v < n; v++) { adj[v] = 0; for (int u = v + 1; u < n; u++) adj[v] |= 1u << u; }
  unsigned char *st = calloc(HMAX + 1, 1); char (*wit)[400] = calloc(HMAX + 1, 400);
  unsigned char *rev = calloc(no, 1);
  u64 h = hpaths(); long acc = 0, newv = 0;
  for (long s = 0; s < STEPS; s++) {
    int r = (int)(rnd() % (u64)no); flip(r);
    u64 h2 = hpaths();
    if (h2 > (u64)HMAX) { flip(r); continue; }
    rev[r] ^= 1; h = h2; acc++;
    if ((long)h >= HLO && !st[h] && strongq()) {
      st[h] = 1; newv++; int p = snprintf(wit[h], 400, "rev");
      for (int t = 0; t < no && p < 390; t++) if (rev[t]) p += snprintf(wit[h] + p, 400 - p, " %d-%d", oi[t], oj[t]);
    }
    if ((s & 0xFFFFF) == 0) fprintf(stderr, "n=%d step %ld accepted %ld strong values %ld\n", n, s, acc, newv);
  }
  for (long v = HLO; v <= HMAX; v++) if (st[v]) printf("W %ld %d 1 %s\n", v, n, wit[v]);
  fprintf(stderr, "done n=%d HLO=%ld HMAX=%ld steps=%ld accepted=%ld strong values=%ld (orbits %d)\n", n, HLO, HMAX, STEPS, acc, newv, no);
  return 0;
}
