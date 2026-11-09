/* sc_odd11.c - mac-mini-2026-11-08-firstrealizers. All self-converse tournaments of order 11 with H < LIM (default 625).
 * A self-converse tournament of odd order has an anti-involution sigma fixing exactly one vertex v0, and all of these
 * are conjugate, so WLOG sigma = (0 5)(1 6)(2 7)(3 8)(4 9) with v0 = 10. Then T - v0 is anti-invariant under the
 * fixed-point-free involution on {0..9} (enumerated as in sc_even10.c, 2^25 tournaments), and H(T - v0) <= H(T)
 * (insertion monotonicity, THM-1370). So every such T with H(T) < LIM extends some order-10 anti-invariant T' with
 * H(T') < LIM by v0, where v0 -> a_i iff b_i -> v0 (5 bits). Prints "E h strong" for every H < LIM attained at order 11. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned long long u64;
#define LIMDEF 625
static u64 f[1 << 11][11];
static u64 hp(const uint16_t *adj, int n){
  int F = 1 << n;
  for (int S = 0; S < F; S++) for (int v = 0; v < n; v++) f[S][v] = 0;
  for (int v = 0; v < n; v++) f[1 << v][v] = 1;
  for (int S = 1; S < F; S++) for (int v = 0; v < n; v++) if (f[S][v]) { u64 x = f[S][v]; uint16_t out = adj[v] & ~S;
    while (out) { int u = __builtin_ctz(out); out &= out - 1; f[S | 1 << u][u] += x; } }
  u64 h = 0; for (int v = 0; v < n; v++) h += f[F - 1][v]; return h;
}
static int strongq(const uint16_t *adj, int n){
  uint16_t all = (1u << n) - 1, sn = 1, fr = 1;
  while (fr) { uint16_t nx = 0; for (int v = 0; v < n; v++) if (fr >> v & 1) nx |= adj[v]; nx &= ~sn; sn |= nx; fr = nx; }
  if (sn != all) return 0;
  sn = 1; fr = 1;
  while (fr) { uint16_t nx = 0; for (int v = 0; v < n; v++) if (!(sn >> v & 1)) for (int u = 0; u < n; u++) if ((fr >> u & 1) && (adj[v] >> u & 1)) { nx |= 1u << v; break; } sn |= nx; fr = nx; }
  return sn == all;
}
int main(int argc, char **argv){
  int LIM = argc > 1 ? atoi(argv[1]) : LIMDEF;
  static unsigned char seen11[1 << 16], strong11[1 << 16];
  int pairs[10][2], np = 0; for (int i = 0; i < 5; i++) for (int j = i + 1; j < 5; j++) { pairs[np][0] = i; pairs[np][1] = j; np++; }
  long small10 = 0;
  for (u64 code = 0; code < (1ULL << 25); code++) {
    uint16_t adj[11] = {0};
    for (int i = 0; i < 5; i++) { if (code >> i & 1) adj[i] |= 1u << (i + 5); else adj[i + 5] |= 1u << i; }
    for (int p = 0; p < np; p++) {
      int i = pairs[p][0], j = pairs[p][1]; int x = code >> (5 + 2 * p) & 1, y = code >> (6 + 2 * p) & 1;
      int ai = i, bi = i + 5, aj = j, bj = j + 5;
      if (x) { adj[ai] |= 1u << aj; adj[bj] |= 1u << bi; } else { adj[aj] |= 1u << ai; adj[bi] |= 1u << bj; }
      if (y) { adj[ai] |= 1u << bj; adj[aj] |= 1u << bi; } else { adj[bj] |= 1u << ai; adj[bi] |= 1u << aj; }
    }
    u64 h10 = hp(adj, 10);
    if (h10 >= (u64)LIM) continue;
    small10++;
    for (int e = 0; e < 32; e++) {
      uint16_t a2[11]; for (int v = 0; v < 10; v++) a2[v] = adj[v]; a2[10] = 0;
      for (int i = 0; i < 5; i++) {
        if (e >> i & 1) { a2[10] |= 1u << i; a2[i + 5] |= 1u << 10; }     /* v0 -> a_i and b_i -> v0 */
        else { a2[i] |= 1u << 10; a2[10] |= 1u << (i + 5); }              /* a_i -> v0 and v0 -> b_i */
      }
      u64 h = hp(a2, 11);
      if (h < (u64)LIM) { seen11[h] = 1; if (!strong11[h] && strongq(a2, 11)) strong11[h] = 1; }
    }
    if ((code & ((1ULL << 22) - 1)) == 0) fprintf(stderr, "%llu/33554432, order-10 parents with H < %d so far %ld\n", code, LIM, small10);
  }
  fprintf(stderr, "order-10 anti-invariant tournaments with H < %d: %ld\n", LIM, small10);
  for (int h = 1; h < LIM; h++) if (seen11[h]) printf("E %d %d\n", h, strong11[h]);
  return 0;
}
