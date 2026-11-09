/* sc_even10.c - mac-mini-2026-11-08-firstrealizers. All self-converse tournaments of even order N (default 10; -DN=6, 8 for checks), labelled and anti-invariant
 * under the fixed-point-free involution s0 = (0 5)(1 6)(2 7)(3 8)(4 9). Every self-converse tournament has an
 * anti-involution (Aut has odd order, so a^(odd) is an involution), at even order it is fixed-point-free, and all of these
 * are conjugate, so every self-converse class of order 10 occurs. 2^25 tournaments: 5 within-pair bits and, for each pair
 * of pairs i<j, x = [a_i->a_j] (then b_j->b_i iff x) and y = [a_i->b_j] (then a_j->b_i iff y). We record the H values
 * of self-converse tournaments and of strong ones (H < HM). Output: lines "S h strong_flag" for every attained h. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned long long u64;
#ifndef N
#define N 10
#endif
#define M (N/2)
#define HM 20000
static u64 f[1 << N][N];
int main(void){
  static unsigned char seen[HM], seens[HM];
  int pairs[64][2], np = 0; for (int i = 0; i < M; i++) for (int j = i + 1; j < M; j++) { pairs[np][0] = i; pairs[np][1] = j; np++; }
  u64 NB = (u64)M + 2ULL * (u64)np;
  for (u64 code = 0; code < (1ULL << NB); code++) {
    uint16_t adj[N] = {0};
    for (int i = 0; i < M; i++) { if (code >> i & 1) adj[i] |= 1u << (i + M); else adj[i + M] |= 1u << i; }
    for (int p = 0; p < np; p++) {
      int i = pairs[p][0], j = pairs[p][1]; int x = code >> (M + 2 * p) & 1, y = code >> (M + 1 + 2 * p) & 1;
      int ai = i, bi = i + M, aj = j, bj = j + M;
      if (x) { adj[ai] |= 1u << aj; adj[bj] |= 1u << bi; } else { adj[aj] |= 1u << ai; adj[bi] |= 1u << bj; }
      if (y) { adj[ai] |= 1u << bj; adj[aj] |= 1u << bi; } else { adj[bj] |= 1u << ai; adj[bi] |= 1u << aj; }
    }
    int F = 1 << N;
    for (int S = 0; S < F; S++) for (int v = 0; v < N; v++) f[S][v] = 0;
    for (int v = 0; v < N; v++) f[1 << v][v] = 1;
    for (int S = 1; S < F; S++) for (int v = 0; v < N; v++) if (f[S][v]) { u64 xx = f[S][v]; uint16_t out = adj[v] & ~S;
      while (out) { int u = __builtin_ctz(out); out &= out - 1; f[S | 1 << u][u] += xx; } }
    u64 h = 0; for (int v = 0; v < N; v++) h += f[F - 1][v];
    if (h >= HM) { fprintf(stderr, "H too large\n"); return 1; }
    if (!seen[h] || !seens[h]) {
      /* strong check */
      uint16_t all = (1u << N) - 1, sn = 1, fr = 1;
      while (fr) { uint16_t nx = 0; for (int v = 0; v < N; v++) if (fr >> v & 1) nx |= adj[v]; nx &= ~sn; sn |= nx; fr = nx; }
      int strong = 0;
      if (sn == all) { sn = 1; fr = 1;
        while (fr) { uint16_t nx = 0; for (int v = 0; v < N; v++) if (!(sn >> v & 1)) for (int u = 0; u < N; u++) if ((fr >> u & 1) && (adj[v] >> u & 1)) { nx |= 1u << v; break; } sn |= nx; fr = nx; }
        strong = (sn == all); }
      seen[h] = 1; if (strong) seens[h] = 1;
    }
    if ((code & ((1ULL << 22) - 1)) == 0) fprintf(stderr, "%llu/%llu\n", code, 1ULL << NB);
  }
  for (int h = 1; h < HM; h++) if (seen[h]) printf("S %d %d\n", h, seens[h]);
  return 0;
}
