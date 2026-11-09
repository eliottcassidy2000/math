/* sc_even_pair.c - mac-mini-2026-11-08-firstrealizers. All self-converse tournaments of even order 2M+2 with H <= TMAX.
 * Every self-converse tournament of even order has a fixed-point-free anti-involution sigma, all such are conjugate, so
 * WLOG sigma = (0 M)(1 M+1)...(M-1 2M-1)(a b) with the new pair a = 2M, b = 2M+1. Removing the pair leaves an
 * anti-invariant T' of order 2M, and H(T') <= H(T' + a) <= H(T) (insertion monotonicity, THM-1370). So we enumerate all
 * 2^(M + 2 C(M,2)) anti-invariant T' with H(T') <= TMAX, then a's arcs (bits x_i = [a_i -> a], y_i = [a -> b_i] for each
 * old pair i; anti-invariance then gives b -> b_i iff x_i and a_i -> b iff y_i), keeping H(T' + a) <= TMAX, then the
 * a-b arc. H is evaluated by splitting paths at the newest vertex (forward/backward subset tables).
 * Output: "P h strong" for every H <= TMAX of a self-converse tournament of order 2M+2.
 * usage: sc_even_pair TMAX [C0 C1]   (compile with -DM=3 for order 8 validation, -DM=5 for order 12) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned long long u64;
#ifndef M
#define M 5
#endif
#define N0 (2*M)          /* order of T' */
static u64 F0[1 << N0][N0], G0[1 << N0][N0];
static u64 F1[1 << (N0 + 1)][N0 + 1], G1[1 << (N0 + 1)][N0 + 1];
static void tables(const uint16_t *adj, int n, u64 (*F)[N0 + 1], u64 (*G)[N0 + 1], int stride_dummy){ (void)stride_dummy;
  int Fm = 1 << n;
  uint16_t inn[N0 + 2];
  for (int v = 0; v < n; v++) { inn[v] = 0; for (int u = 0; u < n; u++) if (adj[u] >> v & 1) inn[v] |= 1u << u; }
  for (int S = 0; S < Fm; S++) for (int v = 0; v < n; v++) { F[S][v] = 0; G[S][v] = 0; }
  for (int v = 0; v < n; v++) { F[1 << v][v] = 1; G[1 << v][v] = 1; }
  for (int S = 1; S < Fm; S++) for (int v = 0; v < n; v++) {
    if (F[S][v]) { u64 x = F[S][v]; uint16_t out = adj[v] & ~S; while (out) { int u = __builtin_ctz(out); out &= out - 1; F[S | 1 << u][u] += x; } }
    if (G[S][v]) { u64 x = G[S][v]; uint16_t pre = inn[v] & ~S; while (pre) { int u = __builtin_ctz(pre); pre &= pre - 1; G[S | 1 << u][u] += x; } } }
}
static void tables0(const uint16_t *adj){                     /* order N0 tables in F0/G0 */
  int n = N0, Fm = 1 << n; uint16_t inn[N0];
  for (int v = 0; v < n; v++) { inn[v] = 0; for (int u = 0; u < n; u++) if (adj[u] >> v & 1) inn[v] |= 1u << u; }
  for (int S = 0; S < Fm; S++) for (int v = 0; v < n; v++) { F0[S][v] = 0; G0[S][v] = 0; }
  for (int v = 0; v < n; v++) { F0[1 << v][v] = 1; G0[1 << v][v] = 1; }
  for (int S = 1; S < Fm; S++) for (int v = 0; v < n; v++)
    if (F0[S][v]) { u64 x = F0[S][v]; uint16_t out = adj[v] & ~S; while (out) { int u = __builtin_ctz(out); out &= out - 1; F0[S | 1 << u][u] += x; } }
}
static void tables0_back(const uint16_t *adj){                /* backward table, only for kept bases */
  int n = N0, Fm = 1 << n; uint16_t inn[N0];
  for (int v = 0; v < n; v++) { inn[v] = 0; for (int u = 0; u < n; u++) if (adj[u] >> v & 1) inn[v] |= 1u << u; }
  for (int S = 1; S < Fm; S++) for (int v = 0; v < n; v++)
    if (G0[S][v]) { u64 x = G0[S][v]; uint16_t pre = inn[v] & ~S; while (pre) { int u = __builtin_ctz(pre); pre &= pre - 1; G0[S | 1 << u][u] += x; } }
}
/* H of (base + new vertex z) from base tables: z's out-set Nout and in-set Nin over base vertices */
static u64 split0(uint16_t Nout){
  int Fm = 1 << N0; static u64 P[1 << N0], Q[1 << N0];
  for (int S = 1; S < Fm; S++) { u64 p = 0, q = 0;
    for (int v = 0; v < N0; v++) if (S >> v & 1) { if (Nout >> v & 1) q += G0[S][v]; else p += F0[S][v]; }
    P[S] = p; Q[S] = q; }
  u64 h = P[Fm - 1] + Q[Fm - 1];
  for (int S = 1; S < Fm - 1; S++) h += P[S] * Q[(Fm - 1) & ~S];
  return h;
}
static u64 split1(uint16_t Nout){
  int n = N0 + 1, Fm = 1 << n; static u64 P[1 << (N0 + 1)], Q[1 << (N0 + 1)];
  for (int S = 1; S < Fm; S++) { u64 p = 0, q = 0;
    for (int v = 0; v < n; v++) if (S >> v & 1) { if (Nout >> v & 1) q += G1[S][v]; else p += F1[S][v]; }
    P[S] = p; Q[S] = q; }
  u64 h = P[Fm - 1] + Q[Fm - 1];
  for (int S = 1; S < Fm - 1; S++) h += P[S] * Q[(Fm - 1) & ~S];
  return h;
}
static u64 hp_direct(const uint16_t *adj, int n){ /* direct DP for strongness-independent checks */
  static u64 f[1 << (N0 + 2)][N0 + 2]; int Fm = 1 << n;
  for (int S = 0; S < Fm; S++) for (int v = 0; v < n; v++) f[S][v] = 0;
  for (int v = 0; v < n; v++) f[1 << v][v] = 1;
  for (int S = 1; S < Fm; S++) for (int v = 0; v < n; v++) if (f[S][v]) { u64 x = f[S][v]; uint16_t out = adj[v] & ~S; while (out) { int u = __builtin_ctz(out); out &= out - 1; f[S | 1 << u][u] += x; } }
  u64 h = 0; for (int v = 0; v < n; v++) h += f[Fm - 1][v]; return h;
}
static int strongq(const uint16_t *adj, int n){
  uint16_t all = (uint16_t)((1u << n) - 1), sn = 1, fr = 1;
  while (fr) { uint16_t nx = 0; for (int v = 0; v < n; v++) if (fr >> v & 1) nx |= adj[v]; nx &= ~sn; sn |= nx; fr = nx; }
  if (sn != all) return 0;
  sn = 1; fr = 1;
  while (fr) { uint16_t nx = 0; for (int v = 0; v < n; v++) if (!(sn >> v & 1)) for (int u = 0; u < n; u++) if ((fr >> u & 1) && (adj[v] >> u & 1)) { nx |= 1u << v; break; } sn |= nx; fr = nx; }
  return sn == all;
}
int main(int argc, char **argv){
  long TMAX = atol(argv[1]);
  int pairs[64][2], np = 0; for (int i = 0; i < M; i++) for (int j = i + 1; j < M; j++) { pairs[np][0] = i; pairs[np][1] = j; np++; }
  u64 NB = (u64)M + 2ULL * np;
  u64 C0 = argc > 2 ? strtoull(argv[2], NULL, 10) : 0, C1 = argc > 3 ? strtoull(argv[3], NULL, 10) : (1ULL << NB);
  unsigned char *seen = calloc(TMAX + 1, 1), *seens = calloc(TMAX + 1, 1);
  long nbase = 0, nmid = 0, check_fail = 0, nchecked = 0;
  for (u64 code = C0; code < C1; code++) {
    uint16_t adj[N0 + 2] = {0};
    for (int i = 0; i < M; i++) { if (code >> i & 1) adj[i] |= 1u << (i + M); else adj[i + M] |= 1u << i; }
    for (int p = 0; p < np; p++) {
      int i = pairs[p][0], j = pairs[p][1]; int x = code >> (M + 2 * p) & 1, y = code >> (M + 1 + 2 * p) & 1;
      int ai = i, bi = i + M, aj = j, bj = j + M;
      if (x) { adj[ai] |= 1u << aj; adj[bj] |= 1u << bi; } else { adj[aj] |= 1u << ai; adj[bi] |= 1u << bj; }
      if (y) { adj[ai] |= 1u << bj; adj[aj] |= 1u << bi; } else { adj[bj] |= 1u << ai; adj[bi] |= 1u << aj; }
    }
    tables0(adj);
    u64 h0 = 0; for (int v = 0; v < N0; v++) h0 += F0[(1 << N0) - 1][v];
    if (h0 > (u64)TMAX) continue;
    tables0_back(adj);
    nbase++;
    const int A = N0, B = N0 + 1;
    for (int bits = 0; bits < (1 << (2 * M)); bits++) {
      /* a's arcs: a_i -> a iff x_i ; a -> b_i iff y_i */
      uint16_t aout = 0;
      for (int i = 0; i < M; i++) { int x = bits >> (2 * i) & 1, y = bits >> (2 * i + 1) & 1;
        if (!x) aout |= 1u << i;          /* a -> a_i */
        if (y) aout |= 1u << (i + M); }   /* a -> b_i */
      u64 ha = split0(aout);
      if (ha > (u64)TMAX) continue;
      nmid++;
      uint16_t a1[N0 + 2]; for (int v = 0; v < N0; v++) a1[v] = adj[v];
      a1[A] = aout; for (int v = 0; v < N0; v++) if (!(aout >> v & 1)) a1[v] |= 1u << A;
      tables(a1, N0 + 1, F1, G1, 0);
      /* b's arcs over the N0 old vertices: b -> b_i iff x_i ; a_i -> b iff y_i */
      uint16_t bout0 = 0;
      for (int i = 0; i < M; i++) { int x = bits >> (2 * i) & 1, y = bits >> (2 * i + 1) & 1;
        if (x) bout0 |= 1u << (i + M);    /* b -> b_i */
        if (!y) bout0 |= 1u << i; }       /* b -> a_i */
      for (int ab = 0; ab < 2; ab++) {   /* ab = 1: a -> b ; 0: b -> a */
        uint16_t bout = bout0 | (ab ? 0 : (1u << A));
        u64 h = split1(bout);
        if (h <= (u64)TMAX) {
          uint16_t t[N0 + 2]; for (int v = 0; v <= A; v++) t[v] = a1[v];
          t[B] = bout; for (int v = 0; v <= A; v++) if (!(bout >> v & 1)) t[v] |= 1u << B;
          if (nchecked < 2000) { nchecked++; if (hp_direct(t, N0 + 2) != h) check_fail++; }
          seen[h] = 1; if (!seens[h] && strongq(t, N0 + 2)) seens[h] = 1;
        }
      }
    }
    if ((code & ((1ULL << 20) - 1)) == 0) fprintf(stderr, "%llu/%llu base<=TMAX %ld mid %ld\n", code, 1ULL << NB, nbase, nmid);
  }
  fprintf(stderr, "done: base tournaments with H <= %ld: %ld; one-vertex extensions kept: %ld; direct-DP checks %ld, failures %ld\n", TMAX, nbase, nmid, nchecked, check_fail);
  for (long h = 1; h <= TMAX; h++) if (seen[h]) printf("P %ld %d\n", h, seens[h]);
  return check_fail ? 3 : 0;
}
