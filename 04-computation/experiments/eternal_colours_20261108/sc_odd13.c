/* sc_odd13.c - mac-mini-2026-11-08-firstrealizers. All self-converse tournaments of orders 2M+2 and 2M+3 with H <= TMAX
 * (M = 5: orders 12 and 13).
 *
 * Model. Every self-converse tournament has an anti-involution (|Aut| odd), fixed-point-free at even order and with one
 * fixed vertex at odd order, and all of one type are conjugate. So take sigma = (a_0 b_0)...(a_M b_M)(v0) with
 * a_i = i, b_i = i + M for i < M (the base), a = a_M = 2M, b = b_M = 2M+1, v0 = 2M+2. The transpositions (a_p b_p)
 * commute with sigma, so conjugating by them (the reduction of audit K1) every within-pair arc is a_p -> b_p. Then
 * T - v0 (order 2M+2) and T - v0 - {a, b} (order 2M) are sigma-anti-invariant, and H does not increase on deleting a
 * vertex (insertion monotonicity, THM-1370). So, pruning at H <= TMAX at every level:
 *   level 0: the 2^(2 C(M,2)) reduced base codes (two bits x, y per base pair i < j, as in sc_even_pair.c);
 *   level 1: vertex a, 2 bits per base pair (x_i = [a_i -> a], y_i = [a -> b_i]), Gray code, H by splitting at a;
 *   level 2: vertex b = sigma(a), forced arcs (b -> b_i iff x_i, a_i -> b iff y_i) and a -> b, H by splitting at b;
 *   level 3: vertex v0, one bit per pair (z_p = [v0 -> a_p], then b_p -> v0; else a_p -> v0 -> b_p), Gray code, split.
 * Output: "E h strong" (order 2M+2) and "O h strong" (order 2M+3) for every H <= TMAX attained; "strong" = some
 * realisation is strong. With -DCHECK the split values are compared with a direct DP on a sample.
 * usage: sc_odd13 TMAX [C0 C1]   (reduced base codes C0 <= code < C1) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
#ifndef M
#define M 5
#endif
#define N0 (2 * M)
#define N1 (N0 + 1)
#define N2 (N0 + 2)
#define N3 (N0 + 3)
static u64 F0[1 << N0][N0], G0[1 << N0][N0], F1[1 << N1][N1], G1[1 << N1][N1], F2[1 << N2][N2], G2[1 << N2][N2];
static i64 P0[1 << N0], Q0[1 << N0], P1[1 << N1], Q1[1 << N1], P2[1 << N2], Q2[1 << N2];   /* one pair per level */
static void build(const uint16_t *adj, int n, u64 *F, u64 *G){      /* F[S*n+v]: paths covering S ending at v; G: starting */
  int Fm = 1 << n; uint16_t inn[16];
  for (int v = 0; v < n; v++) { inn[v] = 0; for (int u = 0; u < n; u++) if (adj[u] >> v & 1) inn[v] |= 1u << u; }
  memset(F, 0, sizeof(u64) * (size_t)Fm * n); memset(G, 0, sizeof(u64) * (size_t)Fm * n);
  for (int v = 0; v < n; v++) { F[(1 << v) * n + v] = 1; G[(1 << v) * n + v] = 1; }
  for (int S = 1; S < Fm; S++) { uint16_t Sv = (uint16_t)S;
    while (Sv) { int v = __builtin_ctz(Sv); Sv &= Sv - 1;
      u64 x = F[S * n + v]; if (x) { uint16_t o = adj[v] & ~S; while (o) { int u = __builtin_ctz(o); o &= o - 1; F[(S | 1 << u) * n + u] += x; } }
      u64 y = G[S * n + v]; if (y) { uint16_t p = inn[v] & ~S; while (p) { int u = __builtin_ctz(p); p &= p - 1; G[(S | 1 << u) * n + u] += y; } } } }
}
static u64 conv(int n, const i64 *P, const i64 *Q){                 /* H from P (ends before z) and Q (starts after z) */
  int Fm = 1 << n, all = Fm - 1; i64 h = P[all] + Q[all];
  for (int S = 1; S < all; S++) if (P[S]) h += P[S] * Q[all & ~S];
  return (u64)h;
}
static void initPQ(int n, const u64 *F, const u64 *G, uint16_t in, uint16_t out, i64 *P, i64 *Q){
  int Fm = 1 << n;
  for (int S = 1; S < Fm; S++) { i64 p = 0, q = 0; uint16_t a = in & S, b = out & S;
    while (a) { int v = __builtin_ctz(a); a &= a - 1; p += (i64)F[S * n + v]; }
    while (b) { int v = __builtin_ctz(b); b &= b - 1; q += (i64)G[S * n + v]; }
    P[S] = p; Q[S] = q; }
}
static void move(int n, const u64 *F, const u64 *G, int v, int to_in, i64 *P, i64 *Q){ /* v changes side: out -> in (to_in) or back */
  int Fm = 1 << n; int s = to_in ? 1 : -1;
  for (int S = 1 << v; S < Fm; S = (S + 1) | (1 << v)) { P[S] += s * (i64)F[S * n + v]; Q[S] -= s * (i64)G[S * n + v]; }
}
static u64 hdirect(const uint16_t *adj, int n){
  static u64 f[1 << N3][N3]; int Fm = 1 << n;
  for (int S = 0; S < Fm; S++) for (int v = 0; v < n; v++) f[S][v] = 0;
  for (int v = 0; v < n; v++) f[1 << v][v] = 1;
  for (int S = 1; S < Fm; S++) for (int v = 0; v < n; v++) if (f[S][v]) { u64 x = f[S][v]; uint16_t o = adj[v] & ~S; while (o) { int u = __builtin_ctz(o); o &= o - 1; f[S | 1 << u][u] += x; } }
  u64 h = 0; for (int v = 0; v < n; v++) h += f[Fm - 1][v]; return h;
}
static int strongq(const uint16_t *adj, int n){
  uint16_t all = (uint16_t)((1u << n) - 1), sn = 1, fr = 1;
  while (fr) { uint16_t nx = 0; for (int v = 0; v < n; v++) if (fr >> v & 1) nx |= adj[v]; nx &= ~sn; sn |= nx; fr = nx; }
  if (sn != all) return 0;
  sn = 1; fr = 1;
  while (fr) { uint16_t nx = 0; for (int v = 0; v < n; v++) if (!(sn >> v & 1) && (adj[v] & fr)) nx |= 1u << v; sn |= nx; fr = nx; }
  return sn == all;
}
static int anti_invariant(const uint16_t *adj, int n){              /* u -> v iff sigma(v) -> sigma(u) */
  int sg[16]; for (int p = 0; p <= M; p++) { int a = p < M ? p : N0, b = p < M ? p + M : N0 + 1; sg[a] = b; sg[b] = a; }
  if (n > N2) sg[N2] = N2;
  for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (u != v && (adj[u] >> v & 1) != (adj[sg[v]] >> sg[u] & 1)) return 0;
  return 1;
}
static void addv(uint16_t *adj, int n, uint16_t out){                /* new vertex n with out-set `out` over 0..n-1 */
  adj[n] = out; for (int v = 0; v < n; v++) if (!(out >> v & 1)) adj[v] |= 1u << n;
}
int main(int argc, char **argv){
  long TMAX = atol(argv[1]);
  int pr[64][2], np = 0; for (int i = 0; i < M; i++) for (int j = i + 1; j < M; j++) { pr[np][0] = i; pr[np][1] = j; np++; }
  u64 NB = 2ULL * np, C0 = argc > 2 ? strtoull(argv[2], NULL, 10) : 0, C1 = argc > 3 ? strtoull(argv[3], NULL, 10) : (1ULL << NB);
  unsigned char *se = calloc(TMAX + 1, 1), *ss = calloc(TMAX + 1, 1), *oe = calloc(TMAX + 1, 1), *os = calloc(TMAX + 1, 1);
  long n0 = 0, n1 = 0, n2 = 0, n3 = 0, nchk = 0, bad = 0;
  for (u64 code = C0; code < C1; code++) {
    uint16_t adj[N3 + 1]; memset(adj, 0, sizeof adj);
    for (int i = 0; i < M; i++) adj[i] |= 1u << (i + M);                       /* reduced: a_i -> b_i */
    for (int p = 0; p < np; p++) { int i = pr[p][0], j = pr[p][1], x = code >> (2 * p) & 1, y = code >> (2 * p + 1) & 1;
      int ai = i, bi = i + M, aj = j, bj = j + M;
      if (x) { adj[ai] |= 1u << aj; adj[bj] |= 1u << bi; } else { adj[aj] |= 1u << ai; adj[bi] |= 1u << bj; }
      if (y) { adj[ai] |= 1u << bj; adj[aj] |= 1u << bi; } else { adj[bj] |= 1u << ai; adj[bi] |= 1u << aj; } }
    build(adj, N0, &F0[0][0], &G0[0][0]);
    u64 h0 = 0; for (int v = 0; v < N0; v++) h0 += F0[(1 << N0) - 1][v];
    if (h0 > (u64)TMAX) continue;
    n0++;
    /* level 1: Gray code over the 2M bits of a; bit 2i = x_i (a_i -> a), bit 2i+1 = y_i (a -> b_i). Start: all 0, i.e.
       a -> a_i (a_i on the out side) and b_i -> a (b_i on the in side). */
    uint16_t in1 = 0, out1 = 0; for (int i = 0; i < M; i++) { out1 |= 1u << i; in1 |= 1u << (i + M); }
    initPQ(N0, &F0[0][0], &G0[0][0], in1, out1, P0, Q0);
    unsigned bits = 0;
    for (unsigned g = 0; g < (1u << (2 * M)); g++) {
      if (g) { int t = __builtin_ctz(g); bits ^= 1u << t; int i = t >> 1;
        if (!(t & 1)) move(N0, &F0[0][0], &G0[0][0], i, (bits >> t) & 1, P0, Q0);          /* x_i = 1: a_i -> a, a_i goes in */
        else move(N0, &F0[0][0], &G0[0][0], i + M, !((bits >> t) & 1), P0, Q0); }       /* y_i = 1: a -> b_i, b_i goes out */
      u64 h1 = conv(N0, P0, Q0);
      if (h1 > (u64)TMAX) continue;
      n1++;
      uint16_t a1[N3 + 1]; memcpy(a1, adj, sizeof a1);
      uint16_t aout = 0, bout = 0;
      for (int i = 0; i < M; i++) { int x = bits >> (2 * i) & 1, y = bits >> (2 * i + 1) & 1;
        if (!x) aout |= 1u << i; if (y) aout |= 1u << (i + M);                  /* a -> a_i iff !x_i ; a -> b_i iff y_i */
        if (x) bout |= 1u << (i + M); if (!y) bout |= 1u << i; }                /* b -> b_i iff x_i ; b -> a_i iff !y_i */
      addv(a1, N0, aout);
#ifdef CHECK
      if (nchk < 3000 && (n1 % 97) == 1) { nchk++; if (hdirect(a1, N1) != h1) bad++; }
#endif
      /* level 2: b with forced arcs and a -> b */
      build(a1, N1, &F1[0][0], &G1[0][0]);
      uint16_t bin = (uint16_t)(((1u << N1) - 1) & ~bout);                    /* includes a (a -> b) */
      initPQ(N1, &F1[0][0], &G1[0][0], bin, bout, P1, Q1);
      u64 h2 = conv(N1, P1, Q1);
      if (h2 > (u64)TMAX) continue;
      n2++;
      uint16_t a2[N3 + 1]; memcpy(a2, a1, sizeof a2); addv(a2, N1, bout);
#ifdef CHECK
      if ((n2 % 13) == 1 && nchk < 6000) { nchk++; if (hdirect(a2, N2) != h2 || !anti_invariant(a2, N2)) bad++; }
#endif
      se[h2] = 1; if (!ss[h2] && strongq(a2, N2)) ss[h2] = 1;
      /* level 3: v0, one bit per pair p = 0..M (a_M = N0, b_M = N0+1). Start z = 0: a_p -> v0 (in), v0 -> b_p (out). */
      build(a2, N2, &F2[0][0], &G2[0][0]);
      uint16_t in3 = 0, out3 = 0;
      for (int p = 0; p <= M; p++) { int ap = p < M ? p : N0, bp = p < M ? p + M : N0 + 1; in3 |= 1u << ap; out3 |= 1u << bp; }
      initPQ(N2, &F2[0][0], &G2[0][0], in3, out3, P2, Q2);
      unsigned z = 0;
      for (unsigned g = 0; g < (1u << (M + 1)); g++) {
        if (g) { int p = __builtin_ctz(g); z ^= 1u << p; int ap = p < M ? p : N0, bp = p < M ? p + M : N0 + 1;
          if (z >> p & 1) { move(N2, &F2[0][0], &G2[0][0], ap, 0, P2, Q2); move(N2, &F2[0][0], &G2[0][0], bp, 1, P2, Q2); }  /* v0 -> a_p, b_p -> v0 */
          else { move(N2, &F2[0][0], &G2[0][0], ap, 1, P2, Q2); move(N2, &F2[0][0], &G2[0][0], bp, 0, P2, Q2); } }
        u64 h3 = conv(N2, P2, Q2);
        if (h3 > (u64)TMAX) continue;
        n3++;
        oe[h3] = 1;
        int need = !os[h3];
#ifdef CHECK
        if (nchk < 9000 && (n3 % 7) == 1) need = 2;
#endif
        if (need) {
          uint16_t a3[N3 + 1]; memcpy(a3, a2, sizeof a3); uint16_t vout = 0;
          for (int p = 0; p <= M; p++) { int ap = p < M ? p : N0, bp = p < M ? p + M : N0 + 1; vout |= (z >> p & 1) ? (1u << ap) : (1u << bp); }
          addv(a3, N2, vout);
#ifdef CHECK
          if (need == 2) { nchk++; if (hdirect(a3, N3) != h3 || !anti_invariant(a3, N3)) bad++; }
#endif
          if (!os[h3] && strongq(a3, N3)) os[h3] = 1;
        }
      }
    }
    if ((code & ((1ULL << 16) - 1)) == 0) fprintf(stderr, "code %llu/%llu: bases %ld, order-%d %ld, order-%d %ld, order-%d %ld\n", code, 1ULL << NB, n0, N1, n1, N2, n2, N3, n3);
  }
  fprintf(stderr, "done [%llu,%llu) TMAX %ld: bases %ld; order %d kept %ld; order %d kept %ld; order %d found %ld; checks %ld, failures %ld\n",
          C0, C1, TMAX, n0, N1, n1, N2, n2, N3, n3, nchk, bad);
  for (long h = 1; h <= TMAX; h++) if (se[h]) printf("E %ld %d\n", h, ss[h]);
  for (long h = 1; h <= TMAX; h++) if (oe[h]) printf("O %ld %d\n", h, os[h]);
  return bad ? 3 : 0;
}
