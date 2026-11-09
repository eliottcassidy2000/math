/* aud_scpair.c - audit K1 (2026-11-08). Validation copy of aud_sc12.c with the order as a parameter: usage aud_scpair n LIM C0 C1
 * (n even, 4 <= n <= 12). With LIM above max H(n) it must reproduce aud_sc.c's reduced order-n counts exactly (bijection).
 * Original header of aud_sc12.c follows.
 * All self-converse tournaments of order 12 with H < LIM (default 625), to decide h = 377 and to test THM-791 there.
 *
 * Model: sigma(v) = v ^ 1 on {0..11}. Every self-converse 12-tournament has a fixed-point-free anti-involution
 * (|Aut| odd), all such involutions are conjugate, so it is isomorphic to a sigma-anti-invariant T; conjugating by the
 * transpositions (2k 2k+1) (which commute with sigma) we may take every within-pair arc as 2k -> 2k+1.
 * Then T'' = T - {10, 11} is sigma-anti-invariant on {0..9} with forward within-pair arcs, i.e. one of the 2^20
 * "reduced" order-10 codes of aud_sc.c, and by insertion monotonicity (H(T - v) <= H(T): insert v at its first feasible
 * position, an injection) H(T'') <= H(T - 11) <= H(T) < LIM.
 * So: for each reduced order-10 code with H < LIM, for each of the 2^10 out-masks N of vertex 10 on {0..9}, compute
 * H(T'' + 10) by splitting each path at vertex 10 (forward/backward subset tables of T''); keep it if < LIM; then add
 * vertex 11 = sigma(10) with the forced arcs (w -> 11 iff 10 -> sigma(w)) and 10 -> 11, and compute H(T) by a direct
 * pull DP. Prints "h count strongcount" for H < LIM. With CHECK the split value is compared with a direct DP for a
 * sample of extensions, and anti-invariance of every order-12 tournament is verified.
 * usage: aud_sc12 LIM C0 C1   (C0 <= reduced order-10 code < C1, codes in [0, 2^20))
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef uint32_t u32; typedef uint64_t u64;
#define MAXN 12
static int n, sig[MAXN];
static int oi[128], oj[128], nother, pi_[16], pj_[16], npair;
static uint16_t outm[MAXN], inm[MAXN];
static u32 *F;
static void setarc(int a, int b){ outm[a] |= (uint16_t)(1u << b); inm[b] |= (uint16_t)(1u << a); }
static void orbits(int nn){
  n = nn; int m = n / 2; for (int v = 0; v < n; v++) sig[v] = (v < 2 * m) ? (v ^ 1) : v;
  static int vis[MAXN][MAXN]; memset(vis, 0, sizeof vis); nother = npair = 0;
  for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++){
    if (vis[u][v]) continue;
    int su = sig[u], sv = sig[v]; int a = su < sv ? su : sv, b = su < sv ? sv : su;
    vis[u][v] = 1; vis[a][b] = 1;
    if (a == u && b == v){ pi_[npair] = u; pj_[npair] = v; npair++; } else { oi[nother] = u; oj[nother] = v; nother++; }
  }
}
static void build10(u64 c){ /* identical construction to aud_sc.c (reduced mode), n = 10 */
  u64 code = c | (((1ULL << npair) - 1) << nother);
  memset(outm, 0, sizeof outm); memset(inm, 0, sizeof inm);
  for (int k = 0; k < nother; k++){ int a = oi[k], b = oj[k]; if (!((code >> k) & 1)){ int t = a; a = b; b = t; } setarc(a, b); setarc(sig[b], sig[a]); }
  for (int k = 0; k < npair; k++){ int a = pi_[k], b = pj_[k]; if (!((code >> (nother + k)) & 1)){ int t = a; a = b; b = t; } setarc(a, b); }
}
static u32 hpaths(int nn){
  int FN = 1 << nn;
  for (int v = 0; v < nn; v++) F[(1u << v) * nn + v] = 1;
  for (int S = 3; S < FN; S++){ if (!(S & (S - 1))) continue;
    for (int v = 0; v < nn; v++) if (S >> v & 1){ int R = S & ~(1 << v); uint16_t m = inm[v] & (uint16_t)R; u32 s = 0;
      while (m){ int u = __builtin_ctz(m); m &= m - 1; s += F[(u32)R * nn + u]; } F[(u32)S * nn + v] = s; } }
  u32 h = 0; for (int v = 0; v < nn; v++) h += F[(u32)(FN - 1) * nn + v]; return h;
}
static int landau_strong(int nn){
  int sc[MAXN]; for (int v = 0; v < nn; v++) sc[v] = __builtin_popcount(outm[v]);
  for (int i = 1; i < nn; i++){ int x = sc[i], j = i - 1; while (j >= 0 && sc[j] > x){ sc[j + 1] = sc[j]; j--; } sc[j + 1] = x; }
  int s = 0; for (int k = 1; k < nn; k++){ s += sc[k - 1]; if (s == k * (k - 1) / 2) return 0; }
  return 1;
}
static u32 fw[1024][10], bw[1024][10]; static int NN, PP, PF;
int main(int argc, char **argv){
  NN = atoi(argv[1]); PP = NN - 2; PF = (1 << PP) - 1;
  int LIM = atoi(argv[2]);
  orbits(PP);
  u64 C0 = argc > 3 ? strtoull(argv[3], 0, 10) : 0, C1 = argc > 4 ? strtoull(argv[4], 0, 10) : (1ULL << nother);
  F = calloc((size_t)(1u << 12) * 12, sizeof(u32));
  static u64 cnt[1 << 18], scnt[1 << 18];
  long parents = 0, kept11 = 0, nchk = 0;
  int sig12[12]; for (int v = 0; v < NN; v++) sig12[v] = v ^ 1;
  for (u64 c = C0; c < C1; c++){
    orbits(PP); build10(c);
    u32 h10 = hpaths(PP);
    if (h10 >= (u32)LIM) continue;
    parents++;
    uint16_t o10[12], i10[12]; memcpy(o10, outm, sizeof o10); memcpy(i10, inm, sizeof i10);
    /* forward table fw[S][v] (paths covering S ending at v) and backward table bw[S][v] (paths covering S starting at v) */
    for (int S = 1; S <= PF; S++) for (int v = 0; v < PP; v++){
      if (!(S >> v & 1)){ fw[S][v] = bw[S][v] = 0; continue; }
      if (S == (1 << v)){ fw[S][v] = bw[S][v] = 1; continue; }
      int R = S & ~(1 << v); u32 a = 0, b = 0;
      uint16_t mi = i10[v] & R, mo = o10[v] & R;
      while (mi){ int u = __builtin_ctz(mi); mi &= mi - 1; a += fw[R][u]; }
      while (mo){ int u = __builtin_ctz(mo); mo &= mo - 1; b += bw[R][u]; }
      fw[S][v] = a; bw[S][v] = b;
    }
    /* N = out-neighbours of vertex 10 in {0..9}, visited in Gray-code order; P[S] = paths over S ending at a vertex
       v -> 10, Q[S] = paths over S starting at a vertex w with 10 -> w, updated incrementally when one bit of N flips */
    static u64 P[1024], Q[1024];
    P[0] = Q[0] = 0;
    for (int S = 1; S <= PF; S++){ u64 p = 0; for (int v = 0; v < PP; v++) if (S >> v & 1) p += fw[S][v]; P[S] = p; Q[S] = 0; }
    for (int k = 0; k <= PF; k++){
      int N = k ^ (k >> 1);
      if (k){ int v = __builtin_ctz(k); int gain = N >> v & 1;   /* bit v of N flipped */
        for (int S = 1; S <= PF; S++) if (S >> v & 1){
          if (gain){ P[S] -= fw[S][v]; Q[S] += bw[S][v]; } else { P[S] += fw[S][v]; Q[S] -= bw[S][v]; } } }
      u64 h11 = P[PF] + Q[PF];
      for (int S = 1; S < PF; S++) h11 += P[S] * Q[PF & ~S];
#ifdef CHECK
      if ((parents % 97) == 1 && (N % 37) == 0){
        memcpy(outm, o10, sizeof o10); memcpy(inm, i10, sizeof i10); outm[PP] = inm[PP] = 0;
        for (int v = 0; v < PP; v++){ if (N >> v & 1) setarc(PP, v); else setarc(v, PP); }
        u32 hd = hpaths(PP + 1); nchk++;
        if (hd != h11){ fprintf(stderr, "SPLIT MISMATCH code %llu N %d: %llu vs %u\n", (unsigned long long)c, N, (unsigned long long)h11, hd); return 2; }
      }
#endif
      if (h11 >= (u64)LIM) continue;
      kept11++;
      memcpy(outm, o10, sizeof o10); memcpy(inm, i10, sizeof i10); outm[PP] = inm[PP] = outm[PP + 1] = inm[PP + 1] = 0;
      for (int v = 0; v < PP; v++){ if (N >> v & 1) setarc(PP, v); else setarc(v, PP); }
      for (int w = 0; w < PP; w++){ if (N >> (w ^ 1) & 1) setarc(w, PP + 1); else setarc(PP + 1, w); }
      setarc(PP, PP + 1);
#ifdef CHECK
      for (int u = 0; u < NN; u++) for (int v = 0; v < NN; v++) if (u != v){
        int uv = outm[u] >> v & 1, vu = outm[v] >> u & 1;
        if (uv + vu != 1 || uv != ((outm[sig12[v]] >> sig12[u]) & 1)){ fprintf(stderr, "ANTI FAIL\n"); return 3; } }
#endif
      u32 h12 = hpaths(NN);
      if (h12 < (u32)LIM){ cnt[h12]++; if (landau_strong(NN)) scnt[h12]++; }
    }
  }
  fprintf(stderr, "codes [%llu,%llu): parents with H < %d: %ld; order-11 intermediates kept: %ld; split checks: %ld\n",
          (unsigned long long)C0, (unsigned long long)C1, LIM, parents, kept11, nchk);
  for (int h = 0; h < (1 << 18) && h < LIM; h++) if (cnt[h]) printf("%d %llu %llu\n", h, (unsigned long long)cnt[h], (unsigned long long)scnt[h]);
  (void)sig12;
  return 0;
}
