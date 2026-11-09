/* aud_sc13.c - audit K2 of THM-4930 statement 5 (2026-11-08). Written from scratch by the auditor (extends audit K1's
 * aud_sc12.c; shares no code with the packet's sc_odd13.c).
 *
 * All self-converse tournaments of the even order NE = NO-1 and of the odd order NO with H < LIM.
 * Labelling: sigma(v) = v ^ 1 on {0, ..., NO-2} (pairs {2k, 2k+1}) and sigma(NO-1) = NO-1 (the fixed vertex v0).
 * Reduction (audit K1): conjugating by the transpositions (2k 2k+1), which commute with sigma, every within-pair arc is
 * taken as 2k -> 2k+1; every class still occurs.
 * Decomposition (base order PB = NO-3, pair a = PB, b = PB+1, fixed vertex v0 = NO-1):
 *   T_NO - v0 = T_NE is sigma-anti-invariant of even order NE; T_NE - {a, b} = T_PB is a reduced code of order PB.
 *   Insertion monotonicity (H(T - v) <= H(T)) gives H(T_PB) <= H(T_NE - b) <= H(T_NE) <= H(T_NO), so pruning every level
 *   at H < LIM loses nothing.
 *   level 1: each reduced base code with H < LIM (aud_sc.c's construction);
 *   level 2: vertex a, all 2^PB out-masks N, Gray-code order, H(T_PB + a) by splitting paths at a (forward/backward
 *            subset tables of the base);
 *   level 3: vertex b = sigma(a) with forced arcs (w -> b iff a -> sigma(w)) and a -> b; H(T_NE) by a full pull DP, whose
 *            table F is the forward table of T_NE;
 *   level 4: the fixed vertex v0, all 2^(NE/2) orbit masks (v0 -> 2k iff 2k+1 -> v0), Gray-code order, H(T_NO) by
 *            splitting at v0. The backward table of T_NE is read off F by anti-invariance: G[S][w] = F[sigma S][sigma w].
 * Strongness by Landau's score criterion. With CHECK, sampled level-2 and level-4 values are recomputed by direct DP and
 * the order-NE and order-NO tournaments are checked to be tournaments and sigma-anti-invariant.
 * Output: "E h count strongcount" (order NE) and "O h count strongcount" (order NO), counts of reduced labelled codes.
 * Validation: with LIM above max H, the E and O lines must equal aud_sc.c's reduced censuses of orders NE and NO.
 * usage: aud_sc13 NO LIM [C0 C1]     (NO odd, 5 <= NO <= 13; C0 <= base code < C1)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef uint32_t u32; typedef uint64_t u64;
#define MAXN 14
#define HCAP (1 << 17)
static int n, sig[MAXN];
static int oi[128], oj[128], nother, pi_[16], pj_[16], npair;
static uint16_t outm[MAXN], inm[MAXN];
static u32 *F, *F2;
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
static void buildbase(u64 c){
  u64 code = c | (((1ULL << npair) - 1) << nother);
  memset(outm, 0, sizeof outm); memset(inm, 0, sizeof inm);
  for (int k = 0; k < nother; k++){ int a = oi[k], b = oj[k]; if (!((code >> k) & 1)){ int t = a; a = b; b = t; } setarc(a, b); setarc(sig[b], sig[a]); }
  for (int k = 0; k < npair; k++){ int a = pi_[k], b = pj_[k]; if (!((code >> (nother + k)) & 1)){ int t = a; a = b; b = t; } setarc(a, b); }
}
static u32 hpaths(u32 *T, int nn){            /* pull DP: T[S*nn+v] = paths covering S ending at v */
  int FN = 1 << nn;
  for (int v = 0; v < nn; v++) T[(1u << v) * nn + v] = 1;
  for (int S = 3; S < FN; S++){ if (!(S & (S - 1))) continue;
    for (int v = 0; v < nn; v++){
      if (!(S >> v & 1)){ T[(u32)S * nn + v] = 0; continue; }
      int R = S & ~(1 << v); uint16_t m = inm[v] & (uint16_t)R; u32 s = 0;
      while (m){ int u = __builtin_ctz(m); m &= m - 1; s += T[(u32)R * nn + u]; }
      T[(u32)S * nn + v] = s; } }
  u32 h = 0; for (int v = 0; v < nn; v++) h += T[(u32)(FN - 1) * nn + v]; return h;
}
static int landau(const int *sc0, int nn){
  int sc[MAXN]; memcpy(sc, sc0, sizeof(int) * nn);
  for (int i = 1; i < nn; i++){ int x = sc[i], j = i - 1; while (j >= 0 && sc[j] > x){ sc[j + 1] = sc[j]; j--; } sc[j + 1] = x; }
  int s = 0; for (int k = 1; k < nn; k++){ s += sc[k - 1]; if (s == k * (k - 1) / 2) return 0; }
  return 1;
}
static int anti_ok(int nn, const int *sg){
  for (int u = 0; u < nn; u++) for (int v = 0; v < nn; v++) if (u != v){
    int uv = outm[u] >> v & 1, vu = outm[v] >> u & 1;
    if (uv + vu != 1 || uv != ((outm[sg[v]] >> sg[u]) & 1)) return 0; }
  return 1;
}
static u32 fw[1 << 10][10], bw[1 << 10][10];
static u64 P[1 << 10], Q[1 << 10], P13[1 << 12], Q13[1 << 12];
int main(int argc, char **argv){
  if (argc < 3){ fprintf(stderr, "usage: aud_sc13 NO LIM [C0 C1]\n"); return 1; }
  int NO = atoi(argv[1]); long LIM = atol(argv[2]);
  if (NO % 2 == 0 || NO < 5 || NO > 13){ fprintf(stderr, "NO must be odd, 5..13\n"); return 1; }
  int NE = NO - 1, PB = NO - 3, A = PB, B = PB + 1, V0 = NO - 1;
  int PF = (1 << PB) - 1, EF = (1 << NE) - 1, MP = NE / 2;
  orbits(PB);
  u64 C0 = argc > 3 ? strtoull(argv[3], 0, 10) : 0, C1 = argc > 4 ? strtoull(argv[4], 0, 10) : (1ULL << nother);
  F = calloc((size_t)(1u << 13) * 13, sizeof(u32)); F2 = calloc((size_t)(1u << 13) * 13, sizeof(u32));
  u64 *cE = calloc(HCAP, sizeof(u64)), *sE = calloc(HCAP, sizeof(u64)), *cO = calloc(HCAP, sizeof(u64)), *sO = calloc(HCAP, sizeof(u64));
  int sgE[MAXN], sgO[MAXN]; for (int v = 0; v < NE; v++) sgE[v] = sgO[v] = v ^ 1; sgO[V0] = V0;
  u32 alt = 0; for (int k = 0; k < MP; k++) alt |= 1u << (2 * k);           /* even positions mask within NE bits */
  long parents = 0, lev2 = 0, lev3 = 0, chk2 = 0, chk4 = 0;
  for (u64 c = C0; c < C1; c++){
    orbits(PB); buildbase(c);
    u32 hb = hpaths(F, PB);
    if (hb >= (u64)LIM) continue;
    parents++;
    uint16_t ob[MAXN], ib[MAXN]; memcpy(ob, outm, sizeof ob); memcpy(ib, inm, sizeof ib);
    for (int S = 1; S <= PF; S++) for (int v = 0; v < PB; v++){
      if (!(S >> v & 1)){ fw[S][v] = bw[S][v] = 0; continue; }
      if (S == (1 << v)){ fw[S][v] = bw[S][v] = 1; continue; }
      int R = S & ~(1 << v); u32 x = 0, y = 0;
      uint16_t mi = ib[v] & R, mo = ob[v] & R;
      while (mi){ int u = __builtin_ctz(mi); mi &= mi - 1; x += fw[R][u]; }
      while (mo){ int u = __builtin_ctz(mo); mo &= mo - 1; y += bw[R][u]; }
      fw[S][v] = x; bw[S][v] = y;
    }
    P[0] = Q[0] = 0;
    for (int S = 1; S <= PF; S++){ u64 p = 0; for (int v = 0; v < PB; v++) if (S >> v & 1) p += fw[S][v]; P[S] = p; Q[S] = 0; }
    for (int k = 0; k <= PF; k++){
      int N = k ^ (k >> 1);                      /* out-neighbours of a in the base */
      if (k){ int v = __builtin_ctz(k); int gain = N >> v & 1;
        for (int S = 1; S <= PF; S++) if (S >> v & 1){
          if (gain){ P[S] -= fw[S][v]; Q[S] += bw[S][v]; } else { P[S] += fw[S][v]; Q[S] -= bw[S][v]; } } }
      u64 h11 = P[PF] + Q[PF];
      for (int S = 1; S < PF; S++) h11 += P[S] * Q[PF & ~S];
      memcpy(outm, ob, sizeof ob); memcpy(inm, ib, sizeof ib);
      for (int v = 0; v < PB; v++){ if (N >> v & 1) setarc(A, v); else setarc(v, A); }
#ifdef CHECK
      if ((parents % 53) == 1 && (N % 29) == 0){
        u32 hd = hpaths(F2, PB + 1); chk2++;
        if (hd != h11){ fprintf(stderr, "LEVEL-2 MISMATCH base %llu N %d: %llu vs %u\n", (unsigned long long)c, N, (unsigned long long)h11, hd); return 2; }
      }
#endif
      if (h11 >= (u64)LIM) continue;
      lev2++;
      for (int w = 0; w < PB; w++){ if (N >> (w ^ 1) & 1) setarc(w, B); else setarc(B, w); }
      setarc(A, B);
#ifdef CHECK
      if (!anti_ok(NE, sgE)){ fprintf(stderr, "ORDER-NE ANTI FAIL\n"); return 3; }
#endif
      u32 h12 = hpaths(F, NE);
      if (h12 >= (u64)LIM) continue;
      if (h12 >= HCAP){ fprintf(stderr, "HCAP (order NE)\n"); return 6; }
      lev3++;
      int scE[MAXN]; for (int v = 0; v < NE; v++) scE[v] = __builtin_popcount(outm[v]);
      cE[h12]++; if (landau(scE, NE)) sE[h12]++;
      /* level 4: fixed vertex v0. Mask z (MP bits): z_k = 1 means v0 -> 2k and 2k+1 -> v0; z_k = 0: 2k -> v0, v0 -> 2k+1 */
      P13[0] = Q13[0] = 0;
      for (int S = 1; S <= EF; S++){
        u32 sS = ((S & alt) << 1) | ((S >> 1) & alt);
        u64 p = 0, q = 0;
        for (int v = 0; v < NE; v++) if (S >> v & 1){
          if (!(v & 1)) p += F[(u32)S * NE + v];                   /* even v -> v0 */
          else q += F[(u32)sS * NE + (v ^ 1)];                     /* v0 -> odd w: G[S][w] = F[sigma S][sigma w] */
        }
        P13[S] = p; Q13[S] = q;
      }
      uint16_t oE[MAXN], iE[MAXN]; memcpy(oE, outm, sizeof oE); memcpy(iE, inm, sizeof iE);
      for (int k = 0; k < (1 << MP); k++){
        int z = k ^ (k >> 1);
        if (k){ int t = __builtin_ctz(k); int on = z >> t & 1; int ev = 2 * t, od = 2 * t + 1;
          for (int S = 1; S <= EF; S++){
            u32 sS = 0; int hasE = S >> ev & 1, hasO = S >> od & 1;
            if (!hasE && !hasO) continue;
            sS = ((S & alt) << 1) | ((S >> 1) & alt);
            if (hasE){ u64 f = F[(u32)S * NE + ev], g = F[(u32)sS * NE + od];       /* G[S][ev] = F[sS][od] */
              if (on){ P13[S] -= f; Q13[S] += g; } else { P13[S] += f; Q13[S] -= g; } }
            if (hasO){ u64 f = F[(u32)S * NE + od], g = F[(u32)sS * NE + ev];       /* G[S][od] = F[sS][ev] */
              if (on){ P13[S] += f; Q13[S] -= g; } else { P13[S] -= f; Q13[S] += g; } }
          } }
        u64 h13 = P13[EF] + Q13[EF];
        for (int S = 1; S < EF; S++) h13 += P13[S] * Q13[EF & ~S];
#ifdef CHECK
        if ((lev3 % 211) == 1){
          memcpy(outm, oE, sizeof oE); memcpy(inm, iE, sizeof iE); outm[V0] = inm[V0] = 0;
          for (int t = 0; t < MP; t++){ if (z >> t & 1){ setarc(V0, 2 * t); setarc(2 * t + 1, V0); } else { setarc(2 * t, V0); setarc(V0, 2 * t + 1); } }
          if (!anti_ok(NO, sgO)){ fprintf(stderr, "ORDER-NO ANTI FAIL\n"); return 4; }
          u32 hd = hpaths(F2, NO); chk4++;
          if (hd != h13){ fprintf(stderr, "LEVEL-4 MISMATCH base %llu N %d z %d: %llu vs %u\n", (unsigned long long)c, N, z, (unsigned long long)h13, hd); return 5; }
        }
#endif
        if (h13 >= (u64)LIM) continue;
        int scO[MAXN]; int s0 = 0;
        for (int t = 0; t < MP; t++){
          if (z >> t & 1){ scO[2 * t] = scE[2 * t]; scO[2 * t + 1] = scE[2 * t + 1] + 1; s0++; }
          else { scO[2 * t] = scE[2 * t] + 1; scO[2 * t + 1] = scE[2 * t + 1]; s0++; } }
        scO[V0] = s0;
        if (h13 >= HCAP){ fprintf(stderr, "HCAP\n"); return 6; }
        cO[h13]++; if (landau(scO, NO)) sO[h13]++;
      }
    }
  }
  fprintf(stderr, "NO=%d LIM=%ld codes [%llu,%llu): parents %ld, order-%d intermediates %ld, order-%d tournaments %ld, checks %ld/%ld\n",
          NO, LIM, (unsigned long long)C0, (unsigned long long)C1, parents, PB + 1, lev2, NE, lev3, chk2, chk4);
  for (long h = 0; h < HCAP && h < LIM; h++) if (cE[h]) printf("E %ld %llu %llu\n", h, (unsigned long long)cE[h], (unsigned long long)sE[h]);
  for (long h = 0; h < HCAP && h < LIM; h++) if (cO[h]) printf("O %ld %llu %llu\n", h, (unsigned long long)cO[h], (unsigned long long)sO[h]);
  return 0;
}
