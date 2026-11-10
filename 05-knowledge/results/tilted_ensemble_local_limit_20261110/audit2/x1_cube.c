// x1_cube.c -- AUDIT 2 (auditor's own code; nothing copied from the packet's biasweights.c / oddedge.c / *.py).
// Exact mod-4 statistics of every vertex cube {R + v_X : X subset of V(R)} of the tournaments R read from stdin,
// the parity law (4.5), the cut quadratic (4.2), Theorem C (3), the polar form (4.4), THM-4177 (16), the gate formula
// (late item L1), the digit-degree law (late item L3), with labelled weights m!/|Aut(R)|.
//
// usage:  gentourng -q m [res/mod] | x1_cube m direct [flags]     H(R+v_X) for EVERY X by a fresh subset DP on the
//                                                                 (m+1)-vertex tournament (no cut formula anywhere);
//                                                                 mu, a from second/first differences of those counts
//         gentourng -q m [res/mod] | x1_cube m cut K [flags]      H(R+v_X) from the insertion decomposition (4.1) with
//                                                                 exact integers; K random X per class are recounted
//                                                                 by the fresh (m+1)-vertex DP as a control
//         x1_cube m labelled r k [flags]                          all 2^C(m,2) labelled tournaments with code = r mod k,
//                                                                 direct mode, weight 1 (no automorphism groups used)
// flags (one word, any of): p  check (4.4) beta = P_ij + P_ji mod 2 and THM-4177 (16) mu = 2 sum_(even paths) H(R - P)
//                           g  check the gate formula H(Gamma R) = H(R) + P_R(2) and P_R(2) = 0 mod 4 (direct DP, m+2 vertices)
//                           d  F_2-degree of the binary digits 1, 2, 3 of X -> H(R + v_X)
// Output: one line of raw integer totals "key value ..." (parts can be added; see x1_merge.py).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64; typedef unsigned __int128 u128;
#define MAXM 10
static int m, FULL;
static int adj[MAXM][MAXM], inR[MAXM], outR[MAXM];
static i64 Ee[1 << MAXM][MAXM], Bs[1 << MAXM][MAXM], HS[1 << MAXM], Hx[1 << MAXM];
static i64 Dp[1 << (MAXM + 2)][MAXM + 2];
static i64 Fse[1 << MAXM][MAXM];

static void dpR(void) {                       // Ee[S][v]: Ham paths of R[S] ending at v; Bs: starting at v; HS[S] = H(R[S])
  HS[0] = 1;
  for (int S = 1; S <= FULL; S++) { i64 tot = 0;
    for (int v = 0; v < m; v++) { if (!(S >> v & 1)) { Ee[S][v] = Bs[S][v] = 0; continue; }
      int rest = S ^ (1 << v); if (!rest) { Ee[S][v] = Bs[S][v] = 1; tot += 1; continue; }
      i64 e = 0, b = 0;
      for (int pp = rest & inR[v]; pp; pp &= pp - 1) e += Ee[rest][__builtin_ctz(pp)];
      for (int pp = rest & outR[v]; pp; pp &= pp - 1) b += Bs[rest][__builtin_ctz(pp)];
      Ee[S][v] = e; Bs[S][v] = b; tot += e; }
    HS[S] = tot; }
}
static i64 ham_fresh(int n1, const int *inm) { // Hamiltonian paths of an n1-vertex tournament given by in-neighbour masks
  int full1 = (1 << n1) - 1;
  for (int S = 1; S <= full1; S++) for (int vv = S; vv; vv &= vv - 1) { int v = __builtin_ctz(vv), rest = S ^ (1 << v);
      if (!rest) { Dp[S][v] = 1; continue; } i64 e = 0;
      for (int pp = rest & inm[v]; pp; pp &= pp - 1) e += Dp[rest][__builtin_ctz(pp)]; Dp[S][v] = e; }
  i64 h = 0; for (int v = 0; v < n1; v++) h += Dp[full1][v]; return h;
}
static i64 direct_H(int X) {                  // fresh count for R + v_X (v = vertex m, v -> X, V \ X -> v)
  int inm[MAXM + 2];
  for (int a = 0; a < m; a++) inm[a] = inR[a] | ((X >> a & 1) ? (1 << m) : 0);
  inm[m] = (~X) & FULL; return ham_fresh(m + 1, inm);
}
static i64 gate_H(void) {                     // Gamma(R): a = vertex m, b = vertex m+1 ; a -> b, b -> R, R -> a
  int inm[MAXM + 2];
  for (int x = 0; x < m; x++) inm[x] = inR[x] | (1 << (m + 1));
  inm[m] = FULL; inm[m + 1] = 1 << m; return ham_fresh(m + 2, inm);
}
static void fill_cut(void) {                  // Hx[X] from (4.1): position of v in the path (first / last / between two paths)
  static i64 Mx[MAXM][MAXM];
  for (int j = 0; j < m; j++) for (int i = 0; i < m; i++) { Mx[j][i] = 0; if (i == j) continue;
      int rest = FULL & ~(1 << i) & ~(1 << j); i64 s = 0;
      for (int S0 = rest;; S0 = (S0 - 1) & rest) { int S = S0 | (1 << j); s += Ee[S][j] * Bs[FULL ^ S][i]; if (!S0) break; }
      Mx[j][i] = s; }
  Hx[0] = HS[FULL];
  for (int X = 1; X <= FULL; X++) { int t = __builtin_ctz(X), Y = X & (X - 1);
    i64 d = Bs[FULL][t] - Ee[FULL][t];
    for (int j = 0; j < m; j++) { if (j == t) continue; if (Y >> j & 1) d -= Mx[t][j]; else d += Mx[j][t]; }
    Hx[X] = Hx[Y] + d; }
}
static int sc[MAXM], pm[MAXM], usedm;
static i64 aut_rec(int k) { if (k == m) return 1; i64 c = 0;
  for (int w = 0; w < m; w++) { if ((usedm >> w & 1) || sc[w] != sc[k]) continue; int ok = 1;
    for (int j = 0; j < k && ok; j++) if (adj[j][k] != adj[pm[j]][w]) ok = 0;
    if (ok) { pm[k] = w; usedm |= 1 << w; c += aut_rec(k + 1); usedm &= ~(1 << w); } }
  return c; }

// affine spans over F_2 (meaningful in labelled mode: the set of labelled tournaments is closed under relabelling)
static u64 baseB, baseLB, pivB[64], pivLB[64]; static int haveBase = 0, dimB = 0, dimLB = 0;
static void span_add(u64 v, u64 *base, u64 *piv, int *dim, int first) { if (first) { *base = v; return; } u64 d = v ^ *base;
  while (d) { int h = 63 - __builtin_clzll(d); if (piv[h]) d ^= piv[h]; else { piv[h] = d; (*dim)++; return; } } }

static u64 nt = 0, cls0[MAXM + 2], cls1[MAXM + 2], edgehist[64], deghist[4][MAXM + 2];
static i64 Wtot = 0, S1 = 0, Usum = 0, Asum = 0, lab0[MAXM + 2], lab1[MAXM + 2], w8[8], wEven0 = 0, wEven1 = 0, minmu = (1LL << 62), maxmu = 0;
static u128 Qsum = 0;
static u64 f_h0 = 0, f_aodd = 0, f_muodd = 0, f_mulow = 0, f_42 = 0, f_c3 = 0, f_45a = 0, f_45b = 0, f_45c = 0, f_edges = 0, f_ell = 0,
           f_44 = 0, f_4177 = 0, f_ctrl = 0, f_evenH = 0, f_aut = 0, n_ctrl = 0, n_44 = 0, f_evenpaths = 0, f_anti = 0,
           f_gate = 0, f_gate4 = 0, n_gate = 0, f_deg = 0, n_deg = 0;
static u64 rng = 88172645463325252ULL;
static u64 xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

static void process(int mode, int K, int do44, int dogate, int dodeg, i64 weight) {   // mode 0 direct, 1 cut ; weight < 0 : use m!/|Aut|
  for (int i = 0; i < m; i++) { inR[i] = outR[i] = 0; sc[i] = 0; }
  for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) if (i != j && adj[i][j]) { outR[i] |= 1 << j; inR[j] |= 1 << i; sc[i]++; }
  dpR(); nt++;
  if (mode == 0) { for (int X = 0; X <= FULL; X++) Hx[X] = direct_H(X); }
  else { fill_cut(); for (int c = 0; c < K; c++) { int X = (int)(xr() >> 20) & FULL; n_ctrl++; if (direct_H(X) != Hx[X]) f_ctrl++; } }
  i64 h0 = Hx[0]; if (h0 != HS[FULL] || Hx[FULL] != h0) f_h0++;
  i64 a[MAXM], mu[MAXM][MAXM], summu = 0, suma = 0, split = 0; int rows[MAXM], ell[MAXM], edges = 0, ellw = 0;
  for (int i = 0; i < m; i++) { a[i] = Hx[1 << i] - h0; suma += a[i]; if (a[i] & 1) f_aodd++; ell[i] = (int)((a[i] >> 1) & 1); ellw += ell[i]; rows[i] = 0; }
  for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { i64 u = Hx[1 << i] + Hx[1 << j] - Hx[(1 << i) | (1 << j)] - h0; mu[i][j] = mu[j][i] = u;
      summu += u; if (u & 1) f_muodd++; if (u < 2 * HS[FULL & ~(1 << i) & ~(1 << j)]) f_mulow++; if (u < minmu) minmu = u; if (u > maxmu) maxmu = u;
      if ((u >> 1) & 1) { rows[i] |= 1 << j; rows[j] |= 1 << i; edges++; } }
  for (int S = 1; S < FULL; S++) split += HS[S] * HS[FULL ^ S];
  { i64 ev = 0; for (int S = 1; S <= FULL; S++) if (!(__builtin_popcount(S) & 1)) ev += HS[S]; if (m >= 2 && !(ev & 1)) f_evenpaths++; }
  if (summu != suma) f_45a++; if (suma != split) f_45b++; if ((summu & 3) != 2) f_45c++; if (!(edges & 1)) f_edges++; if (!(ellw & 1)) f_ell++;
  edgehist[edges]++;
  if (mode == 0) {                             // cut quadratic (4.2) on every X, from the directly counted table
    for (int X = 0; X <= FULL; X++) { i64 p = h0; for (int i = 0; i < m; i++) if (X >> i & 1) { p += a[i]; for (int j = i + 1; j < m; j++) if (X >> j & 1) p -= mu[i][j]; }
      if (p != Hx[X]) { f_42++; break; } } }
  // mod-4 tallies from the table of counts (no quadratic form used)
  i64 c4[4] = {0, 0, 0, 0}, c8[8] = {0}; for (int X = 0; X <= FULL; X++) { c4[Hx[X] & 3]++; c8[Hx[X] & 7]++; }
  if (c4[0] || c4[2]) f_evenH++;
  i64 diff = c4[h0 & 3] - c4[(h0 + 2) & 3], ad = diff < 0 ? -diff : diff;
  // rank and radical of beta
  int rr[MAXM], id[MAXM], rk = 0; for (int i = 0; i < m; i++) { rr[i] = rows[i]; id[i] = 1 << i; }
  for (int c = 0; c < m; c++) { int p = -1; for (int i = rk; i < m; i++) if (rr[i] >> c & 1) { p = i; break; } if (p < 0) continue;
    int t = rr[p]; rr[p] = rr[rk]; rr[rk] = t; t = id[p]; id[p] = id[rk]; id[rk] = t;
    for (int i = 0; i < m; i++) if (i != rk && (rr[i] >> c & 1)) { rr[i] ^= rr[rk]; id[i] ^= id[rk]; } rk++; }
  int radzero = 1;                              // q read off the table: q(X) = (Hx[X] - h0)/2 mod 2 ; test a basis of the radical
  for (int i = rk; i < m; i++) { int X = id[i]; int bad = 0; for (int k2 = 0; k2 < m; k2++) if (__builtin_popcount(rows[k2] & X) & 1) bad = 1;
    if (bad) f_c3 += 1000000; if (((Hx[X] - h0) >> 1) & 1) radzero = 0; }
  if ((rk & 1) || !((ad == 0 && !radzero) || (ad == (1LL << (m - rk / 2)) && radzero))) f_c3++;
  int alleven = 1; for (int i = 0; i < m; i++) if (__builtin_popcount(rows[i]) & 1) alleven = 0;
  { i64 s = 0; for (int X = 0; X <= FULL; X++) s += (((Hx[X] - Hx[FULL ^ X]) >> 1) & 1) ? -1 : 1;   // antipodal pair X, V \ X
    if (s != (alleven ? (i64)FULL + 1 : 0)) f_anti++; }
  if (do44) {                                   // (4.4) and THM-4177 (16): paths between i and j on an even number of vertices
    static i64 P[MAXM][MAXM], Pw[MAXM][MAXM]; memset(P, 0, sizeof P); memset(Pw, 0, sizeof Pw);
    for (int s = 0; s < m; s++) { for (int S = 1; S <= FULL; S++) { if (!(S >> s & 1)) continue;
        if (S == (1 << s)) { for (int v = 0; v < m; v++) Fse[S][v] = 0; Fse[S][s] = 1; continue; }
        int ev = !(__builtin_popcount(S) & 1);
        for (int v = 0; v < m; v++) { if (!(S >> v & 1) || v == s) { Fse[S][v] = 0; continue; } int rest = S ^ (1 << v); i64 e = 0;
          for (int pp = rest & inR[v]; pp; pp &= pp - 1) e += Fse[rest][__builtin_ctz(pp)]; Fse[S][v] = e; if (ev) { P[s][v] += e; Pw[s][v] += e * HS[FULL ^ S]; } } } }
    n_44++; int bad = 0, bad2 = 0; for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if (((P[i][j] + P[j][i]) & 1) != (rows[i] >> j & 1)) bad = 1;
        if (mu[i][j] != 2 * (Pw[i][j] + Pw[j][i])) bad2 = 1; }
    for (int i = 0; i < m; i++) { i64 l = 0; for (int j = 0; j < m; j++) if (j != i) l += P[i][j]; if ((l & 1) != ell[i]) bad = 1; }
    if (bad) f_44++; if (bad2) f_4177++; }
  if (dogate) { i64 P2 = 0; for (int S = 0; S <= FULL; S++) P2 += HS[S] * HS[FULL ^ S]; n_gate++;
    if (gate_H() != h0 + P2) f_gate++; if (m >= 2 && (P2 & 3)) f_gate4++; }
  if (dodeg) { static unsigned char bt[1 << MAXM]; n_deg++;
    for (int k = 1; k <= 3; k++) { for (int X = 0; X <= FULL; X++) bt[X] = (unsigned char)((Hx[X] >> k) & 1);
      for (int i = 0; i < m; i++) for (int X = 0; X <= FULL; X++) if (X >> i & 1) bt[X] ^= bt[X ^ (1 << i)];
      int dg = -1; for (int X = 0; X <= FULL; X++) if (bt[X] && __builtin_popcount(X) > dg) dg = __builtin_popcount(X);
      deghist[k][dg + 1]++; if (dg > (1 << k)) f_deg++; } }
  { u64 vb = 0, vl = 0; int b = 0; for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if (rows[i] >> j & 1) vb |= 1ULL << b; b++; }
    for (int i = 0; i < m; i++) if (ell[i]) vl |= 1ULL << (b + i);
    span_add(vb, &baseB, pivB, &dimB, !haveBase); span_add(vb | vl, &baseLB, pivLB, &dimLB, !haveBase); haveBase = 1; }
  i64 w = weight;
  if (weight < 0) { i64 mf = 1; for (int i = 2; i <= m; i++) mf *= i; usedm = 0; i64 au = aut_rec(0); if (au <= 0 || mf % au) { f_aut++; return; } w = mf / au; }
  Wtot += w; S1 += w * (c4[1] - c4[3]); Usum += w * diff; Asum += w * ad; Qsum += (u128)w * (u128)(ad * ad);
  for (int r = 0; r < 8; r++) w8[r] += w * c8[r];
  if (ad) { cls1[rk]++; lab1[rk] += w; if (alleven) wEven1 += w; } else { cls0[rk]++; lab0[rk] += w; if (alleven) wEven0 += w; }
}
static void pr128(u128 x) { char buf[64]; int k = 0; if (!x) { printf("0"); return; } while (x) { buf[k++] = (char)('0' + (int)(x % 10)); x /= 10; } while (k) putchar(buf[--k]); }

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: x1_cube m direct [flags] | m cut K [flags] | m labelled r k [flags]\n"); return 2; }
  m = atoi(argv[1]); FULL = (1 << m) - 1; if (m < 2 || m > MAXM) return 2;
  int islab = !strcmp(argv[2], "labelled"), iscut = !strcmp(argv[2], "cut"); int nfix = islab ? 5 : (iscut ? 4 : 3);
  const char *flags = argc > nfix ? argv[nfix] : ""; int do44 = strchr(flags, 'p') != 0, dogate = strchr(flags, 'g') != 0, dodeg = strchr(flags, 'd') != 0;
  if (argc < nfix) return 2;
  if (islab) { int r = atoi(argv[3]), k = atoi(argv[4]);
    int np = m * (m - 1) / 2; for (u64 code = 0; code < (1ULL << np); code++) { if ((int)(code % (u64)k) != r) continue; int b = 0;
      for (int i = 0; i < m; i++) { adj[i][i] = 0; for (int j = i + 1; j < m; j++) { int x = (int)(code >> b & 1); b++; adj[i][j] = x; adj[j][i] = !x; } }
      process(0, 0, do44, dogate, dodeg, 1); } }
  else { int K = iscut ? atoi(argv[3]) : 0; char line[256];
    while (fgets(line, sizeof line, stdin)) { if (line[0] != '0' && line[0] != '1') continue; int b = 0;
      for (int i = 0; i < m; i++) { adj[i][i] = 0; for (int j = i + 1; j < m; j++) { int x = line[b++] == '1'; adj[i][j] = x; adj[j][i] = !x; } }
      process(iscut, K, do44, dogate, dodeg, -1); } }
  printf("order %d mode %s classes %llu W %lld S1 %lld U %lld A %lld Q ", m, argv[2], nt, Wtot, S1, Usum, Asum); pr128(Qsum);
  printf(" wEven0 %lld wEven1 %lld minmu %lld maxmu %lld affdimB %d affdimLB %d", wEven0, wEven1, minmu, maxmu, dimB, dimLB);
  printf(" fail_h0 %llu fail_aodd %llu fail_muodd %llu fail_mulow %llu fail_42 %llu fail_c3 %llu fail_45a %llu fail_45b %llu fail_45c %llu fail_edges %llu fail_ell %llu fail_evenpaths %llu fail_anti %llu fail_44 %llu fail_4177 %llu n_44 %llu fail_ctrl %llu n_ctrl %llu fail_evenH %llu fail_aut %llu fail_gate %llu fail_gate4 %llu n_gate %llu fail_deg %llu n_deg %llu",
         f_h0, f_aodd, f_muodd, f_mulow, f_42, f_c3, f_45a, f_45b, f_45c, f_edges, f_ell, f_evenpaths, f_anti, f_44, f_4177, n_44, f_ctrl, n_ctrl, f_evenH, f_aut, f_gate, f_gate4, n_gate, f_deg, n_deg);
  printf(" mod8"); for (int r = 0; r < 8; r++) printf(" %lld", w8[r]);
  for (int r = 0; r <= m; r++) if (cls0[r] || cls1[r]) printf(" rank %d %llu %llu %lld %lld", r, cls0[r], cls1[r], lab0[r], lab1[r]);
  printf(" edges"); for (int e = 0; e < 64; e++) if (edgehist[e]) printf(" %d:%llu", e, edgehist[e]);
  for (int k = 1; k <= 3; k++) { printf(" deg%d", k); for (int d = 0; d <= m + 1; d++) if (deghist[k][d]) printf(" %d:%llu", d - 1, deghist[k][d]); }
  printf("\n"); return 0;
}
