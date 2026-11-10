// x6_affine.c -- AUDIT 2, late item L2 (auditor's own code).  Affine span over F_2 of the vectors beta_R in F_2^C(m,2)
// and of the pairs (ell_R, beta_R) in F_2^(m + C(m,2)), R over labelled tournaments of order m.
// The upper bounds C(m,2) - 1 and m + C(m,2) - 2 are Lemma 4.5 (checked on every tournament used); a lower bound on a
// span is a certificate and needs no complete sample.
// usage: x6_affine m all            every labelled tournament (m <= 7)
//        x6_affine m random N seed  N pseudo-random labelled tournaments (m <= 11)
// (ell, beta) from the insertion decomposition (4.1) (validated against direct counts by x1_cube on all classes).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64; typedef unsigned __int128 u128;
#define MAXM 11
static int m, FULL, inR[MAXM], outR[MAXM];
static i64 Ee[1 << MAXM][MAXM], Bs[1 << MAXM][MAXM];
static u128 baseB, baseLB, pivB[128], pivLB[128]; static int have = 0, dimB = 0, dimLB = 0;
static int hibit(u128 d) { u64 hi = (u64)(d >> 64); return hi ? 127 - __builtin_clzll(hi) : 63 - __builtin_clzll((u64)d); }
static void span_add(u128 v, u128 *base, u128 *piv, int *dim, int first) { if (first) { *base = v; return; } u128 d = v ^ *base;
  while (d) { int h = hibit(d); if (piv[h]) d ^= piv[h]; else { piv[h] = d; (*dim)++; return; } } }
static u64 viol = 0, cnt = 0;
static void one(u64 code) {
  int b = 0; for (int i = 0; i < m; i++) inR[i] = outR[i] = 0;
  for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if (code >> b & 1) { outR[i] |= 1 << j; inR[j] |= 1 << i; } else { outR[j] |= 1 << i; inR[i] |= 1 << j; } b++; }
  for (int S = 1; S <= FULL; S++) for (int vv = S; vv; vv &= vv - 1) { int v = __builtin_ctz(vv), rest = S ^ (1 << v);
      if (!rest) { Ee[S][v] = Bs[S][v] = 1; continue; } i64 e = 0, bb = 0;
      for (int pp = rest & inR[v]; pp; pp &= pp - 1) e += Ee[rest][__builtin_ctz(pp)];
      for (int pp = rest & outR[v]; pp; pp &= pp - 1) bb += Bs[rest][__builtin_ctz(pp)]; Ee[S][v] = e; Bs[S][v] = bb; }
  static i64 Mx[MAXM][MAXM];
  for (int j = 0; j < m; j++) for (int i = 0; i < m; i++) { Mx[j][i] = 0; if (i == j) continue; int rest = FULL & ~(1 << i) & ~(1 << j); i64 s = 0;
      for (int S0 = rest;; S0 = (S0 - 1) & rest) { int S = S0 | (1 << j); s += Ee[S][j] * Bs[FULL ^ S][i]; if (!S0) break; } Mx[j][i] = s; }
  u128 vb = 0, vl = 0; int pos = 0, eb = 0, el = 0;
  for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { i64 mu = Mx[i][j] + Mx[j][i]; if (mu & 1) viol += 1000000; if ((mu >> 1) & 1) { vb |= (u128)1 << pos; eb ^= 1; } pos++; }
  for (int i = 0; i < m; i++) { i64 a = Bs[FULL][i] - Ee[FULL][i]; for (int j = 0; j < m; j++) if (j != i) a += Mx[j][i]; if (a & 1) viol += 1000000;
    if ((a >> 1) & 1) { vl |= (u128)1 << (pos + i); el ^= 1; } }
  if (!eb || !el) viol++;
  span_add(vb, &baseB, pivB, &dimB, !have); span_add(vb | vl, &baseLB, pivLB, &dimLB, !have); have = 1; cnt++;
}
int main(int argc, char **argv) { m = atoi(argv[1]); FULL = (1 << m) - 1; int np = m * (m - 1) / 2; u64 firstfull = 0;
  if (!strcmp(argv[2], "all")) { for (u64 code = 0; code < (1ULL << np); code++) one(code); }
  else { u64 N = strtoull(argv[3], 0, 10), s = strtoull(argv[4], 0, 10) * 0x9E3779B97F4A7C15ULL + 0x1234567ULL;
    for (u64 t = 0; t < N; t++) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; u64 x = s; s ^= s << 13; s ^= s >> 7; s ^= s << 17; one((x ^ (s << 1)) & ((np >= 64) ? ~0ULL : ((1ULL << np) - 1)));
      if (!firstfull && dimB == np - 1 && dimLB == np + m - 2) firstfull = t + 1; } }
  printf("m %d %s tournaments %llu violations %llu affine dim beta %d (C(m,2)-1 = %d) affine dim (ell,beta) %d (m+C(m,2)-2 = %d)",
         m, argv[2], cnt, viol, dimB, np - 1, dimLB, np + m - 2);
  if (firstfull) printf(" full after %llu tournaments", firstfull);
  printf("\n"); return 0; }
