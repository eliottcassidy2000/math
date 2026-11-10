// x2_box.c -- AUDIT 2 (auditor's own code).  The deletion-pair side of (4.3) and the box-norm chain of Proposition 4.3.
//
// For a tournament W on N vertices and a set of k disjoint marked pairs let N_k = number of directed odd cycles meeting
// every marked pair in exactly one vertex, and B(N - k, k) = E (-1)^(N_k) over uniform labelled W.
//   route (i)  : N_k by explicit counting of the directed Hamiltonian cycles of W[S] over all odd S meeting every pair once;
//   route (ii) : (-1)^(N_k) = prod over the 2^k transversals D (one vertex of each pair) of s(W - D), s = +1 / -1 for
//                H = 1 / 3 mod 4.   (A cycle is avoided by 2^(#pairs it misses) transversals if it contains no pair, by none
//                otherwise; with H = 1 + 2 c_odd mod 4 this gives the parity of N_k.)   For k = 1 this is the middle term
//                of (4.3): E (-1)^((H(T-v) - H(T-v'))/2).
// usage:  x2_box N labelled                all 2^C(N,2) labelled tournaments, marked pairs {0,1},{2,3},...; both routes,
//                                          pointwise comparison of the two routes (N <= 7)
//         gentourng -q N [res/mod] | x2_box N classes C     class representatives with weights N!/|Aut|; route (ii) over ALL
//                                          placements of k disjoint pairs, k = 0..floor(N/2); route (i) is compared on every
//                                          placement if C = 1 (use for N <= 8), on a random 1/64 of the placements if C = 2,
//                                          on none if C = 0
// Output: raw integer sums; B(N-k,k) = sum_k / (matchings(N,k) * W) in class mode, sum_k / 2^C(N,2) in labelled mode.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64; typedef __int128 i128;
#define MAXN 10
static int N, FULL, adj[MAXN][MAXN], inm[MAXN], sc[MAXN], pm[MAXN], usedm;
static i64 Ee[1 << MAXN][MAXN], HS[1 << MAXN], Cp[1 << MAXN][MAXN], HC[1 << MAXN];
static signed char sg[1 << MAXN];
static void tables(int needcyc) {
  HS[0] = 1; sg[0] = 1;
  for (int S = 1; S <= FULL; S++) { i64 tot = 0;
    for (int vv = S; vv; vv &= vv - 1) { int v = __builtin_ctz(vv), rest = S ^ (1 << v); if (!rest) { Ee[S][v] = 1; tot += 1; continue; }
      i64 e = 0; for (int pp = rest & inm[v]; pp; pp &= pp - 1) e += Ee[rest][__builtin_ctz(pp)]; Ee[S][v] = e; tot += e; }
    HS[S] = tot; sg[S] = ((tot & 3) == 1) ? 1 : -1; if (!(tot & 1)) { fprintf(stderr, "even H\n"); exit(1); } }
  if (!needcyc) return;
  // Cp[S][u] = number of paths that start at the least vertex of S, cover S and end at u ; HC[S] = directed Hamiltonian cycles of W[S]
  for (int S = 1; S <= FULL; S++) { int s0 = __builtin_ctz(S); i64 hc = 0;
    for (int vv = S; vv; vv &= vv - 1) { int u = __builtin_ctz(vv);
      if (u == s0) { Cp[S][u] = (S == (1 << s0)); continue; }
      int rest = S ^ (1 << u); i64 e = 0; for (int pp = rest & inm[u]; pp; pp &= pp - 1) e += Cp[rest][__builtin_ctz(pp)]; Cp[S][u] = e;
      if (adj[u][s0]) hc += e; }
    HC[S] = (__builtin_popcount(S) >= 3) ? hc : 0; }
}
static i64 aut_rec(int k) { if (k == N) return 1; i64 c = 0;
  for (int w = 0; w < N; w++) { if ((usedm >> w & 1) || sc[w] != sc[k]) continue; int ok = 1;
    for (int j = 0; j < k && ok; j++) if (adj[j][k] != adj[pm[j]][w]) ok = 0;
    if (ok) { pm[k] = w; usedm |= 1 << w; c += aut_rec(k + 1); usedm &= ~(1 << w); } }
  return c; }
static int pa[MAXN], pb[MAXN];                 // current marked pairs
static int sign_route2(int k) {                // product over the 2^k transversals D of s(W - D) (Gray code over the pairs)
  int D = 0; for (int t = 0; t < k; t++) D |= 1 << pa[t];
  int s = sg[FULL ^ D];
  for (unsigned c = 1; c < (1u << k); c++) { int t = __builtin_ctz(c); D ^= (1 << pa[t]) | (1 << pb[t]); s *= sg[FULL ^ D]; }
  return s; }
static int sign_route1(int k) {                // explicit count of odd cycles meeting every marked pair exactly once
  i64 nk = 0; for (int S = 1; S <= FULL; S++) { if (!(__builtin_popcount(S) & 1) || !HC[S]) continue; int ok = 1;
    for (int t = 0; t < k && ok; t++) if (((S >> pa[t]) & 1) + ((S >> pb[t]) & 1) != 1) ok = 0; if (ok) nk += HC[S]; }
  return (nk & 1) ? -1 : 1; }
static u64 nplace[MAXN / 2 + 1], ncmp = 0, fcmp = 0; static int cmpmode = 0;
static u64 rng = 0x9E3779B97F4A7C15ULL; static u64 xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static i64 cur[MAXN / 2 + 1];
static void rec(int start, int used, int k) {  // all sets of disjoint pairs, each set once (pairs in increasing order of least element)
  { int s = sign_route2(k); cur[k] += s; nplace[k]++;
    if (cmpmode == 1 || (cmpmode == 2 && (xr() >> 58) == 0)) { ncmp++; if (sign_route1(k) != s) fcmp++; } }
  if (2 * (k + 1) > N) return;
  for (int a = start; a < N; a++) { if (used >> a & 1) continue;
    for (int b = a + 1; b < N; b++) { if (used >> b & 1) continue; pa[k] = a; pb[k] = b; rec(a + 1, used | (1 << a) | (1 << b), k + 1); } }
}
static void pr128(i128 x) { if (x < 0) { putchar('-'); x = -x; } char buf[64]; int k = 0; if (!x) { printf("0"); return; } while (x) { buf[k++] = (char)('0' + (int)(x % 10)); x /= 10; } while (k) putchar(buf[--k]); }

int main(int argc, char **argv) {
  if (argc < 3) return 2; N = atoi(argv[1]); FULL = (1 << N) - 1; if (N < 2 || N > MAXN) return 2;
  int kmax = N / 2;
  if (!strcmp(argv[2], "labelled")) { int np = N * (N - 1) / 2; i64 s1[MAXN / 2 + 1] = {0}, s2[MAXN / 2 + 1] = {0}; u64 mism = 0, cnt = 0;
    for (int t = 0; t < kmax; t++) { pa[t] = 2 * t; pb[t] = 2 * t + 1; }
    for (u64 code = 0; code < (1ULL << np); code++) { int b = 0;
      for (int i = 0; i < N; i++) { inm[i] = 0; adj[i][i] = 0; }
      for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { int x = (int)(code >> b & 1); b++; adj[i][j] = x; adj[j][i] = !x; if (x) inm[j] |= 1 << i; else inm[i] |= 1 << j; }
      tables(1); cnt++;
      for (int k = 0; k <= kmax; k++) { int a = sign_route1(k), c = sign_route2(k); s1[k] += a; s2[k] += c; if (a != c) mism++; } }
    printf("N %d labelled tournaments %llu pointwise mismatches route1/route2 %llu ;", N, cnt, mism);
    for (int k = 0; k <= kmax; k++) printf(" k %d n %d sum_route1 %lld sum_route2 %lld", k, N - k, s1[k], s2[k]);
    printf("\n"); return 0; }
  cmpmode = argc > 3 ? atoi(argv[3]) : 0; char line[256]; u64 nt = 0; i64 Wtot = 0, sD = 0, sDel = 0; i128 tot[MAXN / 2 + 1]; for (int k = 0; k <= kmax; k++) tot[k] = 0;
  i64 nf = 1; for (int i = 2; i <= N; i++) nf *= i;
  while (fgets(line, sizeof line, stdin)) { if (line[0] != '0' && line[0] != '1') continue; int b = 0;
    for (int i = 0; i < N; i++) { inm[i] = 0; sc[i] = 0; adj[i][i] = 0; }
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { int x = line[b++] == '1'; adj[i][j] = x; adj[j][i] = !x; if (x) { inm[j] |= 1 << i; sc[i]++; } else { inm[i] |= 1 << j; sc[j]++; } }
    tables(cmpmode != 0); usedm = 0; i64 au = aut_rec(0); if (au <= 0 || nf % au) { fprintf(stderr, "bad aut\n"); return 1; } i64 w = nf / au;
    for (int k = 0; k <= kmax; k++) cur[k] = 0; memset(nplace, 0, sizeof nplace);
    rec(0, 0, 0);
    for (int k = 0; k <= kmax; k++) tot[k] += (i128)w * cur[k];
    Wtot += w; sD += w * sg[FULL]; { i64 s = 0; for (int v = 0; v < N; v++) s += sg[FULL ^ (1 << v)]; sDel += w * s; } nt++; }
  printf("N %d classes %llu W %lld sumD %lld sumDel %lld route1_compared %llu route1_mismatch %llu ;", N, nt, Wtot, sD, sDel, ncmp, fcmp);
  for (int k = 0; k <= kmax; k++) { printf(" k %d n %d placements %llu sum ", k, N - k, nplace[k]); pr128(tot[k]); }
  printf("\n"); return 0;
}
