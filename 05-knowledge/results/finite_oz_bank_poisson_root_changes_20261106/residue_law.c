/* residue_law.c - kind-pasteur-2026-11-06-S13 (HYP-9480 (A), prediction P6: residue-resolved Poisson law).
 * For each seed "h R e d" (= D U E d) of SEEDFILE, the exact O/Z tail image H_s(m) = h + R m - d Q(m) - e (m mod 2)
 * (Session 205) inside the window [LO, HI) is enumerated by a DFS over m (Q updated exactly, pruned by
 * H_s(m) >= h + (R - d) m). For every residue class c of the odd N in the window modulo M it prints
 *     h R e d M c n_c p_c v_c
 * n_c = number of odd N in the window with N = c (mod M), p_c = number of parameters m >= 0 with H_s(m) in the window
 * and H_s(m) = c (mod M) (with multiplicity), v_c = number of distinct such values.
 * usage: residue_law SEEDFILE LO HI M      (LO even, M odd)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll; typedef unsigned char u8;
static ll LO, HI; static int M; static u8 *bits; static ll *par; static ll H0, R0, E0, D0;
static inline void mark(ll v){
  if (v >= LO && v < HI && (v & 1)){ ll i = (v - LO) >> 1; bits[i >> 3] |= (u8)(1u << (i & 7)); par[v % M]++; }
}
static void dfs(ll m, ll q){
  if (H0 + (R0 - D0) * m >= HI) return;
  mark(H0 + R0 * m - D0 * q - E0 * (m & 1));
  dfs(2 * m, q);
  dfs(2 * m + 1, q + m + (m & 1));
}
int main(int argc, char **argv){
  if (argc < 5){ fprintf(stderr, "usage: residue_law SEEDFILE LO HI M\n"); return 1; }
  FILE *f = fopen(argv[1], "r"); LO = atoll(argv[2]); HI = atoll(argv[3]); M = atoi(argv[4]);
  if ((LO & 1) || (M % 2) == 0){ fprintf(stderr, "LO must be even and M odd\n"); return 1; }
  ll nodd = (HI - LO) / 2, nb = nodd / 8 + 2;
  bits = calloc(nb, 1); par = calloc(M, sizeof(ll));
  ll *n = calloc(M, sizeof(ll)), *v = calloc(M, sizeof(ll));
  for (ll x = LO | 1; x < HI; x += 2) n[x % M]++;
  ll h, R, e, d;
  while (fscanf(f, "%lld %lld %lld %lld", &h, &R, &e, &d) == 4){
    memset(bits, 0, nb); memset(par, 0, M * sizeof(ll)); memset(v, 0, M * sizeof(ll));
    H0 = h; R0 = R; E0 = e; D0 = d;
    mark(h); dfs(1, 0);
    for (ll x = LO | 1; x < HI; x += 2){ ll i = (x - LO) >> 1; if ((bits[i >> 3] >> (i & 7)) & 1) v[x % M]++; }
    for (int c = 0; c < M; c++) printf("%lld %lld %lld %lld %d %d %lld %lld %lld\n", h, R, e, d, M, c, n[c], par[c], v[c]);
  }
  return 0;
}
