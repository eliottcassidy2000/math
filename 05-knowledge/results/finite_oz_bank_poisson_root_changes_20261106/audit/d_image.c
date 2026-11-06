/* d_image.c - auditor's own exact O/Z tail-image engine (HYP-9480 (A)-(C) evidence).
 * Algorithm (different from the packet's pruned DFS): build Q[0..mmax] as an array from the defining recurrence
 * (Q(2k) = Q(k), Q(2k+1) = Q(k) + k + (k mod 2); cross-checked against the pair-sum definition in a_identities.c),
 * then scan EVERY m = 0..mmax linearly, H = h + R m - d Q[m] - e (m & 1), where mmax = floor((HI - 1 - h)/(R - d))
 * is the complete horizon from H >= h + (R - d) m (Session 205).  Values are marked in a per-seed bitmap of the odd
 * integers of the window [LO, HI).
 * usage: d_image SEEDFILE LO HI M [U]
 *   prints for every seed and every residue class c mod M that contains odd N of the window:
 *       S h R e d M c n_c p_c v_c        (n_c odd N in class, p_c parameters with multiplicity, v_c distinct values)
 *   with U: also the union of all seeds:   U M c n_c union_hits_c
 *   and pairwise intersection counts are NOT computed here (see d_pairs mode below).
 * usage: d_image SEEDFILE LO HI 0 P i j  -> |Im_i cap Im_j|, |Im_i|, |Im_j|, #odd in window (seed lines i, j 0-based)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll; typedef unsigned long long u64;
static uint32_t *Q; static ll qn = 0;
static void build_Q(ll mmax){
  if (mmax + 1 <= qn) return;
  free(Q); qn = mmax + 1; Q = malloc((size_t)qn * sizeof(uint32_t)); if (!Q){ fprintf(stderr, "oom Q\n"); exit(1); }
  Q[0] = 0; for (ll m = 1; m < qn; m++){ ll k = m >> 1; Q[m] = (m & 1) ? (uint32_t)(Q[k] + k + (k & 1)) : Q[k]; }
}
static ll LO, HI, NODD; static uint8_t *bm, *ubm;
static inline void setbit(uint8_t *b, ll i){ b[i >> 3] |= (uint8_t)(1u << (i & 7)); }
static inline int getbit(const uint8_t *b, ll i){ return (b[i >> 3] >> (i & 7)) & 1; }
static ll image(ll h, ll R, ll e, ll d, ll M, ll *p){   /* fills bm; returns #parameters in window */
  ll mmax = (HI - 1 - h) / (R - d); if (mmax < 0) mmax = -1; build_Q(mmax < 0 ? 0 : mmax);
  memset(bm, 0, (size_t)(NODD / 8 + 2)); ll np = 0;
  for (ll m = 0; m <= mmax; m++){
    ll H = h + R * m - d * (ll)Q[m] - e * (m & 1);
    if (H >= LO && H < HI){ if (!(H & 1)){ fprintf(stderr, "even value!\n"); exit(3); } setbit(bm, (H - LO) >> 1); np++; if (M > 0) p[H % M]++; }
  }
  return np;
}
int main(int argc, char **argv){
  if (argc < 5){ fprintf(stderr, "usage\n"); return 1; }
  FILE *fp = fopen(argv[1], "r"); LO = atoll(argv[2]); HI = atoll(argv[3]); ll M = atoll(argv[4]);
  if (LO & 1){ fprintf(stderr, "LO must be even\n"); return 1; }
  NODD = (HI - LO) / 2; bm = calloc((size_t)(NODD / 8 + 2), 1); ubm = calloc((size_t)(NODD / 8 + 2), 1);
  ll hs[512], Rs[512], es[512], ds[512]; int ns = 0;
  while (ns < 512 && fscanf(fp, "%lld %lld %lld %lld", &hs[ns], &Rs[ns], &es[ns], &ds[ns]) == 4) ns++;
  if (M == 0){ /* pair mode */
    int i = atoi(argv[6]), j = atoi(argv[7]);
    image(hs[i], Rs[i], es[i], ds[i], 0, NULL); memcpy(ubm, bm, (size_t)(NODD / 8 + 2));
    image(hs[j], Rs[j], es[j], ds[j], 0, NULL);
    ll a = 0, b = 0, ab = 0; for (ll x = 0; x < NODD; x++){ int u = getbit(ubm, x), v = getbit(bm, x); a += u; b += v; ab += u & v; }
    printf("P %lld %lld %lld %lld  %lld %lld %lld %lld  inter %lld  A %lld  B %lld  nodd %lld\n", hs[i], Rs[i], es[i], ds[i], hs[j], Rs[j], es[j], ds[j], ab, a, b, NODD);
    return 0; }
  int want_union = argc > 5 && argv[5][0] == 'U';
  ll *n = calloc(M, sizeof(ll)), *p = calloc(M, sizeof(ll)), *v = calloc(M, sizeof(ll)), *uh = calloc(M, sizeof(ll));
  for (ll x = LO | 1; x < HI; x += 2) n[x % M]++;
  for (int s = 0; s < ns; s++){
    memset(p, 0, M * sizeof(ll)); memset(v, 0, M * sizeof(ll));
    image(hs[s], Rs[s], es[s], ds[s], M, p);
    for (ll x = 0; x < NODD; x++) if (getbit(bm, x)){ v[(LO + 1 + 2 * x) % M]++; }
    if (want_union) for (ll k = 0; k < NODD / 8 + 2; k++) ubm[k] |= bm[k];
    for (ll c = 0; c < M; c++) if (n[c]) printf("S %lld %lld %lld %lld %lld %lld %lld %lld %lld\n", hs[s], Rs[s], es[s], ds[s], M, c, n[c], p[c], v[c]);
  }
  if (want_union){
    for (ll x = 0; x < NODD; x++) if (getbit(ubm, x)) uh[(LO + 1 + 2 * x) % M]++;
    for (ll c = 0; c < M; c++) if (n[c]) printf("U %lld %lld %lld %lld\n", M, c, n[c], uh[c]);
  }
  return 0;
}
