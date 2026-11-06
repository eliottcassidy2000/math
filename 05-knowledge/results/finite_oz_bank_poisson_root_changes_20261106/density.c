/* density.c - kind-pasteur-2026-11-06-S13 (THM-4888, HYP-9480). Exact O/Z tail images of physical seeds.
 * Tail formula (Session 205): H_s(m) = h + R m - d Q(m) - e (m mod 2), m >= 0; Q(0)=0, Q(2m)=Q(m), Q(2m+1)=Q(m)+m+(m mod 2).
 * Parameters are enumerated by a depth-first search m -> 2m, 2m+1 (Q updated exactly), pruned by H_s(m) >= h + (R-d) m.
 * usage: density SEEDFILE X MODE LO HI
 *   MODE 1: for each seed "h R e d" print the number of distinct odd values in [LO,HI) and the number of odd integers there.
 *   MODE 2: union of all seeds; print the number of missed odd integers in [LO,HI), also per class mod 3.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll; typedef unsigned char u8;
static ll X; static u8 *bits; static ll H0, R0, E0, D0;
static inline void mark(ll v){ if (v >= 1 && v <= X && (v & 1)) bits[v >> 4] |= (u8)(1u << ((v >> 1) & 7)); }
static inline int get(ll v){ return (bits[v >> 4] >> ((v >> 1) & 7)) & 1; }
static void dfs(ll m, ll q){
  if (H0 + (R0 - D0) * m > X) return;
  mark(H0 + R0 * m - D0 * q - E0 * (m & 1));
  dfs(2 * m, q);
  dfs(2 * m + 1, q + m + (m & 1));
}
static void image(ll h, ll R, ll e, ll d){ H0 = h; R0 = R; E0 = e; D0 = d; mark(h); dfs(1, 0); }
int main(int argc, char **argv){
  if (argc < 6){ fprintf(stderr, "usage: density SEEDFILE X MODE LO HI\n"); return 1; }
  FILE *f = fopen(argv[1], "r"); X = atoll(argv[2]); int mode = atoi(argv[3]); ll lo = atoll(argv[4]), hi = atoll(argv[5]);
  ll nb = (X >> 4) + 2; bits = calloc(nb, 1); ll h, R, e, d;
  ll nodd = 0; for (ll v = lo | 1; v < hi; v += 2) nodd++;
  while (fscanf(f, "%lld %lld %lld %lld", &h, &R, &e, &d) == 4){
    if (mode == 1) memset(bits, 0, nb);
    image(h, R, e, d);
    if (mode == 1){ ll c = 0; for (ll v = lo | 1; v < hi; v += 2) c += get(v); printf("%lld %lld %lld %lld %lld %lld\n", h, R, e, d, c, nodd); }
  }
  if (mode == 2){ ll c = 0, c3[3] = {0, 0, 0}, n3[3] = {0, 0, 0};
    for (ll v = lo | 1; v < hi; v += 2){ n3[v % 3]++; if (!get(v)){ c++; c3[v % 3]++; } }
    printf("%lld %lld %lld %lld %lld %lld %lld %lld\n", c, nodd, c3[0], n3[0], c3[1], n3[1], c3[2], n3[2]); }
  return 0;
}
