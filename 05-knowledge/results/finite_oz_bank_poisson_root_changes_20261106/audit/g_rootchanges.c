/* g_rootchanges.c - auditor's independent enumerator of F_k (HYP-9480 (D), P4), k <= 2.
 * NO closed-form tail: every word t . w_1 (J w_2) ... (J w_(k+1)) is walked letter by letter with the physical moves
 *     O(D,E,U,d) = (D,0,2U,d),  Z = (D-E+U,d,2U-d,d),  J = (D-E+U,D-E,D+U,D)   (verified literally in c_literal.c),
 * from t = (1,0,2,1).  H of a word is its D.  Canonical words: trailing O letters deleted (they do not change D), so a
 * word is recorded at the node of its last letter if that letter is Z or J, and the empty word (t itself) at the root.
 * Pruning (complete): D never decreases (O keeps D, Z/J raise it); after an O the next Z or J gives D + 2U, and more O's
 * only raise that, so an O-child is explored only if D + 2U <= XMAX.
 * Output: per exact root-change count j: masses #(canonical words with D <= X) at several X; hit bitmaps to XMAX;
 * odd misses of F_(<=j) up to XA and XMAX; misses above 10^5; per dyadic window counts.
 * usage: g_rootchanges XMAX XA K
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
static ll XMAX; static int K; static uint8_t *hit[3];
static const ll THR[] = {100000, 1000000, 4000000, 10000000, 100000000}; static const int NT = 5;
static ll mass[3][5];
static inline void rec(ll D, int j){
  if (D > XMAX) return;
  hit[j][D >> 4] |= (uint8_t)(1u << ((D >> 1) & 7));
  for (int i = 0; i < NT; i++) if (D <= THR[i]) mass[j][i]++;
}
static void dfs(ll D, ll E, ll U, ll d, int j, int lastO){
  if (!lastO) rec(D, j);
  ll D2 = D - E + U;
  if (D2 <= XMAX){
    dfs(D2, d, 2 * U - d, d, j, 0);                 /* Z */
    if (j < K) dfs(D2, D - E, D + U, D, j + 1, 0);  /* J */
  }
  if (D + 2 * U <= XMAX) dfs(D, 0, 2 * U, d, j, 1); /* O */
}
static inline int get(const uint8_t *b, ll v){ return (b[v >> 4] >> ((v >> 1) & 7)) & 1; }
int main(int argc, char **argv){
  XMAX = atoll(argv[1]); ll XA = atoll(argv[2]); K = atoi(argv[3]);
  for (int i = 0; i <= K; i++) hit[i] = calloc((size_t)(XMAX >> 4) + 2, 1);
  dfs(1, 0, 2, 1, 0, 0);
  for (int j = 0; j <= K; j++){
    printf("exactly %d root changes, masses:", j);
    for (int i = 0; i < NT; i++) if (THR[i] <= XMAX) printf(" N(%lld)=%lld", THR[i], mass[j][i]);
    printf("\n");
  }
  /* cumulative unions */
  for (int i = 1; i <= K; i++) for (ll b = 0; b < (XMAX >> 4) + 2; b++) hit[i][b] |= hit[i - 1][b];
  for (int i = 0; i <= K; i++){
    ll mA = 0, lastA = 0, mX = 0, lastX = 0;
    for (ll v = 1; v <= XMAX; v += 2) if (!get(hit[i], v)){ mX++; lastX = v; if (v <= XA){ mA++; lastA = v; } }
    printf("<=%d root changes: odd misses <= %lld: %lld (largest %lld); <= %lld: %lld (largest %lld)\n", i, XA, mA, lastA, XMAX, mX, lastX);
  }
  printf("misses of F_(<=%d) above 10^5 up to %lld:", K, XMAX); int n1 = 0, nn = 0;
  for (ll v = 100001; v <= XMAX; v += 2) if (!get(hit[K], v)){ printf(" %lld", v); nn++; if (v % 3 == 1) n1++; }
  printf("\n  count %d, of which = 1 (mod 3): %d\n", nn, n1);
  for (ll lo = 1; lo < XMAX; lo *= 2){ ll hi = lo * 2 > XMAX ? XMAX + 1 : lo * 2; printf("  [%lld,%lld):", lo, hi);
    for (int i = 0; i <= K; i++){ ll c = 0; for (ll v = lo | 1; v < hi; v += 2) if (!get(hit[i], v)) c++; printf(" k%d:%lld", i, c); } printf("\n"); }
  return 0;
}
