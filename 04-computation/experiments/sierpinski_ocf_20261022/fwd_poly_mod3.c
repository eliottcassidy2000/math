/* fwd_poly_mod3.c -- check of the completed THM-086: for every tournament T on n vertices,
 *    F(T,x) = sum_sigma x^(#forward steps) = H (x-1)^(n-1)            (mod 3)  if n is odd,
 *                                           = H (x-1)^(n-2) (x+1)      (mod 3)  if n is even,
 * and the exact identity  sum_(X subset V) H(T[X]) H(T[V\X]) = sum_R 4^|R| 2^(n-|V(R)|)  is used in the proof
 * (R over sets of vertex-disjoint odd cycles); here we check the polynomial statement on gentourng streams.
 * F is computed exactly by a DP over (mask, last vertex, forward count).  mac-mini-2026-10-22-sierp. */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64;
static u64 dp[1 << 10][10][10];
int main(void){
  char line[128]; long tested = 0, bad = 0; int n = 0;
  while (fgets(line, sizeof line, stdin)){
    if (line[0] == '>') continue;
    int L = (int)strlen(line); while (L && (line[L-1] == '\n')) line[--L] = 0;
    n = 0; while (n * (n - 1) / 2 < L) n++;
    uint32_t out[16] = {0}; int p = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, p++){ if (line[p] == '1') out[i] |= 1u << j; else out[j] |= 1u << i; }
    uint32_t N = 1u << n;
    memset(dp, 0, sizeof(u64) * N * 10 * 10);
    for (int v = 0; v < n; v++) dp[1u << v][v][0] = 1;
    for (uint32_t m = 1; m < N; m++) for (int v = 0; v < n; v++) if ((m >> v) & 1) for (int f = 0; f < n; f++){
      u64 c = dp[m][v][f]; if (!c) continue;
      for (int w = 0; w < n; w++) if (!((m >> w) & 1)){ int fw = (out[v] >> w) & 1; dp[m | (1u << w)][w][f + fw] += c; } }
    u64 F[16] = {0};
    for (int v = 0; v < n; v++) for (int f = 0; f < n; f++) F[f] += dp[N - 1][v][f];
    u64 H = F[n - 1];
    /* target polynomial mod 3 */
    int poly[16] = {0}; int deg = (n % 2) ? n - 1 : n - 2;
    /* (x-1)^deg mod 3 via binomials */
    int row[16] = {0}; row[0] = 1;
    for (int d = 0; d < deg; d++){ for (int k = d + 1; k >= 1; k--) row[k] = ((row[k - 1] - row[k]) % 3 + 3) % 3; row[0] = (3 - row[0]) % 3; }
    /* row now holds coefficients of (x-1)^deg  (constructed as (x-1)*previous) */
    if (n % 2){ for (int k = 0; k <= deg; k++) poly[k] = row[k]; }
    else { for (int k = 0; k <= deg; k++){ poly[k] = (poly[k] + row[k]) % 3; poly[k + 1] = (poly[k + 1] + row[k]) % 3; } }
    int h3 = (int)(H % 3);
    for (int k = 0; k < n; k++) if ((int)(F[k] % 3) != (poly[k] * h3) % 3) { bad++; break; }
    tested++;
  }
  printf("n=%d: %ld tournaments, %ld violate F = H*(x-1)^(n-1) [odd n] or H*(x-1)^(n-2)(x+1) [even n] mod 3\n", n, tested, bad);
  return 0;
}
