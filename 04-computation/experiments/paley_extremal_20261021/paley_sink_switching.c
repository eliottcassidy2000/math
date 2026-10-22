/* paley_sink_switching.c -- HYP-9357 side test (spark from the wide-net survey, THM-4602): the switching class of
 * "Paley P_q plus a sink" on q+1 points (the Legendre chirotope class of P^1(F_q), on which PSL(2,q) acts 2-transitively).
 * For every subset X of the q+1 vertices containing a fixed vertex (2^q members up to complement) reverse all arcs between
 * X and its complement and count Hamiltonian paths by subset DP.  Reports the maximum, its multiplicity and the H-spectrum
 * of the class, to compare with the global maxima of OEIS A038375 (H_max(4) = 5, H_max(8) = 661, H_max(12) unknown).
 * mac-mini-2026-10-21-fanopsp.  Usage: ./paley_sink_switching q */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned long long u64;
static int n; static int adj[32][32];
static u64 hamcount(void){
  int N = 1 << n; u64 *dp = calloc((size_t)N * n, 8);
  for(int v = 0; v < n; v++) dp[((size_t)1 << v) * n + v] = 1;
  for(int S = 1; S < N; S++) for(int v = 0; v < n; v++){ u64 x = dp[(size_t)S * n + v]; if(!x) continue;
    for(int w = 0; w < n; w++) if(!(S >> w & 1) && adj[v][w]) dp[(size_t)(S | 1 << w) * n + w] += x; }
  u64 t = 0; for(int v = 0; v < n; v++) t += dp[(size_t)(N - 1) * n + v]; free(dp); return t; }
int main(int argc, char **argv){
  int q = atoi(argv[1]); n = q + 1; int base[32][32]; memset(base, 0, sizeof base);
  int qr[64] = {0}; for(int x = 1; x < q; x++) qr[x * x % q] = 1;
  for(int a = 0; a < q; a++) for(int b = 0; b < q; b++) if(a != b && qr[((b - a) % q + q) % q]) base[a][b] = 1;
  for(int a = 0; a < q; a++) base[a][q] = 1;               /* vertex q = infinity is a sink */
  u64 best = 0, cnt = 0; u64 hist_n = 0; u64 *vals = malloc(sizeof(u64) << q);
  for(int X = 0; X < (1 << q); X++){                      /* subsets of {0..q-1}; infinity never switched (complement) */
    for(int a = 0; a < n; a++) for(int b = 0; b < n; b++){ int sa = (a < q) ? (X >> a & 1) : 0, sb = (b < q) ? (X >> b & 1) : 0;
      adj[a][b] = (sa != sb) ? base[b][a] : base[a][b]; }
    u64 h = hamcount(); vals[hist_n++] = h; if(h > best){ best = h; cnt = 1; } else if(h == best) cnt++; }
  /* distinct values */
  for(u64 i = 1; i < hist_n; i++){ u64 v = vals[i]; long j = i - 1; while(j >= 0 && vals[j] > v){ vals[j + 1] = vals[j]; j--; } vals[j + 1] = v; }
  printf("q=%d (n=%d): switching class of P_q + sink, %llu members; max H = %llu (attained %llu times)\n  distinct H values:", q, n, hist_n, best, cnt);
  u64 prev = (u64)-1; int k = 0; for(u64 i = 0; i < hist_n; i++){ if(vals[i] != prev){ u64 c = 0; for(u64 j = i; j < hist_n && vals[j] == vals[i]; j++) c++; printf(" %llu(x%llu)", vals[i], c); prev = vals[i]; k++; } }
  printf("\n"); return 0; }
