/* Audit: the rows r = (0, s2) of the two-run observer M_2 of Paley P19 (connection set QR_19), i.e.
   M_2[(0,s2),(t1,t2)] = sum over U containing 0 and t1, with complement W containing s2 and t2, of
   (#Ham paths of T[U] from 0 to t1) * (#Ham paths of T[W] from s2 to t2).
   These rows determine every block B_j (representatives R = tuples with first coordinate 0).
   Output lines: s2 t1 t2 value (decimal; __int128 accumulation). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
typedef unsigned __int128 u128;
#ifndef n
#define n 19
#endif
static int adj[n][n];
static void pr128(u128 x) { char buf[64]; int i = 63; buf[i] = 0; if (x == 0) { putchar('0'); return; }
  while (x) { buf[--i] = '0' + (int)(x % 10); x /= 10; } fputs(buf + i, stdout); }
int main(void) {
  int qr[n] = {0};
  for (int x = 1; x < n; x++) qr[(x * x) % n] = 1;
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) adj[i][j] = (i != j) && qr[((j - i) % n + n) % n];
  const int B = n - 1;            /* 18 bits for vertices 1..18 */
  const u64 NM = 1ULL << B;
  /* dp0[m][v]: paths from 0 covering {0} u S(m) ending at v (S(m) = vertices j with bit j-1 of m) */
  u64 *dp0 = calloc(NM * n, sizeof(u64));
  dp0[0 * n + 0] = 1;
  for (u64 m = 0; m < NM; m++) for (int v = 0; v < n; v++) { u64 c = dp0[m * n + v]; if (!c) continue;
      for (int w = 1; w < n; w++) if (!((m >> (w - 1)) & 1) && adj[v][w]) dp0[(m | (1ULL << (w - 1))) * n + w] += c; }
  u64 *dq = malloc(NM * n * sizeof(u64));
  u128 *acc = malloc((size_t)n * n * sizeof(u128));
  for (int s2 = 1; s2 < n; s2++) {
    memset(dq, 0, NM * n * sizeof(u64));
    dq[(1ULL << (s2 - 1)) * n + s2] = 1;
    for (u64 m = 0; m < NM; m++) for (int v = 1; v < n; v++) { u64 c = dq[m * n + v]; if (!c) continue;
        for (int w = 1; w < n; w++) if (!((m >> (w - 1)) & 1) && adj[v][w]) dq[(m | (1ULL << (w - 1))) * n + w] += c; }
    memset(acc, 0, (size_t)n * n * sizeof(u128));
    for (u64 m = 0; m < NM; m++) {
      u64 W = (NM - 1) & ~m;
      if (!((W >> (s2 - 1)) & 1)) continue;
      for (int t1 = 0; t1 < n; t1++) { u64 a = dp0[m * n + t1]; if (!a) continue;
        for (int t2 = 1; t2 < n; t2++) { u64 b = dq[W * n + t2]; if (!b) continue;
          acc[t1 * n + t2] += (u128)a * b; } }
    }
    for (int t1 = 0; t1 < n; t1++) for (int t2 = 0; t2 < n; t2++) if (t1 != t2) {
      printf("%d %d %d ", s2, t1, t2); pr128(acc[t1 * n + t2]); putchar('\n'); }
    fprintf(stderr, "s2=%d done\n", s2);
  }
  return 0;
}
