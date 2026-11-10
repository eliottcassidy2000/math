/* a_paths.c (S20 audit, independent of the packet's resp.c)
   Input : one tournament per line: n, then n*(n-1)/2 characters '0'/'1' (no spaces) listing the pairs (i<j)
           in the order (0,1),(0,2),...,(0,n-1),(1,2),...; '1' means i->j.
   Output: one line per tournament:  H Q | e_0 .. e_{n-1} | R_0 .. R_{n-1} | s_0 .. s_{n-1}
           H = # Hamiltonian paths (H(empty)=1), Q = sum_S H(T[S]) H(T[V-S]),
           e_v = # Hamiltonian paths ending at v, s_v = # starting at v,
           R_v = sum over paths P starting at v (P=(v) included) of H(T - V(P)).
   Method: DP over (subset, endpoint) for end-counts; start-counts by a second DP; R_v summed over subsets. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
typedef unsigned __int128 u128;
static void pr(u128 x){ char b[64]; int i = 63; b[i] = 0; if (!x) { putchar('0'); return; }
  while (x) { b[--i] = '0' + (int)(x % 10); x /= 10; } fputs(b + i, stdout); }
int main(void){
  int n; static char bits[256];
  static int adj[20][20];
  while (scanf("%d", &n) == 1) {
    if (n > 18) { fprintf(stderr, "n too big\n"); return 1; }
    if (n >= 2) { if (scanf("%255s", bits) != 1) return 1; } else bits[0] = 0;
    int k = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) adj[i][j] = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
      if (bits[k] == '1') adj[i][j] = 1; else if (bits[k] == '0') adj[j][i] = 1; else { fprintf(stderr, "bad bits\n"); return 1; }
      k++; }
    size_t N = (size_t)1 << n;
    u64 *f = calloc(N * (n ? n : 1), sizeof(u64));   /* f[S*n+v]: paths covering S ending at v */
    u64 *g = calloc(N * (n ? n : 1), sizeof(u64));   /* g[S*n+v]: paths covering S starting at v */
    u64 *H = calloc(N, sizeof(u64));
    H[0] = 1;
    for (size_t S = 1; S < N; S++) {
      u64 tot = 0;
      for (int v = 0; v < n; v++) {
        if (!((S >> v) & 1)) continue;
        size_t P = S & ~((size_t)1 << v);
        u64 a = 0, b = 0;
        if (P == 0) { a = 1; b = 1; }
        else {
          for (int w = 0; w < n; w++) {
            if (!((P >> w) & 1)) continue;
            if (adj[w][v]) a += f[P * n + w];      /* ... w -> v (v last) */
            if (adj[v][w]) b += g[P * n + w];      /* v -> w ... (v first) */
          }
        }
        f[S * n + v] = a; g[S * n + v] = b; tot += a;
      }
      H[S] = tot;
    }
    size_t F = N - 1;
    u128 Q = 0; for (size_t S = 0; S < N; S++) Q += (u128)H[S] * H[F & ~S];
    pr(H[F]); putchar(' '); pr(Q); printf(" |");
    for (int v = 0; v < n; v++) { putchar(' '); pr(f[F * n + v]); }
    printf(" |");
    for (int v = 0; v < n; v++) {
      u128 R = 0;
      for (size_t S = 0; S < N; S++) if ((S >> v) & 1) R += (u128)g[S * n + v] * H[F & ~S];
      putchar(' '); pr(R);
    }
    printf(" |");
    for (int v = 0; v < n; v++) { putchar(' '); pr(g[F * n + v]); }
    putchar('\n');
    free(f); free(g); free(H);
  }
  return 0;
}
