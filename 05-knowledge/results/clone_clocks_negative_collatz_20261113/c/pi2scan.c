/* pi2scan.c -- for single-module ternary clocks (S,u) read as "n u row0..row(n-1)", compute for j = 1..J
   H(S[u<-TT_j]) and Pi_2(S[u<-TT_j]) (ordered pairs of possibly empty paths covering the vertex set).
   Prints both sequences; the driver tests whether Pi_2 is ternary-shaped (rho 3^j + tau).                  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
int main(int argc, char **argv) {
    int J = argc > 1 ? atoi(argv[1]) : 4;
    int n0, u; long long r0[24];
    static ll *f = 0; size_t fcap = 0;
    while (scanf("%d %d", &n0, &u) == 2) {
        for (int i = 0; i < n0; i++) scanf("%lld", &r0[i]);
        printf("%d %d", n0, u);
        for (int j = 1; j <= J; j++) {
            int n = n0 - 1 + j; int adj[32][32]; memset(adj, 0, sizeof adj);
            int map[32], m = 0; for (int x = 0; x < n0; x++) if (x != u) map[x] = m++;
            int blk0 = m;
            for (int x = 0; x < n0; x++) { if (x == u) continue;
                for (int y = 0; y < n0; y++) { if (y == u || y == x) continue; if (r0[x] >> y & 1) adj[map[x]][map[y]] = 1; }
                for (int b = 0; b < j; b++) { if (r0[x] >> u & 1) adj[map[x]][blk0 + b] = 1; else adj[blk0 + b][map[x]] = 1; } }
            for (int a = 0; a < j; a++) for (int b = a + 1; b < j; b++) adj[blk0 + a][blk0 + b] = 1;
            size_t S = (size_t)1 << n; size_t need = S * n;
            if (need > fcap) { free(f); fcap = need; f = malloc(fcap * sizeof(ll)); }
            memset(f, 0, need * sizeof(ll));
            for (int v = 0; v < n; v++) f[((size_t)1 << v) * n + v] = 1;
            ll *h = calloc(S, sizeof(ll)); h[0] = 1;
            for (size_t mask = 1; mask < S; mask++) { ll tot = 0;
                for (int v = 0; v < n; v++) { ll c = f[mask * n + v]; if (!c) continue; tot += c;
                    for (int w = 0; w < n; w++) if (!(mask >> w & 1) && adj[v][w]) f[(mask | ((size_t)1 << w)) * n + w] += c; }
                h[mask] = tot; }
            ll H = h[S - 1], Pi2 = 0; for (size_t mask = 0; mask < S; mask++) Pi2 += h[mask] * h[(S - 1) ^ mask];
            printf(" %lld:%lld", H, Pi2); free(h);
        }
        printf("\n"); fflush(stdout);
    }
    return 0;
}
