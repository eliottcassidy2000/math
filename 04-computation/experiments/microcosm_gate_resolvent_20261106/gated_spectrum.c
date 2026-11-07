/* For every labelled tournament T on n <= 7 vertices: a0 = H(T), a1 = #orders with exactly one backward step,
   via DP over (subset, last vertex, #backward steps <= 1). Output the set of H(T) and of H(G(T)) = (n+2) a0 + a1. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(void) {
    for (int n = 1; n <= 7; n++) {
        int E = n * (n - 1) / 2; long long NT = 1LL << E;
        static char seenH[200000], seenG[2000000];
        memset(seenH, 0, sizeof(seenH)); memset(seenG, 0, sizeof(seenG));
        int maxH = 0, maxG = 0;
        static long long dp0[1 << 7][7], dp1[1 << 7][7];
        int A[7][7];
        for (long long code = 0; code < NT; code++) {
            int e = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { int b = (code >> e) & 1; A[i][j] = b; A[j][i] = !b; e++; }
            int full = (1 << n) - 1;
            for (int S = 0; S <= full; S++) for (int v = 0; v < n; v++) dp0[S][v] = dp1[S][v] = 0;
            for (int v = 0; v < n; v++) dp0[1 << v][v] = 1;
            for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) {
                if (!(S >> v & 1)) continue;
                long long x0 = dp0[S][v], x1 = dp1[S][v]; if (!x0 && !x1) continue;
                for (int w = 0; w < n; w++) { if (S >> w & 1) continue; int T = S | (1 << w);
                    if (A[v][w]) { dp0[T][w] += x0; dp1[T][w] += x1; } else { dp1[T][w] += x0; } }
            }
            long long a0 = 0, a1 = 0; for (int v = 0; v < n; v++) { a0 += dp0[full][v]; a1 += dp1[full][v]; }
            long long g = (n + 2) * a0 + a1;
            seenH[a0] = 1; if (a0 > maxH) maxH = a0;
            seenG[g] = 1; if (g > maxG) maxG = g;
        }
        printf("core order n = %d: H(T) values:", n);
        for (int h = 1; h <= maxH; h++) if (seenH[h]) printf(" %d", h);
        printf("\n    H(G(T)) values (gated strong order %d):", n + 2);
        int c = 0; for (int h = 1; h <= maxG; h++) if (seenG[h]) { if (c < 60) printf(" %d", h); c++; }
        printf("%s  [%d distinct, max %d]\n", c > 60 ? " ..." : "", c, maxG);
    }
    return 0;
}
