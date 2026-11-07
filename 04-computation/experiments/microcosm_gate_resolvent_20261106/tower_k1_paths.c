/* H (Hamiltonian paths) of the nested-gate tower Gamma^k(K1): vertices ordered b_k,...,b_1,v0,a_1,...,a_k (transitive
   spine) with the k symmetric chords a_j -> b_j reversed. Also m (orders with exactly one backward step) for the
   gate recursion H(Gamma T) = (n+2) H(T) + m(T). Bitmask DP. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(void) {
    for (int k = 0; k <= 8; k++) {
        int n = 2 * k + 1; int A[17][17];
        /* positions 0..n-1 along the spine: pos(b_j) = k - j, pos(v0) = k, pos(a_j) = k + j */
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] = (i < j);
        for (int j = 1; j <= k; j++) { int b = k - j, a = k + j; A[b][a] = 0; A[a][b] = 1; }
        long long *dp0 = calloc((size_t)(1 << n) * n, sizeof(long long)), *dp1 = calloc((size_t)(1 << n) * n, sizeof(long long));
        for (int v = 0; v < n; v++) dp0[(1 << v) * n + v] = 1;
        for (int S = 1; S < (1 << n); S++) for (int v = 0; v < n; v++) {
            if (!(S >> v & 1)) continue; long long x0 = dp0[S * n + v], x1 = dp1[S * n + v]; if (!x0 && !x1) continue;
            for (int w = 0; w < n; w++) { if (S >> w & 1) continue; int T = S | (1 << w);
                if (A[v][w]) { dp0[T * n + w] += x0; dp1[T * n + w] += x1; } else dp1[T * n + w] += x0; }
        }
        long long h = 0, m = 0; int full = (1 << n) - 1;
        for (int v = 0; v < n; v++) { h += dp0[full * n + v]; m += dp1[full * n + v]; }
        printf("k = %d (n = %2d): H = %lld, m = %lld, c3 = k^2 = %d\n", k, n, h, m, k * k);
        free(dp0); free(dp1);
    }
    return 0;
}
