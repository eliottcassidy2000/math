/* H (Hamiltonian paths) and hc (directed Hamiltonian cycles) for tournaments up to 20 vertices.
   input lines: n then n rows of 0/1 (adjacency, row i col j = 1 iff i -> j) */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64;
int main(void) {
    int n;
    while (scanf("%d", &n) == 1) {
        unsigned adj[32] = {0};
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { int b; scanf("%d", &b); if (b) adj[i] |= 1u << j; }
        size_t full = ((size_t)1 << n) - 1;
        u64 *dp = calloc((full + 1) * n, sizeof(u64));
        for (int v = 0; v < n; v++) dp[((size_t)1 << v) * n + v] = 1;
        for (size_t S = 1; S <= full; S++) for (int v = 0; v < n; v++) {
            u64 c = dp[S * n + v]; if (!c) continue;
            unsigned m = adj[v] & ~(unsigned)S;
            while (m) { int w = __builtin_ctz(m); m &= m - 1; dp[(S | ((size_t)1 << w)) * n + w] += c; }
        }
        u64 H = 0; for (int v = 0; v < n; v++) H += dp[full * n + v];
        for (size_t z = 0; z < (full + 1) * n; z++) dp[z] = 0;
        dp[1 * n + 0] = 1;
        for (size_t S = 1; S <= full; S += 2) for (int v = 0; v < n; v++) {
            u64 c = dp[S * n + v]; if (!c) continue;
            unsigned m = adj[v] & ~(unsigned)S;
            while (m) { int w = __builtin_ctz(m); m &= m - 1; dp[(S | ((size_t)1 << w)) * n + w] += c; }
        }
        u64 hc = 0; for (int v = 1; v < n; v++) if (adj[v] & 1u) hc += dp[full * n + v];
        printf("%d %llu %llu\n", n, H, hc);
        free(dp);
        fflush(stdout);
    }
    return 0;
}
