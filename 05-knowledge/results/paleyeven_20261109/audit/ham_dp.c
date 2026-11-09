/* ham_dp.c - audit P1 (independent): Hamiltonian path counts by a forward subset DP (n <= 20).
 * dp[S][v] = number of directed paths with vertex set S ending at v.
 * H(T) = sum_v dp[V][v];  H(T - u) = sum_v dp[V \ u][v].
 * Input lines: name n out_0 ... out_{n-1}.  Output: name n H and every deletion count.
 * Compile: cc -O2 -o ham_dp ham_dp.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;

int main(void) {
    char name[64]; int n;
    while (scanf("%63s %d", name, &n) == 2) {
        uint32_t out[32];
        for (int i = 0; i < n; i++) { unsigned long long m; if (scanf("%llu", &m) != 1) return 1; out[i] = (uint32_t)m; }
        if (n > 20) { printf("%s n=%d skipped (n > 20)\n", name, n); continue; }
        size_t NS = (size_t)1 << n;
        u64 *dp = calloc(NS * n, sizeof(u64));
        if (!dp) { fprintf(stderr, "alloc\n"); return 1; }
        for (int v = 0; v < n; v++) dp[((size_t)1 << v) * n + v] = 1;
        for (size_t S = 1; S < NS; S++) {
            u64 *row = dp + S * n;
            for (int v = 0; v < n; v++) {
                u64 c = row[v]; if (!c) continue;
                uint32_t nx = out[v] & ~(uint32_t)S;
                while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dp[(S | ((size_t)1 << w)) * n + w] += c; }
            }
        }
        size_t F = NS - 1;
        u64 H = 0; for (int v = 0; v < n; v++) H += dp[F * n + v];
        printf("%s n=%d H=%llu del:", name, n, H);
        for (int u = 0; u < n; u++) {
            size_t S = F & ~((size_t)1 << u); u64 h = 0;
            for (int v = 0; v < n; v++) h += dp[S * n + v];
            printf(" %llu", h);
        }
        printf("\n"); fflush(stdout);
        free(dp);
    }
    return 0;
}
