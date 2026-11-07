/* audit_H1: brute-force Hamiltonian path (and, with -c, Hamiltonian cycle mod rotation) counts of tournaments read as
   "n upper-triangle-string" lines from stdin (bit (i,j), i<j, '1' means i->j); subset DP, 64-bit counts, n <= 21. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

int main(int argc, char **argv) {
    int cycles = (argc > 1 && strcmp(argv[1], "-c") == 0);
    int n;
    static char buf[1 << 16];
    while (scanf("%d %65535s", &n, buf) == 2) {
        uint32_t out[32];
        memset(out, 0, sizeof out);
        int k = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) {
                if (buf[k] == '1') out[i] |= 1u << j; else out[j] |= 1u << i;
                k++;
            }
        if (n == 1) { printf("%d\n", cycles ? 0 : 1); continue; }
        size_t size = (size_t)1 << n;
        uint64_t *dp = calloc(size * n, sizeof(uint64_t));
        if (!dp) { fprintf(stderr, "oom\n"); return 1; }
        if (cycles) dp[1 * n + 0] = 1;
        else for (int v = 0; v < n; v++) dp[((size_t)1 << v) * n + v] = 1;
        for (size_t S = 1; S < size; S++) {
            if (cycles && !(S & 1)) continue;
            for (int v = 0; v < n; v++) {
                uint64_t c = dp[S * n + v];
                if (!c) continue;
                uint32_t nxt = out[v] & ~(uint32_t)S;
                while (nxt) {
                    int w = __builtin_ctz(nxt);
                    nxt &= nxt - 1;
                    dp[(S | ((size_t)1 << w)) * n + w] += c;
                }
            }
        }
        uint64_t tot = 0;
        size_t full = size - 1;
        for (int v = 0; v < n; v++) {
            if (cycles) { if (v != 0 && (out[v] & 1u)) tot += dp[full * n + v]; }
            else tot += dp[full * n + v];
        }
        printf("%llu\n", (unsigned long long)tot);
        free(dp);
    }
    return 0;
}
