/* hk_big.c -- Held-Karp Hamiltonian path count (unsigned __int128 accumulator) for gentourng-style upper-triangle
 * codes on stdin, n <= 24 (memory 2^n * n * 16 bytes: 1.4 GB at n = 22, 3 GB at n = 23 -- use n <= 22 on 8 GB).
 * mac-mini-2026-10-22-sierp. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static void p128(u128 x) { char b[64]; int i = 63; b[i] = 0; if (!x) { putchar('0'); return; } while (x) { b[--i] = '0' + (int)(x % 10); x /= 10; } fputs(b + i, stdout); }
int main(void) {
    static char line[4096];
    while (fgets(line, sizeof line, stdin)) {
        int L = (int)strlen(line); while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = 0;
        if (!L || line[0] == '>') continue;
        int n = 0; while (n * (n - 1) / 2 < L) n++;
        uint32_t out[32] = {0}; int p = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, p++) { if (line[p] == '1') out[i] |= 1u << j; else out[j] |= 1u << i; }
        size_t N = (size_t)1 << n;
        u128 *dp = calloc(N * n, sizeof(u128));
        if (!dp) { fprintf(stderr, "alloc failed\n"); return 1; }
        for (int v = 0; v < n; v++) dp[((size_t)1 << v) * n + v] = 1;
        for (size_t m = 1; m < N; m++) for (int v = 0; v < n; v++) { u128 c = dp[m * n + v]; if (!c) continue;
            for (uint32_t q = out[v] & ~(uint32_t)m; q; q &= q - 1) { int w = __builtin_ctz(q); dp[(m | ((size_t)1 << w)) * n + w] += c; } }
        u128 H = 0; for (int v = 0; v < n; v++) H += dp[(N - 1) * n + v];
        p128(H); putchar('\n'); fflush(stdout); free(dp);
    }
    return 0;
}
