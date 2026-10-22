/* forest_numbers.c -- the forward-step polynomial F(T, x) = sum over orderings sigma of x^(#forward steps) of a
 * tournament (upper-triangle code on stdin), evaluated modulo a prime p < 2^31 at x = 0, 1, ..., n - 1, by a
 * Held-Karp dynamic program over (set, last vertex):
 *     dp[m + w][w] = T[m] + (x - 1) * sum_(v in m, v -> w) dp[m][v],     T[m] = sum_(v in m) dp[m][v].
 * Interpolation and CRT (drt23_forest_exact.py) give F exactly, and THM-4757 section 5,
 *     sum_L (n - L)! N_L t^L = (1 - t)^(n - 1) F(T, (1 + t)/(1 - t)),
 * gives the linear-forest numbers N_L.  Checks: F(T, 0) = H(T), F(T, 1) = n!.
 * Memory 2^n * n * 4 + 2^n * 4 bytes (0.8 GB at n = 23).  Usage: forest_numbers p < code.  mac-mini-2026-10-22-loop4. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: forest_numbers p < code\n"); return 1; }
    uint64_t p = strtoull(argv[1], 0, 10);
    if (p >= (1ULL << 31)) { fprintf(stderr, "p must be < 2^31\n"); return 1; }
    static char line[4096];
    if (!fgets(line, sizeof line, stdin)) return 1;
    int L = (int)strlen(line); while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = 0;
    int n = 0; while (n * (n - 1) / 2 < L) n++;
    uint32_t out[32] = {0}; int q = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, q++) { if (line[q] == '1') out[i] |= 1u << j; else out[j] |= 1u << i; }
    size_t N = (size_t)1 << n;
    uint32_t *dp = malloc(N * n * sizeof(uint32_t)), *T = malloc(N * sizeof(uint32_t));
    if (!dp || !T) { fprintf(stderr, "alloc failed\n"); return 1; }
    for (int x = 0; x < n; x++) {
        uint64_t xm1 = (x + p - 1) % p;
        memset(dp, 0, N * n * sizeof(uint32_t));
        for (int v = 0; v < n; v++) dp[((size_t)1 << v) * n + v] = 1;
        for (size_t m = 1; m < N; m++) {
            uint32_t *row = dp + m * n;
            uint64_t tot = 0;
            if (m & (m - 1)) {                                   /* |m| >= 2: finalise */
                for (uint32_t r = (uint32_t)m; r; r &= r - 1) {
                    int w = __builtin_ctz(r);
                    uint64_t val = (T[m ^ ((size_t)1 << w)] + xm1 * row[w]) % p;
                    row[w] = (uint32_t)val; tot += val;
                }
            } else tot = 1;
            T[m] = (uint32_t)(tot % p);
            for (uint32_t r = (uint32_t)m; r; r &= r - 1) {        /* scatter */
                int v = __builtin_ctz(r);
                uint32_t c = row[v];
                if (!c) continue;
                for (uint32_t s = out[v] & ~(uint32_t)m; s; s &= s - 1) {
                    int w = __builtin_ctz(s);
                    uint32_t *t = dp + (m | ((size_t)1 << w)) * n + w;
                    uint32_t y = *t + c; if (y >= p) y -= (uint32_t)p; *t = y;
                }
            }
        }
        printf("%d %u\n", x, T[N - 1]); fflush(stdout);
    }
    free(dp); free(T);
    return 0;
}
