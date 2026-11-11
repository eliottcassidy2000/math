/* AUDIT program (independent): number of Hamiltonian paths of a tournament (or any digraph), n <= 24,
 * modulo M, by dynamic programming over (vertex set, last vertex) with compact rows
 * (one 32-bit word per member of each set: 4 * n * 2^(n-1) bytes; n = 23: 386 MB).
 *
 * usage: aud_ham M < input        M = 0: modulo 2^32;  otherwise modulo M (M < 2^32)
 * input: n, then n rows of n characters 0/1 (row i, column j: 1 iff i -> j).
 * Running it with M = 0 and with the prime M = 4294967291 gives the exact count by the Chinese
 * remainder theorem as long as the count is below 2^32 * 4294967291 (about 1.8e19).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

int main(int argc, char **argv) {
    uint64_t M = argc > 1 ? strtoull(argv[1], 0, 10) : 0;
    int n, i, j;
    if (scanf("%d", &n) != 1 || n < 1 || n > 24) { fprintf(stderr, "bad n\n"); return 1; }
    static uint32_t out[32];
    char row[64];
    for (i = 0; i < n; i++) {
        if (scanf("%63s", row) != 1 || (int)strlen(row) != n) { fprintf(stderr, "bad row %d\n", i); return 1; }
        out[i] = 0;
        for (j = 0; j < n; j++) if (row[j] == '1') out[i] |= 1u << j;
    }
    uint32_t N = 1u << n, mask;
    uint32_t *off = malloc((size_t)(N + 1) * sizeof(uint32_t));
    if (!off) { fprintf(stderr, "no memory\n"); return 1; }
    off[0] = 0;
    for (mask = 0; mask < N; mask++) off[mask + 1] = off[mask] + (uint32_t)__builtin_popcount(mask);
    uint32_t *dp = calloc((size_t)off[N], sizeof(uint32_t));
    if (!dp) { fprintf(stderr, "no memory\n"); return 1; }
    for (i = 0; i < n; i++) dp[off[1u << i]] = 1;
    for (mask = 1; mask < N; mask++) {
        uint32_t rem = ~mask & (N - 1), mm, idx = 0;
        const uint32_t *rowp = dp + off[mask];
        for (mm = mask; mm; mm &= mm - 1, idx++) {
            uint32_t x = rowp[idx];
            if (!x) continue;
            int u = __builtin_ctz(mm);
            uint32_t t;
            for (t = out[u] & rem; t; t &= t - 1) {
                int w = __builtin_ctz(t);
                uint32_t nm = mask | (1u << w);
                uint32_t pos = (uint32_t)__builtin_popcount(nm & ((1u << w) - 1));
                uint32_t *cell = dp + off[nm] + pos;
                if (M == 0) *cell += x;
                else { uint64_t s = (uint64_t)*cell + x; if (s >= M) s -= M; *cell = (uint32_t)s; }
            }
        }
    }
    uint64_t h = 0;
    for (i = 0; i < n; i++) {
        h += dp[off[N - 1] + i];
        if (M) h %= M; else h &= 0xffffffffULL;
    }
    printf("%llu\n", (unsigned long long)h);
    return 0;
}
