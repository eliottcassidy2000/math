/* AUDIT program (independent): number of Hamiltonian paths that START at a given vertex, modulo 256,
 * for a tournament on n <= 27 vertices (dynamic programme over subsets of the other n-1 vertices with compact
 * one-byte rows: about (n-1) 2^(n-2) bytes plus a 4 * 2^(n-1) byte offset table; n = 26: 0.55 GB).
 * usage: aud_ham8 START < input     (input: n, then n rows of 0/1)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

int main(int argc, char **argv) {
    int start = argc > 1 ? atoi(argv[1]) : 0;
    int n, i, j;
    char row[64];
    static int adj[32][32];
    if (scanf("%d", &n) != 1 || n < 2 || n > 27) return 1;
    for (i = 0; i < n; i++) {
        if (scanf("%63s", row) != 1) return 1;
        for (j = 0; j < n; j++) adj[i][j] = row[j] == '1';
    }
    int m = n - 1, vs[32], k = 0;
    for (i = 0; i < n; i++) if (i != start) vs[k++] = i;
    uint32_t o[32], first = 0;
    for (i = 0; i < m; i++) {
        o[i] = 0;
        for (j = 0; j < m; j++) if (adj[vs[i]][vs[j]]) o[i] |= 1u << j;
        if (adj[start][vs[i]]) first |= 1u << i;
    }
    uint32_t NM = 1u << m, mask;
    uint32_t *off = malloc((size_t)(NM + 1) * sizeof(uint32_t));
    if (!off) { fprintf(stderr, "no memory\n"); return 1; }
    off[0] = 0;
    for (mask = 0; mask < NM; mask++) off[mask + 1] = off[mask] + (uint32_t)__builtin_popcount(mask);
    uint8_t *dp = calloc((size_t)off[NM], 1);
    if (!dp) { fprintf(stderr, "no memory\n"); return 1; }
    for (i = 0; i < m; i++) if (first >> i & 1) dp[off[1u << i]] = 1;
    for (mask = 1; mask < NM; mask++) {
        uint32_t rem = ~mask & (NM - 1), mm, idx = 0;
        const uint8_t *rowp = dp + off[mask];
        for (mm = mask; mm; mm &= mm - 1, idx++) {
            uint8_t x = rowp[idx];
            if (!x) continue;
            int u = __builtin_ctz(mm);
            uint32_t t;
            for (t = o[u] & rem; t; t &= t - 1) {
                int w = __builtin_ctz(t);
                uint32_t nm = mask | (1u << w);
                dp[off[nm] + (uint32_t)__builtin_popcount(nm & ((1u << w) - 1))] += x;
            }
        }
    }
    unsigned h = 0;
    for (i = 0; i < m; i++) h += dp[off[NM - 1] + i];
    printf("%u\n", h & 255u);
    return 0;
}
