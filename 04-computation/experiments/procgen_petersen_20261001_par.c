/* procgen_petersen_20261001_par.c
   Parities (mod 2) of Hamiltonian-path counts of a digraph: H mod 2 and, for every arc, the parity of the
   number of directed Hamiltonian paths through it.  Bitset subset DP, n <= 24 (memory 2 * 2^n * 4 bytes).

   stdin: n, then n rows of n characters '0'/'1' (row i, column j = 1 iff arc i -> j).
   stdout: "H2 <H mod 2>" then one line per arc "i j <parity>" in row-major order.

   f[S] bit v = parity of the number of directed paths with vertex set S that end at v,
   g[S] bit v = parity of the number of directed paths with vertex set S that start at v (pull recursions);
   c(a -> b) = sum over S with a in S, b not in S of f[S]_a * g[V \ S]_b (mod 2). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(void) {
    int n;
    if (scanf("%d", &n) != 1 || n < 1 || n > 24) return 1;
    uint32_t in[32] = {0}, out[32] = {0};
    char row[64];
    for (int i = 0; i < n; i++) {
        if (scanf("%63s", row) != 1) return 1;
        for (int j = 0; j < n; j++) if (row[j] == '1') { out[i] |= 1u << j; in[j] |= 1u << i; }
    }
    size_t N = (size_t)1 << n;
    uint32_t *f = calloc(N, sizeof(uint32_t)), *g = calloc(N, sizeof(uint32_t));
    if (!f || !g) return 2;
    for (int v = 0; v < n; v++) { f[(size_t)1 << v] = 1u << v; g[(size_t)1 << v] = 1u << v; }
    for (size_t S = 1; S < N; S++) {
        if (!(S & (S - 1))) continue;           /* singletons already set */
        uint32_t fs = 0, gs = 0;
        uint32_t rest = (uint32_t)S;
        while (rest) {
            int w = __builtin_ctz(rest);
            rest &= rest - 1;
            size_t T = S ^ ((size_t)1 << w);
            if (__builtin_popcount(f[T] & in[w]) & 1) fs |= 1u << w;
            if (__builtin_popcount(g[T] & out[w]) & 1) gs |= 1u << w;
        }
        f[S] = fs; g[S] = gs;
    }
    printf("H2 %d\n", __builtin_popcount(f[N - 1]) & 1);
    size_t full = N - 1;
    for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) {
        if (!((out[a] >> b) & 1)) continue;
        int par = 0;
        for (size_t S = 1; S < N; S++) {
            if (!((S >> a) & 1) || ((S >> b) & 1)) continue;
            par ^= (int)((f[S] >> a) & (g[full ^ S] >> b) & 1u);
        }
        printf("%d %d %d\n", a, b, par);
    }
    free(f); free(g);
    return 0;
}
