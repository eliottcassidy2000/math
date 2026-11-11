/* hp.c (AUDIT, independent) - Hamiltonian path counter.
 * stdin: one tournament (or digraph) per line:  n r_0 r_1 ... r_{n-1}   (r_i = bitmask of out-neighbours of i)
 * stdout: H per line (decimal).  DP over subsets: f[mask][v] = number of Hamiltonian paths of mask ending at v. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
int main(void) {
    int n;
    static u64 *f = NULL; static size_t cap = 0;
    while (scanf("%d", &n) == 1) {
        u64 out[24];
        for (int i = 0; i < n; i++) { if (scanf("%llu", &out[i]) != 1) return 1; }
        size_t need = ((size_t)1 << n) * (size_t)n;
        if (need > cap) { free(f); f = (u64 *)malloc(need * sizeof(u64)); cap = need; if (!f) { fprintf(stderr, "oom\n"); return 2; } }
        memset(f, 0, need * sizeof(u64));
        for (int v = 0; v < n; v++) f[((size_t)1 << v) * n + v] = 1;
        u64 full = ((u64)1 << n) - 1;
        for (u64 mask = 1; mask <= full; mask++) {
            for (int v = 0; v < n; v++) {
                u64 c = f[mask * n + v]; if (!c) continue;
                u64 nb = out[v] & ~mask;
                while (nb) { int u = __builtin_ctzll(nb); nb &= nb - 1; f[(mask | ((u64)1 << u)) * n + u] += c; }
            }
        }
        unsigned __int128 tot = 0; for (int v = 0; v < n; v++) tot += f[full * n + v];
        /* print 128-bit */
        char buf[64]; int p = 63; buf[p] = 0; if (tot == 0) buf[--p] = '0';
        while (tot) { buf[--p] = '0' + (int)(tot % 10); tot /= 10; }
        printf("%s\n", buf + p);
    }
    return 0;
}
