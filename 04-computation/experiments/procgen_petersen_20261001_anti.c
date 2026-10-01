/* procgen_petersen_20261001_anti.c
   Arc parities of Hamiltonian-path counts for an ANTI-CIRCULANT tournament on Z/2m (m odd):
   x ↦ x + 1 is an anti-automorphism, i.e. T(x, y) = T(y + 1, x + 1).
   Then g (paths starting at b on S) is a re-indexing of f (paths ending at b + 1 on S + 1), so one bitset array of
   2^(2m) words suffices (2m = 26 needs 256 MB).  Arc orbits under sigma: P -> reverse(P + 1) are indexed by the
   difference d = 1..m (the arc between 0 and d is a representative; d = m is the antipodal orbit).

   stdin: n = 2m, then n rows of n characters '0'/'1' (row i, column j = 1 iff i -> j).
   stdout: "H2 <H mod 2>" and, for d = 1..m, "d <a> <b> <parity of c(a -> b)>" for the arc between 0 and d. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(void) {
    int n;
    if (scanf("%d", &n) != 1 || n < 2 || n > 26 || n % 4 != 2) return 1;
    int m = n / 2;
    uint32_t in[32] = {0}, out[32] = {0};
    int T[32][32];
    char row[64];
    for (int i = 0; i < n; i++) {
        if (scanf("%63s", row) != 1) return 1;
        for (int j = 0; j < n; j++) {
            T[i][j] = (row[j] == '1');
            if (T[i][j]) { out[i] |= 1u << j; in[j] |= 1u << i; }
        }
    }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        if (i == j) { if (T[i][j]) return 3; continue; }
        if (T[i][j] + T[j][i] != 1) return 4;                       /* not a tournament */
        if (T[i][j] != T[(j + 1) % n][(i + 1) % n]) return 5;       /* not anti-circulant */
    }
    size_t N = (size_t)1 << n;
    uint32_t full = (uint32_t)(N - 1);
    uint32_t *f = calloc(N, sizeof(uint32_t));
    if (!f) return 2;
    for (int v = 0; v < n; v++) f[(size_t)1 << v] = 1u << v;
    for (size_t S = 1; S < N; S++) {
        if (!(S & (S - 1))) continue;
        uint32_t fs = 0, rest = (uint32_t)S;
        while (rest) {
            int w = __builtin_ctz(rest);
            rest &= rest - 1;
            if (__builtin_popcount(f[S ^ ((size_t)1 << w)] & in[w]) & 1) fs |= 1u << w;
        }
        f[S] = fs;
    }
    printf("H2 %d\n", __builtin_popcount(f[N - 1]) & 1);
    for (int d = 1; d <= m; d++) {
        int a = 0, b = d;
        if (!T[a][b]) { a = d; b = 0; }
        int b1 = (b + 1) % n;
        int par = 0;
        for (size_t S = 1; S < N; S++) {
            if (!((S >> a) & 1) || ((S >> b) & 1)) continue;
            if (!((f[S] >> a) & 1)) continue;
            uint32_t C = full ^ (uint32_t)S;                         /* V \ S, contains b */
            uint32_t R = ((C << 1) | (C >> (n - 1))) & full;           /* C + 1 */
            par ^= (int)((f[R] >> b1) & 1u);
        }
        printf("%d %d %d %d\n", d, a, b, par);
    }
    free(f);
    return 0;
}
