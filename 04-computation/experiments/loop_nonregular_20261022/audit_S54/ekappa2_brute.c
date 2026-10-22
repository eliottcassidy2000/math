/* Independent auditor code: E kappa^2 for uniformly random tournaments by brute force over ALL 2^C(n,2)
 * tournaments (Held-Karp per tournament), and E kappa by the same sweep.  Prints exact integer sums
 *   sum_T H(T), sum_T H(T)^2  (as unsigned 128-bit split into decimal)  for n = 3..nmax.
 * E kappa^2 = (2^(n-1)/n!)^2 * sum H^2 / 2^C(n,2) is then formed with exact rationals in Python.
 * Usage: ekappa2_brute nmax                                                                          */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static void print_u128(unsigned __int128 x) {
    char s[64]; int k = 0;
    if (x == 0) { printf("0"); return; }
    while (x) { s[k++] = '0' + (int)(x % 10); x /= 10; }
    while (k) putchar(s[--k]);
}

int main(int argc, char **argv) {
    int nmax = argc > 1 ? atoi(argv[1]) : 6;
    for (int n = 3; n <= nmax; n++) {
        int E = n * (n - 1) / 2;
        int ei[64], ej[64], c = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { ei[c] = i; ej[c] = j; c++; }
        unsigned out[16];
        unsigned __int128 s1 = 0, s2 = 0;
        int N = 1 << n;
        unsigned long long *dp = malloc(sizeof(unsigned long long) * N * n);
        for (uint64_t code = 0; code < ((uint64_t)1 << E); code++) {
            for (int v = 0; v < n; v++) out[v] = 0;
            for (int e = 0; e < E; e++) {
                if ((code >> e) & 1) out[ei[e]] |= 1u << ej[e]; else out[ej[e]] |= 1u << ei[e];
            }
            for (int i = 0; i < N * n; i++) dp[i] = 0;
            for (int v = 0; v < n; v++) dp[(1 << v) * n + v] = 1;
            for (int m = 1; m < N; m++)
                for (int v = 0; v < n; v++) {
                    unsigned long long x = dp[m * n + v];
                    if (!x) continue;
                    unsigned nb = out[v] & ~m;
                    while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(m | (1 << w)) * n + w] += x; }
                }
            unsigned long long H = 0;
            for (int v = 0; v < n; v++) H += dp[(N - 1) * n + v];
            s1 += H; s2 += (unsigned __int128)H * H;
        }
        free(dp);
        printf("n=%d tournaments=2^%d sumH=", n, E); print_u128(s1); printf(" sumH2="); print_u128(s2); printf("\n");
        fflush(stdout);
    }
    return 0;
}
