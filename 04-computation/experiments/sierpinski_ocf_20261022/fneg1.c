/* F(T,-1)/n! for the rotational tournament R_n: the alternating forward-step count
   F(T,-1) = sum over orderings sigma of (-1)^(#backward steps) = N_(n-1)(S) (signed Hamiltonian paths of K_n),
   by subset DP in long double.  Usage: fneg1 n.  mac-mini-2026-10-22-sierp. */
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv) {
    int n = atoi(argv[1]), d = (n - 1) / 2;
    int S[32][32];
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) {
        int w = ((v - u) % n + n) % n;
        S[u][v] = (w == 0) ? 0 : (w <= d ? 1 : -1);
    }
    size_t M = (size_t)1 << n;
    long double *dp = calloc(M * n, sizeof(long double));
    for (int v = 0; v < n; v++) dp[((size_t)1 << v) * n + v] = 1;
    for (size_t m = 1; m < M; m++)
        for (int v = 0; v < n; v++) {
            long double c = dp[m * n + v];
            if (c == 0) continue;
            for (int w = 0; w < n; w++) if (!(m >> w & 1)) dp[(m | ((size_t)1 << w)) * n + w] += c * S[v][w];
        }
    long double tot = 0, fact = 1;
    for (int v = 0; v < n; v++) tot += dp[(M - 1) * n + v];
    for (int i = 2; i <= n; i++) fact *= i;
    printf("n=%d F(R_n,-1)=%.0Lf F/n!=%.10Le\n", n, tot, tot / fact);
    return 0;
}
