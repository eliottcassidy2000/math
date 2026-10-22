/* Top of the forest series for the rotational tournament R_n:  T_j = N_(n-1-j)/(n)_(n-1-j), j = 0..JM,
   where N_L is the signed linear-forest number (kappa = sum_L N_L/(n)_L, THM-4757).  A forest with L = n-1-j arcs
   has j+1 components and lies in (j+1)! orderings, so N_(n-1-j) = (1/(j+1)!) sum over orderings w of the sum over
   j-subsets of skipped steps of prod over the other steps of S.  Subset DP, long double.  Usage: forest_top n JM.
   mac-mini-2026-10-22-sierp. */
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv) {
    int n = atoi(argv[1]), JM = atoi(argv[2]), d = (n - 1) / 2, K = JM + 1;
    int S[32][32];
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) { int w = ((v - u) % n + n) % n; S[u][v] = w == 0 ? 0 : (w <= d ? 1 : -1); }
    size_t M = (size_t)1 << n;
    long double *dp = calloc(M * n * K, sizeof(long double));
    if (!dp) { fprintf(stderr, "alloc\n"); return 1; }
    for (int v = 0; v < n; v++) dp[(((size_t)1 << v) * n + v) * K + 0] = 1;
    for (size_t m = 1; m < M; m++) for (int v = 0; v < n; v++) {
        long double *c = &dp[(m * n + v) * K];
        int any = 0; for (int j = 0; j < K; j++) if (c[j] != 0) { any = 1; break; }
        if (!any) continue;
        for (int w = 0; w < n; w++) if (!(m >> w & 1)) {
            long double *t = &dp[((m | ((size_t)1 << w)) * n + w) * K];
            for (int j = 0; j < K; j++) {
                t[j] += c[j] * S[v][w];
                if (j + 1 < K) t[j + 1] += c[j];
            }
        }
    }
    long double fact[40]; fact[0] = 1; for (int i = 1; i < 40; i++) fact[i] = fact[i - 1] * i;
    printf("n=%d", n);
    long double cum = 0;
    for (int j = 0; j < K; j++) {
        long double tot = 0; for (int v = 0; v < n; v++) tot += dp[((M - 1) * n + v) * K + j];
        long double N = tot / fact[j + 1];          /* N_(n-1-j) */
        int L = n - 1 - j;
        long double T = N * fact[n - L] / fact[n]; /* N_L/(n)_L = N_L (n-L)!/n! */
        cum += T;
        printf("  j=%d T=%.6Le cum=%.6Le", j, T, cum);
    }
    printf("\n"); free(dp); return 0;
}
