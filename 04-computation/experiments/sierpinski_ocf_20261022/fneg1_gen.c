/* F(T,-1) = sum over orderings of (-1)^(#backward steps) = N_(n-1)(S) for a tournament given as an upper-triangle
   code on stdin (one per line), subset DP in long double.  mac-mini-2026-10-22-sierp. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(void) {
    static char line[4096];
    while (fgets(line, sizeof line, stdin)) {
        int L = (int)strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
        if (!L) continue;
        int n = 0; while (n * (n - 1) / 2 < L) n++;
        static int S[32][32]; int p = 0;
        for (int i = 0; i < n; i++) { S[i][i] = 0; for (int j = i + 1; j < n; j++, p++) { int s = line[p] == '1' ? 1 : -1; S[i][j] = s; S[j][i] = -s; } }
        size_t M = (size_t)1 << n;
        long double *dp = calloc(M * n, sizeof(long double));
        for (int v = 0; v < n; v++) dp[((size_t)1 << v) * n + v] = 1;
        for (size_t m = 1; m < M; m++) for (int v = 0; v < n; v++) { long double c = dp[m * n + v]; if (c == 0) continue;
            for (int w = 0; w < n; w++) if (!(m >> w & 1)) dp[(m | ((size_t)1 << w)) * n + w] += c * S[v][w]; }
        long double tot = 0, f = 1; for (int v = 0; v < n; v++) tot += dp[(M - 1) * n + v];
        for (int i = 2; i <= n; i++) f *= i;
        printf("n=%d F(T,-1)=%.0Lf F/n!=%.6Le\n", n, tot, tot / f); fflush(stdout); free(dp);
    }
    return 0;
}
