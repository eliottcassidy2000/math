/* mac-mini-2026-11-06-microcosm. Note section 3: Hamiltonian-cycle counts of strong tournaments (gentourng -c), and a gate check for the classes with HC = 7 or 21. */
/* Hamiltonian-cycle counts of strong tournaments read from gentourng -c ascii output on stdin (order n given). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char **argv) {
    int n = atoi(argv[1]); char line[256]; static char seen[400000]; int maxv = 0;
    long long *dp = malloc(sizeof(long long) * (1 << n) * n); int A[16][16];
    long long nonGatedHits[2] = {0, 0};
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == '>' || strlen(line) < 2) continue;
        int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { int b = line[k++] == '1'; A[i][j] = b; A[j][i] = !b; }
        memset(dp, 0, sizeof(long long) * (1 << n) * n); dp[1 * n + 0] = 1;   /* paths starting at vertex 0 */
        for (int S = 1; S < (1 << n); S++) if (S & 1) for (int v = 0; v < n; v++) { long long x = dp[S * n + v]; if (!x) continue;
            for (int w = 0; w < n; w++) if (!(S >> w & 1) && A[v][w]) dp[(S | (1 << w)) * n + w] += x; }
        long long hc = 0; int full = (1 << n) - 1;
        for (int v = 1; v < n; v++) if (A[v][0]) hc += dp[full * n + v];
        if (hc < 400000) seen[hc] = 1; if (hc > maxv && hc < 400000) maxv = hc;
        /* gate check: vertex a with out-degree 1 whose out-neighbour b has in-degree 1 */
        if (hc == 7 || hc == 21) { int gated = 0;
            for (int a = 0; a < n; a++) { int od = 0, bb = -1; for (int y = 0; y < n; y++) if (y != a && A[a][y]) { od++; bb = y; }
                if (od == 1) { int id = 0; for (int x = 0; x < n; x++) if (x != bb && A[x][bb]) id++; if (id == 1) gated = 1; } }
            if (!gated) nonGatedHits[hc == 21]++; }
    }
    printf("n = %d: HC values present among strong tournaments (<= 60):", n);
    for (int h = 1; h <= 60; h++) if (seen[h]) printf(" %d", h);
    printf("\n   HC = 7 present: %s (gate-free examples %lld); HC = 21 present: %s (gate-free examples %lld)\n",
           seen[7] ? "yes" : "no", nonGatedHits[0], seen[21] ? "yes" : "no", nonGatedHits[1]);
    return 0;
}
