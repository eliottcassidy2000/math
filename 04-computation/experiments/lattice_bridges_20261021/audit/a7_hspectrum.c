/* a7_hspectrum.c -- AUDIT of HYP-9362 (d) Phi_3 remarks: exhaustive H-spectrum of all labelled tournaments on
 * n <= 7 vertices (2^C(n,2) of them), with an explicit witness (arc bitmask, upper-triangle order i<j, bit=1 means i->j)
 * for H = 13 = Phi_3(3), 31 = Phi_3(5), 43 = Phi_3(6), and a check that 7 = Phi_3(2) and 21 = Phi_3(4) never occur.
 * Usage: ./a7_hspectrum > a7_hspectrum.out
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
static int n; static uint32_t outm[8];
static long ham(void) {
    static long f[1 << 7][7];
    memset(f, 0, sizeof f);
    for (int v = 0; v < n; v++) f[1 << v][v] = 1;
    for (int mask = 1; mask < (1 << n); mask++)
        for (int v = 0; v < n; v++) {
            long c = f[mask][v]; if (!c) continue;
            uint32_t nx = outm[v] & ~(uint32_t)mask;
            while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; f[mask | (1 << w)][w] += c; }
        }
    long t = 0; for (int v = 0; v < n; v++) t += f[(1 << n) - 1][v];
    return t;
}
int main(void) {
    int targets[5] = {7, 21, 13, 31, 43};
    for (n = 3; n <= 7; n++) {
        int E = n * (n - 1) / 2; long total = 1L << E;
        static long cnt[700]; memset(cnt, 0, sizeof cnt);
        long wit[5]; for (int t = 0; t < 5; t++) wit[t] = -1;
        int pi[32], pj[32], k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { pi[k] = i; pj[k] = j; k++; }
        for (long b = 0; b < total; b++) {
            for (int v = 0; v < n; v++) outm[v] = 0;
            for (int e = 0; e < E; e++) { if ((b >> e) & 1) outm[pi[e]] |= 1u << pj[e]; else outm[pj[e]] |= 1u << pi[e]; }
            long H = ham();
            if (H < 700) cnt[H]++;
            for (int t = 0; t < 5; t++) if (H == targets[t] && wit[t] < 0) wit[t] = b;
        }
        printf("n=%d: %ld labelled tournaments; attained H values:", n, total);
        for (int h = 0; h < 700; h++) if (cnt[h]) printf(" %d", h);
        printf("\n");
        for (int t = 0; t < 5; t++) {
            printf("   H=%d: %s", targets[t], wit[t] < 0 ? "not attained" : "attained, witness arcs:");
            if (wit[t] >= 0) { for (int e = 0; e < E; e++) printf(" %d%s%d", ((wit[t] >> e) & 1) ? pi[e] : pj[e], "->", ((wit[t] >> e) & 1) ? pj[e] : pi[e]); }
            printf("  (labelled count %ld)\n", cnt[targets[t]]);
        }
        fflush(stdout);
    }
    return 0;
}
