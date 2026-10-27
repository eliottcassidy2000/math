/* median_checks.c -- mac-mini-2026-10-27-dipole
 * Median orders (orders minimising the number of backward arcs) are Hamiltonian paths, and by the feedback property
 * their first vertex has score >= (n-1)/2 and their last vertex score <= (n-1)/2.  Tested here, for every class:
 *   (M1) T not regular (n odd) / not almost regular: some median order starts above (n-1)/2 and ends below it;
 *   (M2) some median order starts at a maximum-score vertex and ends at a minimum-score vertex;
 *   (M3) the weaker (C7) by full Hamiltonian-path endpoints is tested in endpoint_checks.c.
 * min backward arcs over orders from a to b = g_a[V \ {b}] + s(b), g_a[S] = min over orders of S starting at a.
 * Build: cc -O3 -march=native median_checks.c -o median_checks ;  Run: gentourng -q n | ./median_checks n
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n; static uint32_t out_[16];
static int g[1 << 13];
int main(int argc, char **argv) {
    n = atoi(argv[1]); int full = (1 << n) - 1;
    char line[256]; long classes = 0, nonreg = 0, m1_bad = 0, m2_bad = 0, m1n = 0;
    static int best_ab[16][16];
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        classes++;
        for (int i = 0; i < n; i++) out_[i] = 0;
        int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') out_[i] |= 1u << j; else out_[j] |= 1u << i; k++; }
        int s[16], mn = 99, mx = -1; for (int v = 0; v < n; v++) { s[v] = __builtin_popcount(out_[v]); if (s[v] < mn) mn = s[v]; if (s[v] > mx) mx = s[v]; }
        if (mn == mx) continue;
        nonreg++;
        int gmin = 1 << 30;
        for (int a = 0; a < n; a++) {
            for (int S = 0; S <= full; S++) g[S] = 1 << 29;
            g[1 << a] = 0;
            for (int S = 1; S <= full; S++) { if (!(S >> a & 1) || g[S] >= (1 << 29)) continue;
                for (int w = 0; w < n; w++) if (!(S >> w & 1)) { int c = g[S] + __builtin_popcount(out_[w] & S); if (c < g[S | 1 << w]) g[S | 1 << w] = c; } }
            for (int b = 0; b < n; b++) { if (b == a) { best_ab[a][b] = 1 << 30; continue; } best_ab[a][b] = g[full & ~(1 << b)] + s[b]; if (best_ab[a][b] < gmin) gmin = best_ab[a][b]; }
        }
        int ok1 = 0, ok2 = 0, haveHL = 0;
        for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if (a != b && best_ab[a][b] == gmin) {
            if (2 * s[a] > n - 1 && 2 * s[b] < n - 1) ok1 = 1;
            if (s[a] == mx && s[b] == mn) ok2 = 1; }
        for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if (2 * s[a] > n - 1 && 2 * s[b] < n - 1) haveHL = 1;
        if (haveHL) { m1n++; if (!ok1) m1_bad++; }
        if (!ok2) m2_bad++;
    }
    printf("n=%d classes=%ld non-regular=%ld\n", n, classes, nonreg);
    printf("(M1) classes with a vertex above and below (n-1)/2: %ld; no median order from above to below: %ld\n", m1n, m1_bad);
    printf("(M2) no median order from a maximum-score vertex to a minimum-score vertex: %ld\n", m2_bad);
    return 0;
}
