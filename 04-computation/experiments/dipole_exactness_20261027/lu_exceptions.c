/* lu_exceptions.c -- structure of the max->min arcs x -> y (x of maximum score, y of minimum score) that lie on no
 * Hamiltonian path from x to y: multiplicities of the extreme scores, connectivity, scores of x, y, and whether the class
 * still has a Hamiltonian path from SOME maximum-score vertex to SOME minimum-score vertex.
 * Build: cc -O3 lu_exceptions.c -o lu_exceptions ; Run: gentourng -q n | ./lu_exceptions n */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n; static uint32_t out_[16]; static uint16_t pi_[16]; static uint16_t dp[1 << 13][13];
static void endpoints(void) { int full = (1 << n) - 1;
    for (int S = 0; S <= full; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1 << v;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) { uint16_t st = dp[S][v]; if (!st) continue; uint32_t nx = out_[v] & ~S; while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dp[S | 1 << w][w] |= st; } }
    for (int a = 0; a < n; a++) pi_[a] = 0;
    for (int b = 0; b < n; b++) { uint16_t st = dp[full][b]; while (st) { int a = __builtin_ctz(st); st &= st - 1; pi_[a] |= 1 << b; } } }
int main(int argc, char **argv) {
    n = atoi(argv[1]); char line[256]; long ncls = 0; long hist[16][16] = {{0}}; long rescued = 0, unique_both = 0;
    long xdeg_hist[16] = {0}; int shown = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        for (int i = 0; i < n; i++) out_[i] = 0;
        int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') out_[i] |= 1u << j; else out_[j] |= 1u << i; k++; }
        int s[16], mn = 99, mx = -1; for (int v = 0; v < n; v++) { s[v] = __builtin_popcount(out_[v]); if (s[v] < mn) mn = s[v]; if (s[v] > mx) mx = s[v]; }
        if (mn == mx) continue;
        endpoints();
        int bad = 0, good = 0, nmax = 0, nmin = 0;
        for (int v = 0; v < n; v++) { nmax += s[v] == mx; nmin += s[v] == mn; }
        for (int x = 0; x < n; x++) for (int y = 0; y < n; y++) if (s[x] == mx && s[y] == mn) { if (pi_[x] >> y & 1) good = 1; else if (out_[x] >> y & 1) bad = 1; }
        if (!bad) continue;
        ncls++; hist[nmax][nmin]++; if (good) rescued++; if (nmax == 1 && nmin == 1) unique_both++;
        xdeg_hist[mx]++;
        if (shown < 6) { shown++; printf("example: scores"); for (int v = 0; v < n; v++) printf(" %d", s[v]); printf("  arcs %s", line); }
    }
    printf("n=%d classes with a max->min arc on no max->min Hamiltonian path: %ld; rescued by another extreme pair: %ld; unique max and unique min: %ld\n", n, ncls, rescued, unique_both);
    printf("  by (#max-score vertices, #min-score vertices):"); for (int a = 0; a < 16; a++) for (int b = 0; b < 16; b++) if (hist[a][b]) printf(" (%d,%d):%ld", a, b, hist[a][b]); printf("\n");
    printf("  by maximum score:"); for (int a = 0; a < 16; a++) if (xdeg_hist[a]) printf(" %d:%ld", a, xdeg_hist[a]); printf("\n");
    return 0;
}
