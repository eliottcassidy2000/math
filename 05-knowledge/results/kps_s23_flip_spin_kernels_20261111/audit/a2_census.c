/* a2_census.c (AUDIT, independent) - visit matrices of all two-port kernels of order m by plain recursive
 * enumeration of walks (no dynamic programming, no cap on the number of visits).
 * stdin: lines "code H" (data/cls_m.txt; code bit b, LSB first, for the pair (i<j) in lexicographic order, 1 iff i -> j)
 * stdout per (class, ordered ports x != y):  m code x y maxi maxj N11 N12 N21 N22 N13.. flag
 *   printed: m code x y maxi maxj N11 N12 N21 N22 det two
 * usage: a2_census m < cls_m.txt */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MV 10
static int m, X, Y, adj[MV][MV], fullmask;
static u64 N[MV + 1][MV + 1];
static void rec(int v, int mask, int i, int j) {
    if (mask == fullmask) N[i][j]++;
    for (int u = 0; u < m; u++) {
        if (!adj[v][u]) continue;
        if (u == X) rec(u, mask, i + 1, j);
        else if (u == Y) rec(u, mask, i, j + 1);
        else if (!(mask >> u & 1)) rec(u, mask | (1 << u), i, j);
    }
}
int main(int argc, char **argv) {
    m = atoi(argv[1]); char line[1 << 12];
    while (fgets(line, sizeof line, stdin)) {
        u64 code = strtoull(line, NULL, 10); memset(adj, 0, sizeof adj); int b = 0;
        for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if (code >> b & 1) adj[i][j] = 1; else adj[j][i] = 1; b++; }
        for (X = 0; X < m; X++) for (Y = 0; Y < m; Y++) { if (X == Y) continue;
            fullmask = 0; for (int v = 0; v < m; v++) if (v != X && v != Y) fullmask |= 1 << v;
            memset(N, 0, sizeof N);
            for (int v = 0; v < m; v++) { if (v == X) rec(v, 0, 1, 0); else if (v == Y) rec(v, 0, 0, 1); else rec(v, 1 << v, 0, 0); }
            int maxi = 0, maxj = 0;
            for (int i = 1; i <= MV; i++) for (int j = 1; j <= MV; j++) if (N[i][j]) { if (i > maxi) maxi = i; if (j > maxj) maxj = j; }
            long long det = (long long)N[1][1] * (long long)N[2][2] - (long long)N[1][2] * (long long)N[2][1];
            printf("%d %llu %d %d %d %d %llu %llu %llu %llu %lld %d\n", m, code, X, Y, maxi, maxj, N[1][1], N[1][2], N[2][1], N[2][2], det, (maxi <= 2 && maxj <= 2));
        }
    }
    return 0;
}
