/* kernels.c - E2 of the S23 pre-registration: the visit matrix of a two-port kernel.
 * K a tournament of order m (code: bit b for the pair (i<j) in lexicographic order, 1 iff i -> j), ports x != y.
 * N[i][j] = number of walks in K that use every arc as K gives it, visit every vertex other than x, y exactly once,
 * visit x exactly i times and y exactly j times (i, j >= 1).  Then for all tournaments A, B
 *     H(K[x <- A, y <- B]) = sum_{i,j} N[i][j] P_i(A) P_j(B),
 * P_i(A) = number of ordered i-tuples of disjoint nonempty paths covering A  (P_1 = H, P_2 = Pi_2 - 2H).
 * Output per (class, x, y): m code x y  maxi maxj  N11 N12 N21 N22  det  [flag two = 1 iff no walk visits a port 3 times]
 * usage: kernels m < cls_m.txt */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MV 8
#define CAP 8
static u64 f[2][1 << (MV - 2)][MV][CAP][CAP];
int main(int argc, char **argv) {
    int m = atoi(argv[1]); char line[1 << 12];
    while (fgets(line, sizeof line, stdin)) {
        u64 code = strtoull(line, NULL, 10); int adj[MV][MV]; memset(adj, 0, sizeof adj); int b = 0;
        for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if (code >> b & 1) adj[i][j] = 1; else adj[j][i] = 1; b++; }
        for (int x = 0; x < m; x++) for (int y = 0; y < m; y++) { if (x == y) continue;
            int idx[MV], np = 0; for (int v = 0; v < m; v++) idx[v] = (v == x || v == y) ? -1 : np++;
            int full = (1 << np) - 1, cap = m;            /* a port is visited at most m - 1 times */
            if (cap > CAP) cap = CAP;
            static u64 N[CAP][CAP]; memset(N, 0, sizeof N);
            memset(f, 0, sizeof f);
            int cur = 0;
            for (int v = 0; v < m; v++) { if (v == x) f[0][0][v][1][0] = 1; else if (v == y) f[0][0][v][0][1] = 1; else f[0][1 << idx[v]][v][0][0] = 1; }
            int maxlen = np + 2 * (cap - 1);
            for (int len = 1; len <= maxlen; len++) {
                for (int i = 1; i < cap; i++) for (int j = 1; j < cap; j++) for (int v = 0; v < m; v++) N[i][j] += f[cur][full][v][i][j];
                int nx = cur ^ 1; memset(f[nx], 0, sizeof f[nx]);
                for (int mask = 0; mask <= full; mask++) for (int v = 0; v < m; v++) for (int i = 0; i < cap; i++) for (int j = 0; j < cap; j++) {
                    u64 c = f[cur][mask][v][i][j]; if (!c) continue;
                    for (int u = 0; u < m; u++) { if (!adj[v][u]) continue;
                        if (u == x) { if (i + 1 < cap) f[nx][mask][u][i + 1][j] += c; }
                        else if (u == y) { if (j + 1 < cap) f[nx][mask][u][i][j + 1] += c; }
                        else if (!(mask >> idx[u] & 1)) f[nx][mask | 1 << idx[u]][u][i][j] += c; } }
                cur = nx;
            }
            int maxi = 0, maxj = 0; for (int i = 1; i < cap; i++) for (int j = 1; j < cap; j++) if (N[i][j]) { if (i > maxi) maxi = i; if (j > maxj) maxj = j; }
            long long det = (long long)N[1][1] * (long long)N[2][2] - (long long)N[1][2] * (long long)N[2][1];
            printf("%d %llu %d %d  %d %d  %llu %llu %llu %llu  %lld  %d", m, code, x, y, maxi, maxj, N[1][1], N[1][2], N[2][1], N[2][2], det, (maxi <= 2 && maxj <= 2));
            if (maxi > 2 || maxj > 2) { printf("  |"); for (int i = 1; i <= maxi; i++) { for (int j = 1; j <= maxj; j++) printf(" %llu", N[i][j]); printf(" ;"); } }
            printf("\n"); }
    }
    return 0;
}
