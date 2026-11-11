/* kernels_sat.c - the visit matrix of two-port kernels with saturating counters (fast version of kernels.c,
 * used for order 8).  For every class K of order m (code as in kernels.c) and every pair of ports x < y:
 *   walks of K that use every vertex other than x, y exactly once; visits to x and to y are counted up to 3
 *   (3 means "three or more").  N[i][j] for i, j <= 2 are exact; flag two = 1 iff no walk visits a port three
 *   times, and then H(K[x <- A, y <- B]) = sum_{i,j <= 2} N[i][j] P_i(A) P_j(B) for all tournaments A, B.
 * Output per (class, x < y):  m code x y  N11 N12 N21 N22  det  two
 * (the pair (y, x) has the transposed matrix).   usage: kernels_sat m < cls_m.txt */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MV 9
static u64 f[2][1 << (MV - 2)][MV][4][4];
int main(int argc, char **argv) {
    int m = atoi(argv[1]); char line[1 << 12];
    while (fgets(line, sizeof line, stdin)) {
        u64 code = strtoull(line, NULL, 10); int adj[MV][MV]; memset(adj, 0, sizeof adj); int b = 0;
        for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if (code >> b & 1) adj[i][j] = 1; else adj[j][i] = 1; b++; }
        for (int x = 0; x < m; x++) for (int y = x + 1; y < m; y++) {
            int idx[MV], np = 0; for (int v = 0; v < m; v++) idx[v] = (v == x || v == y) ? -1 : np++;
            int full = (1 << np) - 1;
            u64 N[4][4]; memset(N, 0, sizeof N);
            size_t sz = sizeof(u64) * (size_t)(1 << np) * MV * 16;
            memset(f[0], 0, sizeof f[0]);
            int cur = 0;
            for (int v = 0; v < m; v++) { if (v == x) f[0][0][v][1][0] = 1; else if (v == y) f[0][0][v][0][1] = 1; else f[0][1 << idx[v]][v][0][0] = 1; }
            int maxlen = np + 2 * (m - 1);          /* a port is visited at most m - 1 times */
            (void)sz;
            for (int len = 1; len <= maxlen; len++) {
                for (int i = 1; i < 4; i++) for (int j = 1; j < 4; j++) for (int v = 0; v < m; v++) N[i][j] += f[cur][full][v][i][j];
                int nx = cur ^ 1; memset(f[nx], 0, sizeof f[nx]);
                int any = 0;
                for (int mask = 0; mask <= full; mask++) for (int v = 0; v < m; v++) for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
                    u64 c = f[cur][mask][v][i][j]; if (!c) continue;
                    for (int u = 0; u < m; u++) { if (!adj[v][u]) continue;
                        if (u == x) { int i2 = i < 3 ? i + 1 : 3; f[nx][mask][u][i2][j] += c; any = 1; }
                        else if (u == y) { int j2 = j < 3 ? j + 1 : 3; f[nx][mask][u][i][j2] += c; any = 1; }
                        else if (!(mask >> idx[u] & 1)) { f[nx][mask | 1 << idx[u]][u][i][j] += c; any = 1; } } }
                cur = nx;
                if (!any) break;
            }
            int two = 1; for (int i = 1; i < 4; i++) for (int j = 1; j < 4; j++) if ((i == 3 || j == 3) && N[i][j]) two = 0;
            long long det = (long long)N[1][1] * (long long)N[2][2] - (long long)N[1][2] * (long long)N[2][1];
            printf("%d %llu %d %d  %llu %llu %llu %llu  %lld  %d\n", m, code, x, y, N[1][1], N[1][2], N[2][1], N[2][2], det, two);
        }
    }
    return 0;
}
