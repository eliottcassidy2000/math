/* a8_census_dp.c (AUDIT, independent) - census of two-port kernels of order m from the audit's OWN class list
 * (aud_cls_m.txt, gen_classes.c).  For every class and every pair of ports x < y: number of walks visiting every
 * other vertex once, x i times and y j times, for i, j in {1, 2, 3+}, by ONE pass over the states
 * (set of visited non-port vertices, last vertex, i, j) in a topological order (no layering by length):
 * within one set, states ending at a non-port vertex come first, then the port that beats the other, then the other.
 * Validated against the plain walk enumeration (a2_census.c) at orders 3..7 by a8_census_run.py.
 * usage: a8_census_dp m file modulus residue [dump]
 * output: "S m all two join ell par hyp" summary line (ordered port pairs = 2 x unordered) and
 *         "C m code x y al be ga de" for every carrier (al = +-1, be = +-2^u same sign, ga = 0, de > 0), port order (x, y);
 *         with 'dump': "K m code x y N11 N12 N21 N22 two" for every x < y. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
typedef long long i64;
#define MV 9
static u64 f[1 << (MV - 2)][MV][4][4];
static int ispow2(i64 x) { return x > 0 && (x & (x - 1)) == 0; }
int main(int argc, char **argv) {
    int m = atoi(argv[1]); FILE *fp = fopen(argv[2], "r"); int modulus = atoi(argv[3]), residue = atoi(argv[4]); int dump = argc > 5;
    char line[1 << 10]; long idx = -1; i64 all = 0, two_ = 0, join = 0, ell = 0, par = 0, hyp = 0;
    while (fgets(line, sizeof line, fp)) {
        idx++; if (idx % modulus != residue) continue;
        u64 code = strtoull(line, NULL, 10); int adj[MV][MV]; memset(adj, 0, sizeof adj); int b = 0;
        for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if (code >> b & 1) adj[i][j] = 1; else adj[j][i] = 1; b++; }
        for (int x = 0; x < m; x++) for (int y = x + 1; y < m; y++) {
            int idxv[MV], np = 0, vert[MV]; for (int v = 0; v < m; v++) { if (v == x || v == y) idxv[v] = -1; else { idxv[v] = np; vert[np++] = v; } }
            int full = (1 << np) - 1; int p1 = adj[x][y] ? x : y, p2 = adj[x][y] ? y : x;
            memset(f, 0, sizeof(u64) * (size_t)(1 << np) * MV * 16);
            f[0][x][1][0] = 1; f[0][y][0][1] = 1; for (int k = 0; k < np; k++) f[1 << k][vert[k]][0][0] = 1;
            for (int mask = 0; mask <= full; mask++) {
                /* non-port last vertices */
                for (int k = 0; k < np; k++) { if (!(mask >> k & 1)) continue; int v = vert[k];
                    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { u64 c = f[mask][v][i][j]; if (!c) continue;
                        if (adj[v][x]) f[mask][x][i < 3 ? i + 1 : 3][j] += c;
                        if (adj[v][y]) f[mask][y][i][j < 3 ? j + 1 : 3] += c;
                        for (int k2 = 0; k2 < np; k2++) if (!(mask >> k2 & 1) && adj[v][vert[k2]]) f[mask | 1 << k2][vert[k2]][i][j] += c; } }
                /* ports: p1 first (p1 -> p2) */
                for (int t = 0; t < 2; t++) { int v = t ? p2 : p1;
                    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { u64 c = f[mask][v][i][j]; if (!c) continue;
                        if (t == 0) { if (p2 == x) f[mask][x][i < 3 ? i + 1 : 3][j] += c; else f[mask][y][i][j < 3 ? j + 1 : 3] += c; }
                        for (int k2 = 0; k2 < np; k2++) if (!(mask >> k2 & 1) && adj[v][vert[k2]]) f[mask | 1 << k2][vert[k2]][i][j] += c; } }
            }
            u64 N[4][4]; memset(N, 0, sizeof N);
            for (int v = 0; v < m; v++) for (int i = 1; i < 4; i++) for (int j = 1; j < 4; j++) N[i][j] += f[full][v][i][j];
            int two = 1; for (int i = 1; i < 4; i++) for (int j = 1; j < 4; j++) if ((i == 3 || j == 3) && N[i][j]) two = 0;
            i64 n11 = (i64)N[1][1], n12 = (i64)N[1][2], n21 = (i64)N[2][1], n22 = (i64)N[2][2]; i64 det = n11 * n22 - n12 * n21;
            if (dump) printf("K %d %llu %d %d %lld %lld %lld %lld %d\n", m, code, x, y, n11, n12, n21, n22, two);
            all += 2; if (!two) continue; two_ += 2;
            if (n12 == 0 && n21 == 0 && n22 == 0) join += 2; else if (det > 0) ell += 2; else if (det == 0) par += 2; else hyp += 2;
            if (n22 <= 0) continue;
            i64 al = n11 - 2 * n12 - 2 * n21 + 4 * n22, de = n22;
            for (int sw = 0; sw < 2; sw++) { i64 be = (sw ? n12 : n21) - 2 * n22, ga = (sw ? n21 : n12) - 2 * n22;
                if (ga == 0 && ((al == 1 && ispow2(be)) || (al == -1 && ispow2(-be))))
                    printf("C %d %llu %d %d %lld %lld %lld %lld\n", m, code, sw ? y : x, sw ? x : y, al, be, ga, de); }
        }
    }
    printf("S %d %lld %lld %lld %lld %lld %lld\n", m, all, two_, join, ell, par, hyp);
    return 0;
}
