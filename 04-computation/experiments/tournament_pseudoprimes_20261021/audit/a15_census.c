/* AUDIT a15 (independent auditor, 2026-10-21): census of all tournaments on n vertices from `gentourng -q n`
 * (upper triangle row by row; '1' at pair (i,j), i<j, read as i->j; the converse convention gives the same H, C, |Aut|).
 * Per tournament prints: H C strong |Aut|, where
 *   H  = Hamiltonian paths (subset DP), C = directed Hamiltonian cycles (paths from 0 closing to 0),
 *   |Aut| = #{Hamiltonian paths Q : the bijection P0[i] -> Q[i] is an automorphism}, P0 a fixed Hamiltonian path
 *          (Aut acts freely on Hamiltonian paths), enumerated by DFS.
 * Usage: gentourng -q n | ./a15_census n */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n, A[16][16], P0[16], have0, seq[16], aut; static u64 Hdfs;
static void dfs(int depth, int used) {
    if (depth == n) {
        Hdfs++;
        if (!have0) { memcpy(P0, seq, sizeof(int) * n); have0 = 1; }
        int ok = 1;
        for (int i = 0; i < n && ok; i++) for (int j = 0; j < n; j++) if (i != j && A[P0[i]][P0[j]] != A[seq[i]][seq[j]]) { ok = 0; break; }
        aut += ok; return;
    }
    int last = seq[depth - 1];
    for (int w = 0; w < n; w++) if (!(used >> w & 1) && A[last][w]) { seq[depth] = w; dfs(depth + 1, used | 1 << w); }
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); char line[256];
    static u64 dp[1 << 12][12];
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        memset(A, 0, sizeof A); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') A[i][j] = 1; else A[j][i] = 1; k++; }
        /* H by subset DP */
        int N = 1 << n; for (int m = 0; m < N; m++) for (int v = 0; v < n; v++) dp[m][v] = 0;
        for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
        for (int m = 1; m < N; m++) for (int v = 0; v < n; v++) if (dp[m][v]) for (int w = 0; w < n; w++) if (!(m >> w & 1) && A[v][w]) dp[m | 1 << w][w] += dp[m][v];
        u64 H = 0; for (int v = 0; v < n; v++) H += dp[N - 1][v];
        /* C: paths from 0 */
        for (int m = 0; m < N; m++) for (int v = 0; v < n; v++) dp[m][v] = 0;
        dp[1][0] = 1;
        for (int m = 1; m < N; m += 2) for (int v = 0; v < n; v++) if (dp[m][v]) for (int w = 0; w < n; w++) if (!(m >> w & 1) && A[v][w]) dp[m | 1 << w][w] += dp[m][v];
        u64 C = 0; for (int v = 1; v < n; v++) if (A[v][0]) C += dp[N - 1][v];
        if (n == 1) C = 0;
        /* strong: reachability both ways from 0 */
        int r1 = 1, r2 = 1, ch = 1;
        while (ch) { ch = 0; for (int v = 0; v < n; v++) if (r1 >> v & 1) for (int w = 0; w < n; w++) if (A[v][w] && !(r1 >> w & 1)) { r1 |= 1 << w; ch = 1; } }
        ch = 1; while (ch) { ch = 0; for (int v = 0; v < n; v++) if (r2 >> v & 1) for (int w = 0; w < n; w++) if (A[w][v] && !(r2 >> w & 1)) { r2 |= 1 << w; ch = 1; } }
        int strong = (r1 == N - 1 && r2 == N - 1);
        /* |Aut| via path enumeration */
        have0 = 0; aut = 0; Hdfs = 0;
        for (int v = 0; v < n; v++) { seq[0] = v; dfs(1, 1 << v); }
        if (Hdfs != H) { fprintf(stderr, "H mismatch\n"); return 1; }
        printf("%llu %llu %d %d\n", H, C, strong, aut);
    }
    return 0; }
