/* Strong tournaments (gentourng -c): for each cycle length k = 3..n, count directed k-cycles with NO balanced arc,
   and the histogram over T of the minimum number of balanced arcs on a Hamiltonian cycle. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n; static int adj[16];
static u64 *f, *g; static u64 cnt[16][16]; static u64 H; static int bal[16][16];
static long long cyc_tot[17], cyc_bad[17]; static int minbalHC;
static int path[16], used, start;
static void dfs(int len) {   /* cycles with minimum vertex 'start', each directed cycle once */
    int last = path[len - 1];
    if (len >= 3 && (adj[last] >> start & 1)) {
        int nb = 0; for (int i = 0; i < len; i++) nb += bal[path[i]][path[(i + 1) % len]];
        cyc_tot[len]++; if (!nb) cyc_bad[len]++;
        if (len == n && nb < minbalHC) minbalHC = nb;
    }
    if (len == n) return;
    for (int u = start + 1; u < n; u++) if (!(used >> u & 1) && (adj[last] >> u & 1)) { used |= 1 << u; path[len] = u; dfs(len + 1); used &= ~(1 << u); }
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); int N = 1 << n, full = N - 1;
    f = calloc((size_t)N * n, sizeof(u64)); g = calloc((size_t)N * n, sizeof(u64));
    char line[256]; long long total = 0; long long minbal_hist[17] = {0};
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == '>' || line[0] == '\n') continue;
        memset(adj, 0, sizeof adj); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') adj[i] |= 1 << j; else adj[j] |= 1 << i; k++; }
        memset(f, 0, sizeof(u64) * (size_t)N * n); memset(g, 0, sizeof(u64) * (size_t)N * n);
        for (int v = 0; v < n; v++) { f[(1 << v) * n + v] = 1; g[(1 << v) * n + v] = 1; }
        for (int S = 1; S < N; S++) { if (!(S & (S - 1))) continue;
            for (int v = 0; v < n; v++) if (S >> v & 1) { int R = S ^ (1 << v); u64 a = 0, b = 0;
                for (int u = 0; u < n; u++) if (R >> u & 1) { if (adj[u] >> v & 1) a += f[R * n + u]; if (adj[v] >> u & 1) b += g[R * n + u]; }
                f[S * n + v] = a; g[S * n + v] = b; } }
        H = 0; for (int v = 0; v < n; v++) H += f[full * n + v];
        memset(cnt, 0, sizeof cnt);
        for (int x = 0; x < n; x++) for (int S = 0; S < N; S++) if (S >> x & 1) { u64 a = f[S * n + x]; if (!a) continue;
            int C = (full ^ S) | (1 << x); u64 prod = a * g[C * n + x]; if (!prod) continue;
            int rest = full ^ S; for (int y = 0; y < n; y++) if (rest >> y & 1) cnt[x][y] += prod; }
        total++;
        for (int x = 0; x < n; x++) for (int y = 0; y < n; y++) bal[x][y] = (x != y) && 3 * cnt[x][y] >= H && 3 * cnt[x][y] <= 2 * H;
        minbalHC = 99;
        for (start = 0; start < n; start++) { path[0] = start; used = 1 << start; dfs(1); }
        minbal_hist[minbalHC < 17 ? minbalHC : 16]++;
    }
    printf("n=%d strong %lld\n", n, total);
    for (int L = 3; L <= n; L++) printf("  %d-cycles: total %lld, with no balanced arc %lld\n", L, cyc_tot[L], cyc_bad[L]);
    printf("  min #balanced arcs on a Ham cycle, histogram over T:"); for (int b = 0; b < 17; b++) if (minbal_hist[b]) printf(" %d:%lld", b, minbal_hist[b]); printf("\n");
    return 0;
}
