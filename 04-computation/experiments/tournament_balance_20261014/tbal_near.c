/* Print every strong tournament (gentourng -c ascii on stdin) with delta <= thr.  Usage: ./tbal_near n thr */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
int main(int argc, char **argv) {
    int n = atoi(argv[1]); double thr = atof(argv[2]); int N = 1 << n, full = N - 1;
    u64 *f = calloc((size_t)N * n, sizeof(u64)), *g = calloc((size_t)N * n, sizeof(u64)); u64 cnt[16][16]; int adj[16];
    char line[256]; long hits = 0;
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
        u64 H = 0; for (int v = 0; v < n; v++) H += f[full * n + v];
        memset(cnt, 0, sizeof cnt);
        for (int x = 0; x < n; x++) for (int S = 0; S < N; S++) if (S >> x & 1) { u64 a = f[S * n + x]; if (!a) continue;
            int C = (full ^ S) | (1 << x); u64 prod = a * g[C * n + x]; if (!prod) continue;
            int rest = full ^ S; for (int y = 0; y < n; y++) if (rest >> y & 1) cnt[x][y] += prod; }
        u64 bm = 0; for (int x = 0; x < n; x++) for (int y = x + 1; y < n; y++) { u64 c = cnt[x][y], m = c < H - c ? c : H - c; if (m > bm) bm = m; }
        if ((double)bm / H <= thr) { hits++; fputs(line, stdout); }
    }
    fprintf(stderr, "n=%d: %ld strong tournaments with delta <= %.4f\n", n, hits, thr);
    return 0;
}
