/* Strong tournaments (gentourng -c ascii): candidate strengthenings of the tournament 1/3-2/3 conjecture.
   (C) every Hamiltonian cycle has a balanced arc  -- reports min #balanced arcs over Ham cycles (histogram) and failures;
   (D) the balanced graph (pairs with H/3 <= c <= 2H/3) is connected;
   (B) vertex-wise constant min_v max_w min(p, 1-p) (exact fraction) and how many tournaments attain it.
   Usage: gentourng -q -c n | ./tbal5 n */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n; static int adj[16];
static u64 *f, *g; static u64 cnt[16][16]; static u64 H; static int bal[16][16];
static int path[16], used, minbal; static long long ncyc;
static void dfs(int len) {   /* Hamiltonian cycles through vertex 0, each once */
    int last = path[len - 1];
    if (len == n) {
        if (adj[last] & 1) { ncyc++; int nb = 0; for (int i = 0; i < n; i++) nb += bal[path[i]][path[(i + 1) % n]]; if (nb < minbal) minbal = nb; }
        return;
    }
    for (int u = 1; u < n; u++) if (!(used >> u & 1) && (adj[last] >> u & 1)) { used |= 1 << u; path[len] = u; dfs(len + 1); used &= ~(1 << u); }
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); int N = 1 << n, full = N - 1;
    f = calloc((size_t)N * n, sizeof(u64)); g = calloc((size_t)N * n, sizeof(u64));
    char line[256]; long long total = 0, cfail = 0, dfail = 0, cyc_all = 0; long long mb_hist[20] = {0};
    u64 vw_num = 1, vw_den = 1; char vwline[256] = ""; int nvw = 0;
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
        minbal = 99; ncyc = 0; path[0] = 0; used = 1; dfs(1); cyc_all += ncyc;
        if (minbal == 0) { cfail++; if (cfail <= 5) printf("C FAILS: %s", line); }
        mb_hist[minbal < 19 ? minbal : 19]++;
        int seen = 1, ch = 1;
        while (ch) { ch = 0; for (int u = 0; u < n; u++) if (seen >> u & 1) for (int w = 0; w < n; w++) if (bal[u][w] && !(seen >> w & 1)) { seen |= 1 << w; ch = 1; } }
        if (seen != full) { dfail++; if (dfail <= 5) printf("D FAILS: %s", line); }
        u64 vmin_num = H, vmin_den = H;
        for (int v = 0; v < n; v++) { u64 bz = 0;
            for (int w = 0; w < n; w++) if (w != v) { u64 c = cnt[v][w], m = c < H - c ? c : H - c; if (m > bz) bz = m; }
            if ((unsigned __int128)bz * vmin_den < (unsigned __int128)vmin_num * H) { vmin_num = bz; vmin_den = H; } }
        if ((unsigned __int128)vmin_num * vw_den < (unsigned __int128)vw_num * vmin_den) { vw_num = vmin_num; vw_den = vmin_den; strcpy(vwline, line); nvw = 1; }
        else if ((unsigned __int128)vmin_num * vw_den == (unsigned __int128)vw_num * vmin_den) nvw++;
        if (total % 100000 == 0) { fprintf(stderr, "progress %lld\n", total); }
    }
    printf("n=%d strong %lld; Ham cycles %lld\n", n, total, cyc_all);
    printf("  (C) tournaments with a Ham cycle carrying no balanced arc: %lld\n", cfail);
    printf("  min #balanced arcs on a Ham cycle, histogram:"); for (int b = 0; b < 20; b++) if (mb_hist[b]) printf(" %d:%lld", b, mb_hist[b]); printf("\n");
    printf("  (D) tournaments with disconnected balanced graph: %lld\n", dfail);
    printf("  (B) vertex-wise constant = %llu/%llu = %.6f, attained by %d tournaments, e.g. %s", vw_num, vw_den, (double)vw_num / vw_den, nvw, vwline);
    return 0;
}
