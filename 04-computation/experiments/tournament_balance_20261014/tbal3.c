/* Exploration (strong tournaments from gentourng -c): structural strengthenings of vertex-wise balance.
   S1 every v has a balanced OUT-neighbour; S2 every v has a balanced IN-neighbour; S3 the balanced graph B is connected;
   S4 for every v the sorted values {0, p(v before w) (w != v), 1} have all gaps <= 1/3;
   S5 every Hamiltonian cycle has a balanced consecutive pair; histogram of the minimum degree of B;
   vertex-wise constant min_T min_v max_w min(p, 1-p). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n; static int adj[16];
static u64 *f, *g; static u64 cnt[16][16]; static u64 H;
static int bal[16][16];
static long long hc_total, hc_nobal; static int hc_found_nobal;
static int path[16], used;
static void dfs_cycle(int len) {
    int last = path[len - 1];
    if (len == n) {
        if (adj[last] >> 0 & 1) {
            hc_total++; int ok = 0;
            for (int i = 0; i < n; i++) { int a = path[i], b = path[(i + 1) % n]; if (bal[a][b]) { ok = 1; break; } }
            if (!ok) { hc_nobal++; hc_found_nobal = 1; }
        }
        return;
    }
    for (int u = 1; u < n; u++) if (!(used >> u & 1) && (adj[last] >> u & 1)) { used |= 1 << u; path[len] = u; dfs_cycle(len + 1); used &= ~(1 << u); }
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); int N = 1 << n, full = N - 1;
    f = calloc((size_t)N * n, sizeof(u64)); g = calloc((size_t)N * n, sizeof(u64));
    char line[256];
    long long total = 0, s1fail = 0, s2fail = 0, s3fail = 0, s4fail = 0, s5fail = 0, s12fail = 0, hc_all = 0, hc_bad = 0;
    u64 vw_num = 1, vw_den = 1; char vwline[256] = "";
    long long degmin_hist[16] = {0};
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
        int f1 = 0, f2 = 0, f12 = 0;
        for (int v = 0; v < n; v++) { int bo = 0, bi = 0;
            for (int w = 0; w < n; w++) if (bal[v][w]) { if (adj[v] >> w & 1) bo = 1; else bi = 1; }
            if (!bo) f1 = 1; if (!bi) f2 = 1; if (!bo && !bi) f12 = 1; }
        s1fail += f1; s2fail += f2; s12fail += f12;
        int seen = 1, ch = 1;
        while (ch) { ch = 0; for (int u = 0; u < n; u++) if (seen >> u & 1) for (int w = 0; w < n; w++) if (bal[u][w] && !(seen >> w & 1)) { seen |= 1 << w; ch = 1; } }
        if (seen != full) s3fail++;
        int dmin = 99; for (int v = 0; v < n; v++) { int d = 0; for (int w = 0; w < n; w++) d += bal[v][w]; if (d < dmin) dmin = d; }
        degmin_hist[dmin]++;
        int f4 = 0;
        for (int v = 0; v < n && !f4; v++) { u64 vals[18]; int m = 0; vals[m++] = 0; vals[m++] = H;
            for (int w = 0; w < n; w++) if (w != v) vals[m++] = cnt[v][w];
            for (int i = 1; i < m; i++) { u64 t = vals[i]; int j = i - 1; while (j >= 0 && vals[j] > t) { vals[j + 1] = vals[j]; j--; } vals[j + 1] = t; }
            for (int i = 0; i + 1 < m; i++) if (3 * (vals[i + 1] - vals[i]) > H) { f4 = 1; break; } }
        s4fail += f4;
        hc_total = 0; hc_nobal = 0; hc_found_nobal = 0; path[0] = 0; used = 1; dfs_cycle(1);
        hc_all += hc_total; hc_bad += hc_nobal; s5fail += hc_found_nobal;
        u64 vmin_num = H, vmin_den = H;
        for (int v = 0; v < n; v++) { u64 bz = 0; for (int w = 0; w < n; w++) if (w != v) { u64 c = cnt[v][w], m = c < H - c ? c : H - c; if (m > bz) bz = m; }
            if ((unsigned __int128)bz * vmin_den < (unsigned __int128)vmin_num * H) { vmin_num = bz; vmin_den = H; } }
        if ((unsigned __int128)vmin_num * vw_den < (unsigned __int128)vw_num * vmin_den) { vw_num = vmin_num; vw_den = vmin_den; strcpy(vwline, line); }
    }
    printf("n=%d strong tournaments %lld\n", n, total);
    printf("  S1 (some v w/o balanced out-nbr) %lld; S2 (some v w/o balanced in-nbr) %lld; both at one vertex %lld; S3 (B disconnected) %lld; S4 (gap > 1/3) %lld; S5 (Ham cycle w/o balanced arc) %lld (cycles %lld of %lld)\n",
           s1fail, s2fail, s12fail, s3fail, s4fail, s5fail, hc_bad, hc_all);
    printf("  min degree of B histogram:"); for (int d = 0; d < n; d++) if (degmin_hist[d]) printf(" %d:%lld", d, degmin_hist[d]); printf("\n");
    printf("  vertex-wise constant = %llu/%llu = %.6f at %s", vw_num, vw_den, (double)vw_num / vw_den, vwline);
    return 0;
}
