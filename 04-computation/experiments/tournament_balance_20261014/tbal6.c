/* Strong tournaments (gentourng -c): local forms of vertex-wise balance (HYP-9279(B)).
   R_tri: every vertex w has a balanced partner x lying on a common directed 3-cycle with w;
   R_maj: every vertex w has a balanced partner x with p in [1/3,2/3] AND the arc between them forward (p(tail first) > 1/2);
   R_both: every vertex w has a partner x satisfying both (a forward balanced arc, 1/2 < p(tail first) <= 2/3, on a 3-cycle);
   also the minimum over vertices of the best balance among 3-cycle partners (exact). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n; static int adj[16];
static u64 *f, *g; static u64 cnt[16][16]; static u64 H;
int main(int argc, char **argv) {
    n = atoi(argv[1]); int N = 1 << n, full = N - 1;
    f = calloc((size_t)N * n, sizeof(u64)); g = calloc((size_t)N * n, sizeof(u64));
    char line[256]; long long total = 0, trifail = 0, majfail = 0, bothfail = 0; int ex = 0;
    u64 tb_num = 1, tb_den = 1; char tbline[256] = "";
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
        int tf = 0, mf = 0, bf = 0;
        for (int w = 0; w < n; w++) {
            int got_tri = 0, got_maj = 0, got_both = 0; u64 best = 0;
            for (int x = 0; x < n; x++) if (x != w) {
                u64 c = cnt[w][x]; int balanced = 3 * c >= H && 3 * c <= 2 * H;
                /* common 3-cycle: w->x->y->w or x->w->y->x for some y */
                int tri = 0;
                for (int y = 0; y < n && !tri; y++) if (y != w && y != x) {
                    if ((adj[w] >> x & 1) && (adj[x] >> y & 1) && (adj[y] >> w & 1)) tri = 1;
                    if ((adj[x] >> w & 1) && (adj[w] >> y & 1) && (adj[y] >> x & 1)) tri = 1;
                }
                if (tri) { u64 m = c < H - c ? c : H - c; if (m > best) best = m; if (balanced) got_tri = 1; }
                int tail = (adj[w] >> x & 1) ? w : x, head = tail == w ? x : w;
                if (balanced && 2 * cnt[tail][head] > H) { got_maj = 1; if (tri) got_both = 1; }
            }
            if (!got_tri) tf = 1; if (!got_maj) mf = 1; if (!got_both) bf = 1;
            if ((unsigned __int128)best * tb_den < (unsigned __int128)tb_num * H) { tb_num = best; tb_den = H; strcpy(tbline, line); }
        }
        if (tf) { trifail++; if (ex++ < 3) printf("R_tri fails: %s", line); }
        if (mf) majfail++; if (bf) bothfail++;
    }
    printf("n=%d strong %lld: R_tri failures %lld; R_maj failures %lld; R_both failures %lld; min over vertices of best 3-cycle-partner balance %llu/%llu = %.6f at %s",
           n, total, trifail, majfail, bothfail, tb_num, tb_den, (double)tb_num / tb_den, tbline);
    return 0;
}
