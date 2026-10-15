/* Minimiser tower for HYP-9279(A): read strong tournaments (gentourng ascii, n vertices) on stdin; for each, try all 2^n ways
   to attach a new vertex (keeping only strong results) and record the minimum delta over all extensions (exact).
   Prints the minimum and up to MAXOUT extension tournaments attaining it (ascii, n+1 vertices), for the next round.
   With a third argument thr, instead print every strong extension with delta <= thr (up to maxout), deduplicated by line.
   Usage: ./tbal_extend n maxout [thr] < in > out */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n; static int adj[16];
static u64 *f, *g; static u64 cnt[16][16];
static int strong(int m) {
    int full = (1 << m) - 1, seen = 1, ch = 1;
    while (ch) { ch = 0; for (int u = 0; u < m; u++) if (seen >> u & 1) { int nf = seen | (adj[u] & full); if (nf != seen) { seen = nf; ch = 1; } } }
    if (seen != full) return 0;
    seen = 1; ch = 1;
    while (ch) { ch = 0; for (int u = 0; u < m; u++) if (!(seen >> u & 1) && (adj[u] & seen)) { seen |= 1 << u; ch = 1; } }
    return seen == full;
}
static void delta(int m, u64 *num, u64 *den) {
    int N = 1 << m, full = N - 1;
    memset(f, 0, sizeof(u64) * (size_t)N * m); memset(g, 0, sizeof(u64) * (size_t)N * m);
    for (int v = 0; v < m; v++) { f[(1 << v) * m + v] = 1; g[(1 << v) * m + v] = 1; }
    for (int S = 1; S < N; S++) { if (!(S & (S - 1))) continue;
        for (int v = 0; v < m; v++) if (S >> v & 1) { int R = S ^ (1 << v); u64 a = 0, b = 0;
            for (int u = 0; u < m; u++) if (R >> u & 1) { if (adj[u] >> v & 1) a += f[R * m + u]; if (adj[v] >> u & 1) b += g[R * m + u]; }
            f[S * m + v] = a; g[S * m + v] = b; } }
    u64 H = 0; for (int v = 0; v < m; v++) H += f[full * m + v];
    memset(cnt, 0, sizeof cnt);
    for (int x = 0; x < m; x++) for (int S = 0; S < N; S++) if (S >> x & 1) { u64 a = f[S * m + x]; if (!a) continue;
        int C = (full ^ S) | (1 << x); u64 prod = a * g[C * m + x]; if (!prod) continue;
        int rest = full ^ S; for (int y = 0; y < m; y++) if (rest >> y & 1) cnt[x][y] += prod; }
    u64 bm = 0; for (int x = 0; x < m; x++) for (int y = x + 1; y < m; y++) { u64 c = cnt[x][y], mm = c < H - c ? c : H - c; if (mm > bm) bm = mm; }
    *num = bm; *den = H;
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); int maxout = atoi(argv[2]); int m = n + 1, N = 1 << m; double thr = argc > 3 ? atof(argv[3]) : -1; long nthr = 0;
    f = malloc(sizeof(u64) * (size_t)N * m); g = malloc(sizeof(u64) * (size_t)N * m);
    char line[256]; u64 bn = 1, bd = 1; char outs[64][256]; int nout = 0; long ext = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == '>' || line[0] == '\n' || line[0] == '#') continue;
        int base[16]; memset(base, 0, sizeof base); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') base[i] |= 1 << j; else base[j] |= 1 << i; k++; }
        for (int pat = 0; pat < (1 << n); pat++) {
            memcpy(adj, base, sizeof base); adj[n] = 0;
            for (int i = 0; i < n; i++) { if (pat >> i & 1) adj[i] |= 1 << n; else adj[n] |= 1 << i; }
            if (!strong(m)) continue;
            ext++;
            u64 num, den; delta(m, &num, &den);
            if (thr > 0 && (double)num / den <= thr && nthr < maxout) {
                nthr++; int q = 0; char o[256];
                for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) o[q++] = (adj[i] >> j & 1) ? '1' : '0';
                o[q] = 0; printf("%s\n", o);
            }
            if ((unsigned __int128)num * bd < (unsigned __int128)bn * den) { bn = num; bd = den; nout = 0; }
            if ((unsigned __int128)num * bd == (unsigned __int128)bn * den && nout < maxout) {
                char *o = outs[nout++]; int q = 0;
                for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) o[q++] = (adj[i] >> j & 1) ? '1' : '0';
                o[q] = 0;
            }
        }
    }
    fprintf(stderr, "n=%d -> %d: %ld strong extensions; minimum delta %llu/%llu = %.6f\n", n, m, ext, bn, bd, (double)bn / bd);
    if (thr <= 0) for (int i = 0; i < nout; i++) printf("%s\n", outs[i]);
    else fprintf(stderr, "  printed %ld extensions with delta <= %.4f\n", nthr, thr);
    return 0;
}
