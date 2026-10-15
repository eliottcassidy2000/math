/* Annealing search beyond the exhaustive range (n = 11..15) for strong tournaments minimising
   mode 0: the vertex-wise constant  min_v max_w min(p, 1-p)   (HYP-9279(B); exactly 1/3 for n = 5..9)
   mode 1: the global balance delta = max_{x,y} min(p, 1-p)     (HYP-9277 / HYP-9279(A))
   Exact counts via the subset DP (u64 is enough for n <= 15).  Start: transitive tournament with a few reversed arcs.
   Moves: reverse one arc (rejecting non-strong results).  Accept if the objective does not increase, or with
   probability exp(-(increase)/temp).  Usage: ./tbal_anneal n mode seed steps */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef unsigned long long u64;
static int n; static int adj[16];
static u64 *f, *g; static u64 cnt[16][16];
static u64 rng_s;
static u64 rnd(void) { rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }
static int strong(void) {
    int full = (1 << n) - 1, seen = 1, ch = 1;
    while (ch) { ch = 0; for (int u = 0; u < n; u++) if (seen >> u & 1) { int nf = seen | adj[u]; if (nf != seen) { seen = nf; ch = 1; } } }
    if (seen != full) return 0;
    seen = 1; ch = 1;
    while (ch) { ch = 0; for (int u = 0; u < n; u++) if (!(seen >> u & 1) && (adj[u] & seen)) { seen |= 1 << u; ch = 1; } }
    return seen == full;
}
static double evaluate(int mode, u64 *Hout) {
    int N = 1 << n, full = N - 1;
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
    *Hout = H;
    if (mode == 1) { u64 bm = 0; for (int x = 0; x < n; x++) for (int y = x + 1; y < n; y++) { u64 c = cnt[x][y], m = c < H - c ? c : H - c; if (m > bm) bm = m; } return (double)bm / H; }
    double vmin = 1;
    for (int v = 0; v < n; v++) { u64 bz = 0; for (int w = 0; w < n; w++) if (w != v) { u64 c = cnt[v][w], m = c < H - c ? c : H - c; if (m > bz) bz = m; }
        double r = (double)bz / H; if (r < vmin) vmin = r; }
    return vmin;
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); int mode = atoi(argv[2]); rng_s = 0x9E3779B97F4A7C15ULL ^ (u64)atoll(argv[3]); long steps = atol(argv[4]);
    int N = 1 << n; f = malloc(sizeof(u64) * (size_t)N * n); g = malloc(sizeof(u64) * (size_t)N * n);
    do {   /* start: transitive plus 1..4 reversed arcs */
        memset(adj, 0, sizeof adj);
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) adj[i] |= 1 << j;
        int r = 1 + rnd() % 4;
        for (int t = 0; t < r; t++) { int i = rnd() % n, j = rnd() % n; if (i == j) continue; if (adj[i] >> j & 1) { adj[i] &= ~(1 << j); adj[j] |= 1 << i; } }
        if (!strong()) { int i = 0, j = n - 1; adj[i] &= ~(1 << j); adj[j] |= 1 << i; }
    } while (!strong());
    u64 H; double cur = evaluate(mode, &H), best = cur; int bestadj[16]; memcpy(bestadj, adj, sizeof adj);
    for (long s = 0; s < steps; s++) {
        double temp = 0.02 * (1.0 - (double)s / steps) + 1e-6;
        int i = rnd() % n, j = rnd() % n; if (i == j) continue;
        int fwd = adj[i] >> j & 1;
        if (fwd) { adj[i] &= ~(1 << j); adj[j] |= 1 << i; } else { adj[j] &= ~(1 << i); adj[i] |= 1 << j; }
        if (!strong()) { if (fwd) { adj[j] &= ~(1 << i); adj[i] |= 1 << j; } else { adj[i] &= ~(1 << j); adj[j] |= 1 << i; } continue; }
        double val = evaluate(mode, &H);
        if (val <= cur || (double)(rnd() % 1000000) / 1e6 < exp(-(val - cur) / temp)) {
            cur = val; if (val < best) { best = val; memcpy(bestadj, adj, sizeof adj); }
        } else { if (fwd) { adj[j] &= ~(1 << i); adj[i] |= 1 << j; } else { adj[i] &= ~(1 << j); adj[j] |= 1 << i; } }
    }
    memcpy(adj, bestadj, sizeof adj); double v = evaluate(mode, &H);
    printf("n=%d mode=%s seed=%s steps=%ld: best %.6f (H=%llu)  T=", n, mode ? "delta" : "vertexwise", argv[3], steps, v, H);
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) putchar((adj[i] >> j & 1) ? '1' : '0');
    printf("\n"); return 0;
}
