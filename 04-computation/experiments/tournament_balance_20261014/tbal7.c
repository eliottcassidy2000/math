/* Cut vertices of strong tournaments (gentourng -c): z with T - z not strong; strong components D_1 => ... => D_k of T - z.
   Tests (exact):  I1: z is balanced against EVERY vertex of the interior blocks D_2 .. D_(k-1) (k >= 3);
                   I2: z has a balanced partner in D_1 u D_2, or in D_(k-1) u D_k;
                   I3: z has a balanced partner at all (HYP-9279(B) at cut vertices).
   Also the histogram of k and, for failures of I1, the block sizes. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n; static int adj[16];
static u64 *f, *g; static u64 cnt[16][16]; static u64 H;
int main(int argc, char **argv) {
    n = atoi(argv[1]); int N = 1 << n, full = N - 1;
    f = calloc((size_t)N * n, sizeof(u64)); g = calloc((size_t)N * n, sizeof(u64));
    char line[256]; long long total = 0, ncut = 0, i1fail = 0, i1tested = 0, i2fail = 0, i3fail = 0, khist[17] = {0};
    int ex = 0;
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
        for (int z = 0; z < n; z++) {
            int R = full ^ (1 << z);
            /* strong components of T - z, ordered: a tournament's condensation is a transitive order; component of v =
               vertices mutually reachable.  Score order: sort components by (number of vertices of R beaten) */
            int comp[16]; int nc = 0; int done = 0;
            for (int v = 0; v < n; v++) comp[v] = -1;
            int reachm[16];
            for (int v = 0; v < n; v++) if (R >> v & 1) {
                int seen = 1 << v, ch = 1;
                while (ch) { ch = 0; for (int u = 0; u < n; u++) if (seen >> u & 1) { int nf = seen | (adj[u] & R); if (nf != seen) { seen = nf; ch = 1; } } }
                reachm[v] = seen;
            }
            int compmask[16];
            for (int v = 0; v < n; v++) if ((R >> v & 1) && comp[v] < 0) {
                int m = 0; for (int u = 0; u < n; u++) if ((R >> u & 1) && (reachm[v] >> u & 1) && (reachm[u] >> v & 1)) m |= 1 << u;
                for (int u = 0; u < n; u++) if (m >> u & 1) comp[u] = nc;
                compmask[nc++] = m;
            }
            (void)done;
            if (nc < 2) continue;      /* T - z strong */
            ncut++;
            /* order components: component c1 before c2 iff some (any) vertex of c1 beats a vertex of c2 */
            int order[16]; for (int c = 0; c < nc; c++) order[c] = c;
            for (int a = 0; a < nc; a++) for (int b = a + 1; b < nc; b++) {
                int ca = order[a], cb = order[b]; int va = __builtin_ctz(compmask[ca]), vb = __builtin_ctz(compmask[cb]);
                if (adj[vb] >> va & 1) { order[a] = cb; order[b] = ca; }
            }
            khist[nc]++;
            int balz[16]; for (int x = 0; x < n; x++) balz[x] = x != z && 3 * cnt[z][x] >= H && 3 * cnt[z][x] <= 2 * H;
            int any = 0; for (int x = 0; x < n; x++) any |= balz[x];
            if (!any) i3fail++;
            if (nc >= 3) {
                i1tested++; int ok = 1;
                for (int c = 1; c < nc - 1; c++) for (int x = 0; x < n; x++) if ((compmask[order[c]] >> x & 1) && !balz[x]) ok = 0;
                if (!ok) { i1fail++; if (ex++ < 4) { printf("I1 fails: %s  z=%d blocks:", line, z); for (int c = 0; c < nc; c++) printf(" %d", __builtin_popcount(compmask[order[c]])); printf("\n"); } }
            }
            int m12 = compmask[order[0]] | compmask[order[1]], mlast = compmask[order[nc - 1]] | compmask[order[nc - 2]];
            int ok2 = 0; for (int x = 0; x < n; x++) if (balz[x] && (((m12 >> x) & 1) || ((mlast >> x) & 1))) ok2 = 1;
            if (!ok2) i2fail++;
        }
    }
    printf("n=%d strong %lld; cut-vertex instances %lld; k-histogram:", n, total, ncut);
    for (int c = 2; c <= n; c++) if (khist[c]) printf(" %d:%lld", c, khist[c]);
    printf("\n  I1 (z balanced vs all interior-block vertices, k>=3): tested %lld, failures %lld\n  I2 (partner in D1uD2 or D(k-1)uDk): failures %lld\n  I3 (some partner): failures %lld\n", i1tested, i1fail, i2fail, i3fail);
    return 0;
}
