/* Exact balance statistics of uniform Hamiltonian paths for tournaments read from gentourng (default ascii format:
   upper triangle row by row, bit (i,j), i<j, is 1 iff i -> j).
   For each tournament: H = #Hamiltonian paths, c[x][y] = #paths with x before y (subset DP: forward f, backward g),
   delta = max_{x<y} min(c,H-c)/H.  Reports: count, #non-transitive with delta < 1/3, minimum delta (exact fraction) with
   witnesses, #with a balanced arc (H/3 <= c[x][y] <= 2H/3 for an arc x->y), minimum over arcs of p(tail first).
   Usage: gentourng [-c] n | ./tbal n [maxwitness] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
static int n; static int adj[16];
static u64 *f, *g;
static u64 cnt[16][16];
int main(int argc, char **argv) {
    n = atoi(argv[1]); int maxw = argc > 2 ? atoi(argv[2]) : 3;
    int N = 1 << n, full = N - 1;
    f = calloc((size_t)N * n, sizeof(u64)); g = calloc((size_t)N * n, sizeof(u64));
    char line[256];
    long long total = 0, nontrans = 0, below = 0, with_bal_arc = 0, noarc_examples = 0;
    long long with_cut = 0, every_cut_bal = 0, some_cut_bal = 0, every_vertex_bal = 0; int nocut_ex = 0;
    u64 best_num = 1, best_den = 1;   /* minimal delta among non-transitive, as fraction num/den (start at 1) */
    u64 arc_num = 1, arc_den = 1;     /* minimal p(tail first) over arcs */
    int nbest = 0; char bestline[8][256]; int nb_store = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == '>' || line[0] == '\n') continue;
        memset(adj, 0, sizeof adj);
        int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            if (line[k] == '1') adj[i] |= 1 << j; else adj[j] |= 1 << i;
            k++;
        }
        /* forward DP */
        memset(f, 0, sizeof(u64) * (size_t)N * n); memset(g, 0, sizeof(u64) * (size_t)N * n);
        for (int v = 0; v < n; v++) { f[(1 << v) * n + v] = 1; g[(1 << v) * n + v] = 1; }
        for (int S = 1; S < N; S++) {
            if (!(S & (S - 1))) continue;
            for (int v = 0; v < n; v++) if (S >> v & 1) {
                int R = S ^ (1 << v); u64 a = 0, b = 0;
                for (int u = 0; u < n; u++) if (R >> u & 1) {
                    if (adj[u] >> v & 1) a += f[R * n + u];      /* path ending at u, then u -> v */
                    if (adj[v] >> u & 1) b += g[R * n + u];      /* v -> u, then path starting at u */
                }
                f[S * n + v] = a; g[S * n + v] = b;
            }
        }
        u64 H = 0; for (int v = 0; v < n; v++) H += f[full * n + v];
        memset(cnt, 0, sizeof cnt);
        for (int x = 0; x < n; x++) {
            for (int S = 0; S < N; S++) if (S >> x & 1) {
                u64 a = f[S * n + x]; if (!a) continue;
                int C = (full ^ S) | (1 << x);
                u64 prod = a * g[C * n + x]; if (!prod) continue;
                int rest = full ^ S;
                for (int y = 0; y < n; y++) if (rest >> y & 1) cnt[x][y] += prod;
            }
        }
        total++;
        /* delta */
        u64 bm = 0; /* best min(c, H-c) */
        int has_bal_arc = 0; u64 amin = H; /* min over arcs of c[tail][head] */
        for (int x = 0; x < n; x++) for (int y = x + 1; y < n; y++) {
            u64 c = cnt[x][y]; u64 m = c < H - c ? c : H - c; if (m > bm) bm = m;
            int tail = (adj[x] >> y & 1) ? x : y, head = tail == x ? y : x;
            u64 ct = cnt[tail][head];
            if (3 * ct >= H && 3 * ct <= 2 * H) has_bal_arc = 1;
            if (ct < amin) amin = ct;
        }
        if (bm == 0) continue;            /* transitive: unique path, every pair determined */
        {   /* per-vertex: does z have a balanced partner? is T - z strong? */
            int ncut = 0, cut_bal = 0, all_bal = 1;
            for (int z = 0; z < n; z++) {
                u64 bz = 0;
                for (int x = 0; x < n; x++) if (x != z) { u64 c = cnt[z][x]; u64 m = c < H - c ? c : H - c; if (m > bz) bz = m; }
                int zbal = 3 * bz >= H;
                if (!zbal) all_bal = 0;
                int R = full ^ (1 << z); int s0 = __builtin_ctz(R);
                int fw = 1 << s0, bw = 1 << s0, ch = 1;
                while (ch) { ch = 0; for (int u = 0; u < n; u++) if (fw >> u & 1) { int nf = fw | (adj[u] & R); if (nf != fw) { fw = nf; ch = 1; } } }
                ch = 1;
                while (ch) { ch = 0; for (int u = 0; u < n; u++) if ((R >> u & 1) && !(bw >> u & 1)) { if (adj[u] & bw) { bw |= 1 << u; ch = 1; } } }
                int strongz = (fw == R) && (bw == R);
                if (!strongz) { ncut++; if (zbal) cut_bal++; }
            }
            if (ncut) { with_cut++; if (cut_bal == ncut) every_cut_bal++; if (cut_bal > 0) some_cut_bal++; else if (nocut_ex < 3) { nocut_ex++; printf("NO CUT VERTEX WITH BALANCED PARTNER: %s", line); } }
            if (all_bal) every_vertex_bal++;
        }
        nontrans++;
        if (3 * bm < H) below++;
        if (has_bal_arc) with_bal_arc++; else if (noarc_examples < 3) { noarc_examples++; printf("NO BALANCED ARC: %s", line); }
        /* compare bm/H with best_num/best_den */
        if ((unsigned __int128)bm * best_den < (unsigned __int128)best_num * H) {
            best_num = bm; best_den = H; nb_store = 0;
        }
        if ((unsigned __int128)bm * best_den == (unsigned __int128)best_num * H && nb_store < maxw) {
            strcpy(bestline[nb_store++], line);
        }
        if ((unsigned __int128)amin * arc_den < (unsigned __int128)arc_num * H) { arc_num = amin; arc_den = H; }
    }
    printf("n=%d: tournaments read %lld, non-transitive %lld, with delta < 1/3: %lld, with a balanced arc: %lld\n",
           n, total, nontrans, below, with_bal_arc);
    printf("  minimum delta = %llu/%llu = %.6f ; minimum over arcs of p(tail first) = %llu/%llu = %.6f\n",
           best_num, best_den, (double)best_num / best_den, arc_num, arc_den, (double)arc_num / arc_den);
    for (int i = 0; i < nb_store; i++) printf("  minimiser: %s", bestline[i]);
    printf("  strong T with a cut vertex z (T - z not strong): %lld; every cut vertex has a balanced partner: %lld; some cut vertex does: %lld; every vertex of T has a balanced partner: %lld\n",
           with_cut, every_cut_bal, some_cut_bal, every_vertex_bal);
    return 0;
}
