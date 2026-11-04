/* gen1.c (klein-2026-11-04) -- pruned one-vertex extension of digit-1-ready strong tournaments.
   stdin: lines of n out-row encodings (order n).  For each sole-exit vertex a -> r of digit 1, form Y = Phi(T) (core conversed,
   a re-attached as sole exit), find the unique Hamiltonian path of Y - a from r (the unique cycle of Y), and try every pattern
   for a new vertex w whose "beats w" bits along that cycle form one cyclic block not containing a (necessary for Y + w to keep
   a unique cycle), converted back to T; plus the n patterns where w is a new sole exit.  Output "H digraph6" for every
   extension that is strong and digit-1-ready. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n; static u64 adj[64]; static u64 *f;
static u64 paths_end(u64 verts, int end) {   /* Hamiltonian paths of T[verts]; if end >= 0 only those ending at end */
    int m = __builtin_popcountll(verts); int map[64], k = 0;
    for (int v = 0; v < 64 && k < m; v++) if (verts >> v & 1) map[k++] = v;
    u64 N = 1ULL << m; memset(f, 0, sizeof(u64) * N * m);
    for (int i = 0; i < m; i++) f[(1ULL << i) * m + i] = 1;
    for (u64 S = 1; S < N; S++) for (int i = 0; i < m; i++) { u64 c = f[S * m + i]; if (!c) continue; u64 out = adj[map[i]];
        for (int j = 0; j < m; j++) if (!(S >> j & 1) && (out >> map[j] & 1)) f[(S | (1ULL << j)) * m + j] += c; }
    u64 h = 0; for (int i = 0; i < m; i++) if (end < 0 || map[i] == end) h += f[(N - 1) * m + i]; return h;
}
static int is_strong(void) {
    u64 all = (1ULL << n) - 1, seen = 1, fr = 1;
    while (fr) { u64 nf = 0; for (int v = 0; v < n; v++) if (fr >> v & 1) nf |= adj[v]; nf &= ~seen; seen |= nf; fr = nf; }
    if (seen != all) return 0;
    u64 radj[64]; for (int v = 0; v < n; v++) radj[v] = 0;
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (adj[u] >> v & 1) radj[v] |= 1ULL << u;
    seen = 1; fr = 1;
    while (fr) { u64 nf = 0; for (int v = 0; v < n; v++) if (fr >> v & 1) nf |= radj[v]; nf &= ~seen; seen |= nf; fr = nf; }
    return seen == all;
}
/* unique Hamiltonian path of the core Yc (vertices 'verts', adjacency 'yadj') starting at r: reconstruct into order[] (length m) */
static int unique_path(u64 *yadj, u64 verts, int r, int *order) {
    int m = __builtin_popcountll(verts); int map[64], k = 0, inv[64];
    for (int v = 0; v < 64 && k < m; v++) if (verts >> v & 1) { inv[v] = k; map[k++] = v; }
    u64 N = 1ULL << m; memset(f, 0, sizeof(u64) * N * m);
    f[(1ULL << inv[r]) * m + inv[r]] = 1;
    for (u64 S = 1; S < N; S++) for (int i = 0; i < m; i++) { u64 c = f[S * m + i]; if (!c) continue; u64 out = yadj[map[i]];
        for (int j = 0; j < m; j++) if (!(S >> j & 1) && (out >> map[j] & 1)) f[(S | (1ULL << j)) * m + j] += c; }
    int endi = -1; u64 tot = 0; for (int i = 0; i < m; i++) { tot += f[(N - 1) * m + i]; if (f[(N - 1) * m + i]) endi = i; }
    if (tot != 1) return 0;
    u64 S = N - 1; int cur = endi; order[m - 1] = map[cur];
    for (int pos = m - 1; pos > 0; pos--) { u64 Sp = S & ~(1ULL << cur); int found = -1;
        for (int i = 0; i < m; i++) if ((Sp >> i & 1) && f[Sp * m + i] && (yadj[map[i]] >> map[cur] & 1)) { found = i; break; }
        cur = found; S = Sp; order[pos - 1] = map[cur]; }
    return 1;
}
static void emit(u64 h) {
    printf("%llu &%c", h, (char)(63 + n)); int acc = 0, cnt = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { acc = acc * 2 + ((adj[i] >> j) & 1); cnt++; if (cnt == 6) { putchar(63 + acc); acc = 0; cnt = 0; } }
    if (cnt) { while (cnt < 6) { acc *= 2; cnt++; } putchar(63 + acc); } putchar('\n');
}
static void test_and_emit(u64 *base, int n0, u64 wbeats /* bits over old vertices: v beats w */) {
    n = n0 + 1; memcpy(adj, base, sizeof(u64) * n0); adj[n0] = 0;
    for (int v = 0; v < n0; v++) { if (wbeats >> v & 1) adj[v] |= 1ULL << n0; else adj[n0] |= 1ULL << v; }
    if (!is_strong()) return;
    u64 all = (1ULL << n) - 1; int ready = 0;
    for (int a = 0; a < n && !ready; a++) if (__builtin_popcountll(adj[a]) == 1) { int r = __builtin_ctzll(adj[a]); if (paths_end(all & ~(1ULL << a), r) == 1) ready = 1; }
    if (ready) emit(paths_end(all, -1));
}
int main(int argc, char **argv) {
    int n0 = atoi(argv[1]); f = malloc(sizeof(u64) * (1ULL << (n0 + 1)) * (n0 + 1));
    char line[4096]; u64 base[64], yadj[64]; int order[64];
    while (fgets(line, sizeof line, stdin)) {
        int k = 0; char *p = line; while (k < n0) { char *e; u64 v = strtoull(p, &e, 10); if (e == p) break; base[k++] = v; p = e; }
        if (k != n0) continue;
        n = n0; memcpy(adj, base, sizeof(u64) * n0);
        u64 all0 = (1ULL << n0) - 1;
        /* collect patterns in a set to avoid duplicates (bitmask over old vertices) */
        static unsigned char seenpat[1 << 20]; memset(seenpat, 0, 1 << n0 > (1 << 20) ? (1 << 20) : (1 << n0));
        for (int a = 0; a < n0; a++) {
            if (__builtin_popcountll(base[a]) != 1) continue; int r = __builtin_ctzll(base[a]);
            n = n0; memcpy(adj, base, sizeof(u64) * n0);
            if (paths_end(all0 & ~(1ULL << a), r) != 1) continue;
            /* Y core: converse of T - a */
            for (int u = 0; u < n0; u++) { yadj[u] = 0; } 
            for (int u = 0; u < n0; u++) for (int v = 0; v < n0; v++) if (u != a && v != a && (base[v] >> u & 1)) yadj[u] |= 1ULL << v;   /* v -> u in T  =>  u -> v in Y */
            int m = n0 - 1; if (!unique_path(yadj, all0 & ~(1ULL << a), r, order)) continue;   /* order[0] = r ... order[m-1] */
            /* cycle of Y: a, order[0..m-1]; w pattern: block of positions [p, q] within 1..m (cycle positions 1..m = order[0..m-1]) where c_i beats w in Y */
            for (int pp = 0; pp < m; pp++) for (int qq = pp; qq < m; qq++) {
                u64 wbeats = 0;   /* in T: old vertex v beats w */
                for (int i = 0; i < m; i++) { int u = order[i]; int ybeats = (i >= pp && i <= qq); if (!ybeats) wbeats |= 1ULL << u; }   /* converse: u beats w in T iff w beats u in Y */
                /* a: w -> a in T (a does not beat w): bit 0 */
                if (wbeats < (1ULL << 20) && seenpat[wbeats]) continue; if (wbeats < (1ULL << 20)) seenpat[wbeats] = 1;
                test_and_emit(base, n0, wbeats);
            }
        }
        /* w as a new sole exit: w beats exactly one old vertex u */
        for (int u = 0; u < n0; u++) { u64 wbeats = all0 & ~(1ULL << u); if (wbeats < (1ULL << 20) && seenpat[wbeats]) continue; if (wbeats < (1ULL << 20)) seenpat[wbeats] = 1; test_and_emit(base, n0, wbeats); }
    }
    return 0;
}
