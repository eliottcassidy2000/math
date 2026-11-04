/* Audit of THM-4870 (C): square-sum graphs Q_N on {1..N}, i ~ j iff i + j is a perfect square (i != j).
 * Written independently of the author's square_sum_vertex_deleted_20261104.c.
 *   - Hamiltonian cycle / path existence by DFS with sound degree pruning;
 *   - for graphs on <= 25 vertices, Hamiltonian path existence ALSO by a subset DP (cross-check);
 *   - exact Hamiltonian cycle counts (DFS, each undirected cycle found twice).
 * Usage: sq NLO NHI
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64;
#define MAXV 64
static u64 adj[MAXV + 1]; /* bit j set iff i ~ j, vertices 1..N */
static int N;
static int issq(int x) { int r = 0; while ((r + 1) * (r + 1) <= x) r++; return r * r == x; }
static void build(int n) {
    N = n;
    for (int i = 0; i <= MAXV; i++) adj[i] = 0;
    for (int i = 1; i <= n; i++) for (int j = 1; j <= n; j++) if (i != j && issq(i + j)) adj[i] |= (u64)1 << j;
}
static int pc(u64 x) { return __builtin_popcountll(x); }

/* ---------- Hamiltonian cycle: DFS from s, vertex set 'vs' ---------- */
static u64 VS; static int S0; static long long cyc_count; static int stop_at_first; static int path_[MAXV + 1], best_cycle[MAXV + 1], plen;
static int ok_cycle_prune(u64 unv, int end) {
    u64 t = unv;
    u64 avail_base = unv | ((u64)1 << end) | ((u64)1 << S0);
    while (t) {
        int u = __builtin_ctzll(t); t &= t - 1;
        if (pc(adj[u] & avail_base) < 2) return 0;
    }
    return 1;
}
static void dfs_cycle(int end, u64 unv) {
    if (stop_at_first && cyc_count) return;
    if (!unv) {
        if (adj[end] >> S0 & 1) { if (!cyc_count) memcpy(best_cycle, path_, sizeof(int) * plen); cyc_count++; }
        return;
    }
    if (!ok_cycle_prune(unv, end)) return;
    u64 c = adj[end] & unv;
    while (c) {
        int w = __builtin_ctzll(c); c &= c - 1;
        path_[plen++] = w;
        dfs_cycle(w, unv & ~((u64)1 << w));
        plen--;
        if (stop_at_first && cyc_count) return;
    }
}
/* returns number of directed Ham cycles from the min-degree vertex (2 x undirected), or 0/1 if stop */
static long long ham_cycles(u64 vs, int first_only) {
    VS = vs; cyc_count = 0; stop_at_first = first_only;
    int nv = pc(vs);
    if (nv < 3) return 0;
    int s = -1, md = 1000;
    u64 t = vs;
    while (t) { int u = __builtin_ctzll(t); t &= t - 1; int d = pc(adj[u] & vs); if (d < 2) return 0; if (d < md) { md = d; s = u; } }
    S0 = s; plen = 0; path_[plen++] = s;
    dfs_cycle(s, vs & ~((u64)1 << s));
    return cyc_count;
}

/* ---------- Hamiltonian path existence by DFS ---------- */
static int found_path;
static void dfs_path(int end, u64 unv) {
    if (found_path) return;
    if (!unv) { found_path = 1; return; }
    /* prune: every unvisited u needs >= 1 neighbour in unv+end; those with exactly 1 must be last: at most 1 */
    u64 t = unv; int ones = 0; u64 base = unv | ((u64)1 << end);
    while (t) {
        int u = __builtin_ctzll(t); t &= t - 1;
        int d = pc(adj[u] & base & ~((u64)1 << u));
        if (d == 0) return;
        if (d == 1) { ones++; if (ones > 1) return; }
    }
    u64 c = adj[end] & unv;
    while (c && !found_path) { int w = __builtin_ctzll(c); c &= c - 1; dfs_path(w, unv & ~((u64)1 << w)); }
}
static int ham_path_dfs(u64 vs) {
    int nv = pc(vs);
    if (nv <= 1) return 1;
    /* vertices of degree 0 -> no; > 2 leaves -> no */
    u64 t = vs; int leaves = 0; int leaf = -1;
    while (t) { int u = __builtin_ctzll(t); t &= t - 1; int d = pc(adj[u] & vs); if (d == 0) return 0; if (d == 1) { leaves++; leaf = u; } }
    if (leaves > 2) return 0;
    found_path = 0;
    if (leaves >= 1) { dfs_path(leaf, vs & ~((u64)1 << leaf)); return found_path; }
    t = vs;
    while (t && !found_path) { int u = __builtin_ctzll(t); t &= t - 1; dfs_path(u, vs & ~((u64)1 << u)); }
    return found_path;
}

/* ---------- Hamiltonian path existence by subset DP (<= 25 vertices) ---------- */
static uint32_t *DPm;
static int ham_path_dp(u64 vs) {
    int idx[MAXV + 1], nv = 0, lab[32];
    for (int u = 1; u <= N; u++) if (vs >> u & 1) { idx[u] = nv; lab[nv++] = u; }
    if (nv <= 1) return 1;
    uint32_t a[32];
    for (int i = 0; i < nv; i++) { a[i] = 0; for (int j = 0; j < nv; j++) if (adj[lab[i]] >> lab[j] & 1) a[i] |= 1u << j; }
    size_t M = (size_t)1 << nv;
    memset(DPm, 0, sizeof(uint32_t) * M);
    for (int i = 0; i < nv; i++) DPm[(size_t)1 << i] = 1u << i;
    for (size_t m = 1; m < M; m++) {
        uint32_t e = DPm[m]; if (!e) continue;
        uint32_t nb = 0, t = e;
        while (t) { int i = __builtin_ctz(t); t &= t - 1; nb |= a[i]; }
        nb &= ~(uint32_t)m;
        while (nb) { int j = __builtin_ctz(nb); nb &= nb - 1; DPm[m | ((size_t)1 << j)] |= 1u << j; }
    }
    return DPm[M - 1] != 0;
}

int main(int argc, char **argv) {
    int lo = atoi(argv[1]), hi = atoi(argv[2]);
    DPm = malloc(sizeof(uint32_t) << 25);
    for (int n = lo; n <= hi; n++) {
        build(n);
        u64 all = 0; for (int i = 1; i <= n; i++) all |= (u64)1 << i;
        int md = 1000, nleaf = 0, niso = 0;
        for (int i = 1; i <= n; i++) { int d = pc(adj[i]); if (d < md) md = d; if (d == 1) nleaf++; if (d == 0) niso++; }
        int tr = ham_path_dfs(all);
        if (n <= 25 && tr != ham_path_dp(all)) printf("  !! DFS/DP disagree on traceability of Q_%d\n", n);
        long long hc2 = ham_cycles(all, 0);
        long long hcu = hc2 / 2;
        int cyc[MAXV + 1]; int clen = n;
        if (hc2) memcpy(cyc, best_cycle, sizeof(int) * n);
        /* vertex deletions */
        int nt = 0, nh = 0, dpdis = 0; char tl[512] = "", hl[512] = "";
        u64 hamset = 0;
        for (int v = 1; v <= n; v++) {
            u64 vs = all & ~((u64)1 << v);
            int t = ham_path_dfs(vs);
            if (n - 1 <= 25 && t != ham_path_dp(vs)) dpdis++;
            int h = ham_cycles(vs, 1) > 0;
            if (t) { nt++; sprintf(tl + strlen(tl), "%d ", v); }
            if (h) { nh++; hamset |= (u64)1 << v; sprintf(hl + strlen(hl), "%d ", v); }
        }
        printf("N=%2d mindeg=%d leaves=%d isolated=%d traceable=%d Hcycles=%lld | #v: Q-v traceable=%d%s | #v: Q-v Hamiltonian=%d {%s}%s\n",
               n, md, nleaf, niso, tr, hcu, nt, (n - 1 <= 25) ? (dpdis ? " (DP DISAGREES)" : " (DP agrees)") : "", nh, hl,
               "");
        if (!tr) printf("        traceable deletions: {%s}\n", tl);
        if (hcu == 1) {
            /* non-critical set: independence in Q_N and non-consecutiveness on the unique cycle */
            int indep = 1, consec = 0;
            for (int a = 1; a <= n; a++) for (int b = a + 1; b <= n; b++)
                if ((hamset >> a & 1) && (hamset >> b & 1) && (adj[a] >> b & 1)) indep = 0;
            for (int i = 0; i < clen; i++) { int a = cyc[i], b = cyc[(i + 1) % clen]; if ((hamset >> a & 1) && (hamset >> b & 1)) consec++; }
            printf("        unique cycle:");
            for (int i = 0; i < clen; i++) printf(" %d", cyc[i]);
            printf("\n        non-critical set size %d, independent in Q_N: %s, cycle-consecutive NC pairs: %d\n",
                   nh, indep ? "yes" : "NO", consec);
        }
        fflush(stdout);
    }
    return 0;
}
