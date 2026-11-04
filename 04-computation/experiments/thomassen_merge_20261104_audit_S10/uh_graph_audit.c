/* Audit of THM-4870 (D) / HYP-9462 (B): labelled graphs  C_n (cycle 0-1-...-(n-1)-0) + chord set S  with exactly
 * one Hamiltonian cycle; their non-critical vertices (G - v Hamiltonian).
 * Written independently of the author's uh_graph_noncritical / uh_graph_downset programs.
 *   brute NLO NHI : every chord subset (2^(n(n-3)/2)), Hamiltonian cycles counted by DFS (stop at 2)
 *   down  NLO NHI : down-set enumeration (removing a chord never creates a Hamiltonian cycle), each extension
 *                   S -> S+{a,b} tested by "no Hamiltonian a-b path in G(S)"; OpenMP over depth-2 subtrees
 * For each UH graph: |NC|, whether two cycle-consecutive vertices are both non-critical, whether two
 * non-critical vertices are joined by a chord, #edges.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif
typedef uint32_t u32; typedef uint64_t u64;
#define MAXN 12
static int n, m, CA[80], CB[80];

/* count Hamiltonian cycles (directed, from vertex s0 = lowest vertex in vs), stop when count > lim */
static int hc_dfs(const u32 *adj, int s0, int end, u32 unv, int lim, int cnt) {
    if (!unv) return cnt + ((adj[end] >> s0) & 1);
    u32 c = adj[end] & unv;
    while (c) {
        int w = __builtin_ctz(c); c &= c - 1;
        cnt = hc_dfs(adj, s0, w, unv & ~(1u << w), lim, cnt);
        if (cnt > lim) return cnt;
    }
    return cnt;
}
static int count_hc(const u32 *adj, u32 vs, int lim) { /* number of DIRECTED cycles, capped */
    if (__builtin_popcount(vs) < 3) return 0;
    int s0 = __builtin_ctz(vs);
    return hc_dfs(adj, s0, s0, vs & ~(1u << s0), lim, 0);
}
static int hamiltonian_minus(const u32 *adj, int v) { /* is G - v Hamiltonian? */
    u32 a2[MAXN]; u32 vs = ((1u << n) - 1) & ~(1u << v);
    for (int i = 0; i < n; i++) a2[i] = adj[i] & ~(1u << v);
    for (int i = 0; i < n; i++) if (i != v && __builtin_popcount(a2[i]) < 2) return 0;
    return count_hc(a2, vs, 0) > 0;
}
/* Hamiltonian path from a to b covering all n vertices */
static int hp_dfs(const u32 *adj, int end, int b, u32 unv) {
    if (unv == (1u << b)) return (adj[end] >> b) & 1;
    u32 c = adj[end] & unv & ~(1u << b);
    while (c) {
        int w = __builtin_ctz(c); c &= c - 1;
        if (hp_dfs(adj, w, b, unv & ~(1u << w))) return 1;
    }
    return 0;
}

typedef struct { long long uh, consec, chordadj_graphs, nc_hist[MAXN + 1], wit_nc_hist[MAXN + 1]; int minedges_wit; u64 minwit_mask; } Stats;
static void stats_init(Stats *s) { memset(s, 0, sizeof(*s)); s->minedges_wit = 1 << 30; }
static void analyse(const u32 *adj, u64 mask, int nch, Stats *st) {
    int nc[MAXN], k = 0;
    for (int v = 0; v < n; v++) { nc[v] = hamiltonian_minus(adj, v); k += nc[v]; }
    st->uh++; st->nc_hist[k]++;
    int consec = 0;
    for (int v = 0; v < n; v++) if (nc[v] && nc[(v + 1) % n]) consec = 1;
    int chordadj = 0;
    for (int c = 0; c < m; c++) if ((mask >> c & 1) && nc[CA[c]] && nc[CB[c]]) chordadj = 1;
    if (chordadj) st->chordadj_graphs++;
    if (consec) {
        st->consec++; st->wit_nc_hist[k]++;
        if (n + nch < st->minedges_wit || (n + nch == st->minedges_wit && mask < st->minwit_mask)) { st->minedges_wit = n + nch; st->minwit_mask = mask; }
    }
}
static void build(u64 mask, u32 *adj) {
    for (int i = 0; i < n; i++) adj[i] = (1u << ((i + 1) % n)) | (1u << ((i + n - 1) % n));
    for (int c = 0; c < m; c++) if (mask >> c & 1) { adj[CA[c]] |= 1u << CB[c]; adj[CB[c]] |= 1u << CA[c]; }
}
static void setup(int nn) {
    n = nn; m = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        if (j == i + 1 || (i == 0 && j == n - 1)) continue;
        CA[m] = i; CB[m] = j; m++;
    }
}
static void report(const char *mode, Stats *st, double secs) {
    printf("[%s] n=%2d chords=%d  UH graphs=%lld  with two cycle-consecutive NC vertices=%lld  with two NC joined by a chord=%lld  (%.1fs)\n",
           mode, n, m, st->uh, st->consec, st->chordadj_graphs, secs);
    printf("      |NC| distribution:");
    for (int k = 0; k <= n; k++) if (st->nc_hist[k]) printf(" %d:%lld", k, st->nc_hist[k]);
    if (st->consec) {
        printf("\n      witnesses' |NC| distribution:");
        for (int k = 0; k <= n; k++) if (st->wit_nc_hist[k]) printf(" %d:%lld", k, st->wit_nc_hist[k]);
        printf("\n      fewest edges among witnesses: %d; lexicographically first such chord set:", st->minedges_wit);
        for (int c = 0; c < m; c++) if (st->minwit_mask >> c & 1) printf(" {%d,%d}", CA[c], CB[c]);
    }
    printf("\n");
    fflush(stdout);
}

static void brute(int nn) {
    setup(nn);
    Stats st; stats_init(&st);
#ifdef _OPENMP
    double t0 = omp_get_wtime();
#endif
    u64 total = (u64)1 << m;
    #pragma omp parallel
    {
        Stats loc; stats_init(&loc);
        #pragma omp for schedule(dynamic, 4096)
        for (long long mask = 0; mask < (long long)total; mask++) {
            u32 adj[MAXN]; build((u64)mask, adj);
            if (count_hc(adj, (1u << n) - 1, 2) != 2) continue; /* exactly one undirected cycle = 2 directed */
            analyse(adj, (u64)mask, __builtin_popcountll(mask), &loc);
        }
        #pragma omp critical
        {
            st.uh += loc.uh; st.consec += loc.consec; st.chordadj_graphs += loc.chordadj_graphs;
            for (int k = 0; k <= MAXN; k++) { st.nc_hist[k] += loc.nc_hist[k]; st.wit_nc_hist[k] += loc.wit_nc_hist[k]; }
            if (loc.minedges_wit < st.minedges_wit || (loc.minedges_wit == st.minedges_wit && loc.minwit_mask < st.minwit_mask)) { st.minedges_wit = loc.minedges_wit; st.minwit_mask = loc.minwit_mask; }
        }
    }
    double secs = 0;
#ifdef _OPENMP
    secs = omp_get_wtime() - t0;
#endif
    report("brute", &st, secs);
}

/* down-set DFS: S is UH; try adding chords with index > last */
static void down_rec(u64 mask, int last, int nch, Stats *st) {
    u32 adj[MAXN]; build(mask, adj);
    analyse(adj, mask, nch, st);
    for (int c = last + 1; c < m; c++) {
        /* G(S + c) has a 2nd Hamiltonian cycle iff G(S) has a Hamiltonian CA[c]-CB[c] path */
        if (hp_dfs(adj, CA[c], CB[c], ((1u << n) - 1) & ~(1u << CA[c]))) continue;
        down_rec(mask | ((u64)1 << c), c, nch + 1, st);
    }
}
static void down(int nn) {
    setup(nn);
    Stats st; stats_init(&st);
#ifdef _OPENMP
    double t0 = omp_get_wtime();
#endif
    /* level 0 and 1 serially (analyse only), level-2 subtrees in parallel */
    { u32 adj[MAXN]; build(0, adj); analyse(adj, 0, 0, &st); }
    int pa[4000], pb[4000], np = 0; int single_ok[80];
    { u32 adj[MAXN]; build(0, adj);
      for (int c = 0; c < m; c++) single_ok[c] = !hp_dfs(adj, CA[c], CB[c], ((1u << n) - 1) & ~(1u << CA[c])); }
    for (int c = 0; c < m; c++) if (single_ok[c]) {
        u32 adj[MAXN]; build((u64)1 << c, adj); analyse(adj, (u64)1 << c, 1, &st);
        for (int d = c + 1; d < m; d++) if (!hp_dfs(adj, CA[d], CB[d], ((1u << n) - 1) & ~(1u << CA[d]))) { pa[np] = c; pb[np] = d; np++; }
    }
    #pragma omp parallel
    {
        Stats loc; stats_init(&loc);
        #pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < np; i++) down_rec(((u64)1 << pa[i]) | ((u64)1 << pb[i]), pb[i], 2, &loc);
        #pragma omp critical
        {
            st.uh += loc.uh; st.consec += loc.consec; st.chordadj_graphs += loc.chordadj_graphs;
            for (int k = 0; k <= MAXN; k++) { st.nc_hist[k] += loc.nc_hist[k]; st.wit_nc_hist[k] += loc.wit_nc_hist[k]; }
            if (loc.minedges_wit < st.minedges_wit || (loc.minedges_wit == st.minedges_wit && loc.minwit_mask < st.minwit_mask)) { st.minedges_wit = loc.minedges_wit; st.minwit_mask = loc.minwit_mask; }
        }
    }
    double secs = 0;
#ifdef _OPENMP
    secs = omp_get_wtime() - t0;
#endif
    report("down ", &st, secs);
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: brute|down NLO NHI\n"); return 1; }
    int lo = atoi(argv[2]), hi = atoi(argv[3]);
    for (int k = lo; k <= hi; k++) {
        if (!strcmp(argv[1], "brute")) brute(k); else down(k);
    }
    return 0;
}
