/* opus-2026-10-31-S69: ninja pairs (= weak B pairs) of the square-sum graph, exhaustive for small n.

   Q_N: vertices 1..N, u ~ v iff u + v is a perfect square.
   A ninja pair at n: a Hamiltonian path P of Q_n (the Shorter) and P' of Q_(n+1) (the Longer), both starting 1,3,
   with equal one-based positional parity on every common vertex 1..n; the even-order member ends at 8 and the
   odd-order member ends at 10 (session124 / square_sum_B_base_interval definition).

   Modes
     pairs LO HI        : for each n in [LO,HI] count (Shorter, Longer) parity-compatible pairs exhaustively
                          (all Shorter paths, all Longer paths, intersect their parity masks on 1..n)
     show N A B END [K] : print up to K Hamiltonian paths of Q_N with prefix A,B (B = 0: none) ending at END (-1 any)
     cycles N K         : print up to K Hamiltonian cycles of Q_N per pair of neighbours of 1
     pairdump n K       : print up to K compatible (Shorter, Longer) path pairs at n, and their total number
     exists LO HI       : per N, does a Hamiltonian cycle through 8-1-3 (N even) / 15-1-3 (N odd) exist (early stop)
     paths LO HI        : census per N: Hamiltonian paths (undirected), Hamiltonian cycles, cycles through 8-1-3 /
                          15-1-3, paths starting 1,3 that end at 8 / 10 / 15
   Output is plain text, one line per n.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned long long u64;
#define MAXN 64

static int N;                 /* current graph order */
static int deg[MAXN + 1], adj[MAXN + 1][16];
static int path[MAXN + 1];
static int pos_of[MAXN + 1];

static int is_square(int x) { int r = 0; while ((r + 1) * (r + 1) <= x) r++; return r * r == x; }

static void build(int n) {
    N = n;
    for (int v = 1; v <= n; v++) {
        deg[v] = 0;
        for (int w = 1; w <= n; w++)
            if (w != v && is_square(v + w)) adj[v][deg[v]++] = w;
    }
}

/* pruning: every unvisited vertex needs >= 1 available neighbour (unvisited or the current end);
   at most one unvisited vertex (other than a forced end) may have exactly 1 available neighbour */
static int feasible(u64 vis, int cur, int endv) {
    int deadends = 0;
    for (int v = 1; v <= N; v++) {
        if (vis >> v & 1ULL) continue;
        int a = 0;
        for (int k = 0; k < deg[v]; k++) {
            int w = adj[v][k];
            if (!(vis >> w & 1ULL) || w == cur) a++;
        }
        if (a == 0) return 0;
        if (a == 1 && v != endv) { deadends++; if (deadends > 1) return 0; }
    }
    return 1;
}

/* ---- generic DFS with a callback at complete paths ---- */
typedef void (*leaf_fn)(int len);
static leaf_fn on_leaf;
static int want_end;          /* -1: any end */
static long long stop_after;  /* stop once this many leaves seen (-1: never) */
static long long leaves;

static void dfs(u64 vis, int len) {
    if (stop_after >= 0 && leaves >= stop_after) return;
    int cur = path[len - 1];
    if (len == N) {
        if (want_end < 0 || cur == want_end) { leaves++; on_leaf(len); }
        return;
    }
    if (cur == want_end) return;           /* the end vertex may only be last */
    if (!feasible(vis, cur, want_end)) return;
    for (int k = 0; k < deg[cur]; k++) {
        int w = adj[cur][k];
        if (vis >> w & 1ULL) continue;
        path[len] = w;
        dfs(vis | (1ULL << w), len + 1);
    }
}

/* ---- pairs mode ---- */
static u64 *masks; static long long nmask, capmask;
static int mask_upto;
static u64 parity_mask(int len) {
    u64 m = 0;
    for (int j = 0; j < len; j++) if (path[j] <= mask_upto && (j % 2 == 0)) m |= 1ULL << path[j];  /* position j+1 odd */
    return m;
}
static void store_mask(int len) {
    if (nmask == capmask) { capmask = capmask ? 2 * capmask : 1024; masks = realloc(masks, capmask * sizeof(u64)); }
    masks[nmask++] = parity_mask(len);
}
static int cmp64(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }
static long long pair_hits; static u64 example_mask;
static void probe_mask(int len) {
    u64 m = parity_mask(len);
    if (bsearch(&m, masks, nmask, sizeof(u64), cmp64)) { if (!pair_hits) example_mask = m; pair_hits++; }
}

static long long enumerate_from13(int order, int endv, leaf_fn f) {
    build(order);
    path[0] = 1; path[1] = 3;
    on_leaf = f; want_end = endv; stop_after = -1; leaves = 0;
    if (order < 3) return 0;
    dfs((1ULL << 1) | (1ULL << 3), 2);
    return leaves;
}

static void mode_pairs(int lo, int hi) {
    for (int n = lo; n <= hi; n++) {
        int endS = (n % 2 == 0) ? 8 : 10;          /* Shorter has order n */
        int endL = ((n + 1) % 2 == 0) ? 8 : 10;    /* Longer has order n+1 */
        nmask = 0; mask_upto = n;
        long long nL = enumerate_from13(n + 1, endL, store_mask);
        qsort(masks, nmask, sizeof(u64), cmp64);
        long long uniq = 0;
        for (long long i = 0; i < nmask; i++) if (i == 0 || masks[i] != masks[i - 1]) masks[uniq++] = masks[i];
        nmask = uniq;
        pair_hits = 0;
        long long nS = enumerate_from13(n, endS, probe_mask);
        printf("n=%d shorter(order %d, end %d) paths=%lld longer(order %d, end %d) paths=%lld distinct_longer_masks=%lld "
               "compatible_pairs(shorter paths matching some longer mask)=%lld NINJA=%s\n",
               n, n, endS, nS, n + 1, endL, nL, nmask, pair_hits, pair_hits ? "YES" : "no");
        fflush(stdout);
    }
}

/* ---- census mode ---- */
static long long cnt;
static void count_leaf(int len) { (void)len; cnt++; }
static long long count_paths_from(int order, int a, int b, int endv) {
    build(order);
    path[0] = a; int len = 1; u64 vis = 1ULL << a;
    if (b > 0) { path[1] = b; len = 2; vis |= 1ULL << b; }
    on_leaf = count_leaf; want_end = endv; stop_after = -1; leaves = 0; cnt = 0;
    dfs(vis, len);
    return cnt;
}
static void mode_paths(int lo, int hi) {
    for (int n = lo; n <= hi; n++) {
        build(n);
        /* all Hamiltonian paths (directed count from every start, then /2) */
        long long directed = 0;
        for (int s = 1; s <= n; s++) directed += count_paths_from(n, s, 0, -1);
        /* Hamiltonian cycles: paths from 1 whose last vertex is a neighbour of 1, each cycle counted twice */
        long long cyc2 = 0;
        build(n);
        for (int k = 0; k < deg[1]; k++) cyc2 += count_paths_from(n, 1, 0, adj[1][k]);
        long long c813 = (n >= 8) ? count_paths_from(n, 1, 3, 8) : 0;     /* cycles through 8-1-3 */
        long long c1513 = (n >= 15) ? count_paths_from(n, 1, 3, 15) : 0;  /* cycles through 15-1-3 */
        long long e10 = (n >= 10) ? count_paths_from(n, 1, 3, 10) : 0;
        printf("N=%d ham_paths=%lld ham_cycles=%lld cycles_through_8-1-3=%lld cycles_through_15-1-3=%lld paths_13_to_10=%lld\n",
               n, directed / 2, cyc2 / 2, c813, c1513, e10);
        fflush(stdout);
    }
}


/* ---- show mode: print paths with a given prefix and end (cycles: prefix 1 alone, end a neighbour of 1) ---- */
static int show_limit;
static void print_leaf(int len) {
    if (leaves > show_limit) return;
    printf("  ");
    for (int j = 0; j < len; j++) printf("%d%s", path[j], j + 1 < len ? "," : "\n");
}
static void mode_show(int order, int a, int b, int endv, int limit) {
    build(order);
    path[0] = a; int len = 1; u64 vis = 1ULL << a;
    if (b > 0) { path[1] = b; len = 2; vis |= 1ULL << b; }
    on_leaf = print_leaf; want_end = endv; stop_after = -1; leaves = 0; show_limit = limit;
    dfs(vis, len);
    printf("total=%lld\n", leaves);
}
static void mode_cycles(int order, int limit) {
    build(order);
    int nb[16], nd = deg[1];
    memcpy(nb, adj[1], sizeof(int) * nd);
    for (int k = 0; k < nd; k++) for (int k2 = 0; k2 < nd; k2++) {
        if (nb[k] >= nb[k2]) continue;          /* cycle through nb[k]-1-nb[k2], canonical orientation */
        build(order);
        path[0] = 1; path[1] = nb[k]; u64 vis = (1ULL << 1) | (1ULL << nb[k]);
        on_leaf = print_leaf; want_end = nb[k2]; stop_after = -1; leaves = 0; show_limit = limit;
        dfs(vis, 2);
        if (leaves) printf("  [through %d-1-%d: %lld]\n", nb[k], nb[k2], leaves);
    }
}

/* ---- exists mode: terminal-form cycles (8-1-3 for even N, 15-1-3 for odd N), early stop ---- */
static void noop_leaf(int len) { (void)len; }
static void mode_exists(int lo, int hi) {
    for (int n = lo; n <= hi; n++) {
        int endv = (n % 2 == 0) ? 8 : 15;
        build(n);
        path[0] = 1; path[1] = 3;
        on_leaf = noop_leaf; want_end = endv; stop_after = 1; leaves = 0;
        dfs((1ULL << 1) | (1ULL << 3), 2);
        long long a = leaves;
        int other = (n % 2 == 0) ? 15 : 8;
        build(n); path[0] = 1; path[1] = 3; want_end = other; stop_after = 1; leaves = 0;
        dfs((1ULL << 1) | (1ULL << 3), 2);
        printf("N=%d parity_terminal(%d-1-3)=%s other(%d-1-3)=%s\n", n, endv, a ? "yes" : "NO", other, leaves ? "yes" : "no");
        fflush(stdout);
    }
}
/* ---- pairshow mode: print the compatible ninja pairs at n (Shorter matched against Longer masks) ---- */
static int pp_n;
static u64 target_mask;
static void show_matching(int len) {
    u64 m = parity_mask(len);
    if (m == target_mask) { printf("    "); for (int j = 0; j < len; j++) printf("%d%s", path[j], j + 1 < len ? "," : "\n"); }
}
static void mode_pairshow(int n, int limit) {
    int endS = (n % 2 == 0) ? 8 : 10, endL = ((n + 1) % 2 == 0) ? 8 : 10;
    nmask = 0; mask_upto = n;
    enumerate_from13(n + 1, endL, store_mask);
    qsort(masks, nmask, sizeof(u64), cmp64);
    long long uniq = 0;
    for (long long i = 0; i < nmask; i++) if (i == 0 || masks[i] != masks[i - 1]) masks[uniq++] = masks[i];
    nmask = uniq;
    /* collect matching shorter masks */
    static u64 hitm[4096]; static int nh; nh = 0;
    pair_hits = 0;
    build(n); path[0] = 1; path[1] = 3;
    /* first pass: list the matching shorter paths */
    printf("n=%d: compatible Shorter paths (order %d, end %d):\n", n, n, endS);
    on_leaf = probe_mask; want_end = endS; stop_after = -1; leaves = 0;
    /* reuse probe to count, then print by re-running with explicit masks */
    dfs((1ULL << 1) | (1ULL << 3), 2);
    printf("  matches=%lld\n", pair_hits);
    (void)hitm; (void)nh; (void)limit; (void)pp_n;
}

/* ---- pairdump mode: print every compatible (Shorter, Longer) ninja pair at n ---- */
#define MAXSTORE 2000000
static unsigned char *Lpaths; static u64 *Lmask; static long long nL_store; static int Lorder;
static void store_longer(int len) {
    if (nL_store >= MAXSTORE) return;
    Lmask[nL_store] = parity_mask(len);
    for (int j = 0; j < len; j++) Lpaths[nL_store * 64 + j] = (unsigned char)path[j];
    nL_store++;
}
static long long dumped; static int dump_limit;
static void dump_matching(int len) {
    u64 m = parity_mask(len);
    for (long long i = 0; i < nL_store; i++) if (Lmask[i] == m) {
        dumped++;
        if (dumped <= dump_limit) {
            printf("  S: "); for (int j = 0; j < len; j++) printf("%d%s", path[j], j + 1 < len ? "," : "\n");
            printf("  L: "); for (int j = 0; j < Lorder; j++) printf("%d%s", Lpaths[i * 64 + j], j + 1 < Lorder ? "," : "\n");
        }
    }
}
static void mode_pairdump(int n, int limit) {
    int endS = (n % 2 == 0) ? 8 : 10, endL = ((n + 1) % 2 == 0) ? 8 : 10;
    Lpaths = malloc((size_t)MAXSTORE * 64); Lmask = malloc(sizeof(u64) * MAXSTORE);
    nL_store = 0; mask_upto = n; Lorder = n + 1;
    enumerate_from13(n + 1, endL, store_longer);
    dumped = 0; dump_limit = limit;
    enumerate_from13(n, endS, dump_matching);
    printf("n=%d compatible (Shorter,Longer) path pairs = %lld\n", n, dumped);
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s pairs|paths LO HI\n", argv[0]); return 1; }
    int lo = atoi(argv[2]), hi = atoi(argv[3]);
    int order_hi = (!strcmp(argv[1], "pairdump") || !strcmp(argv[1], "cycles") || !strcmp(argv[1], "show")) ? lo : hi;
    if (order_hi + 1 > MAXN - 1) { fprintf(stderr, "order too large\n"); return 1; }
    if (!strcmp(argv[1], "pairs")) mode_pairs(lo, hi);
    else if (!strcmp(argv[1], "paths")) mode_paths(lo, hi);
    else if (!strcmp(argv[1], "show")) mode_show(lo, hi, atoi(argv[4]), atoi(argv[5]), argc > 6 ? atoi(argv[6]) : 20);
    else if (!strcmp(argv[1], "cycles")) mode_cycles(lo, hi);
    else if (!strcmp(argv[1], "exists")) mode_exists(lo, hi);
    else if (!strcmp(argv[1], "pairdump")) mode_pairdump(lo, hi);
    else { fprintf(stderr, "unknown mode\n"); return 1; }
    return 0;
}
