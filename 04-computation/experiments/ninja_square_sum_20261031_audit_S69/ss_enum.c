/* S69 audit (independent referee, 2026-10-31): exhaustive square-sum enumerator written from scratch.
   It does not reuse ninja_pairs_square_sum_20261031.c.

   Q_N: vertices 1..N (N <= 62), u ~ v iff u + v is a perfect square. Sets are uint64 bitmasks (bit v = vertex v).

   Search = depth-first extension of a path, state = (unvisited set U, current end cur), with three sound rules:
     R1  every v in U has >= 1 neighbour in S = U + {cur}; a vertex with exactly 1 must be the last vertex of the
         path (so it must be the target when a target is fixed; at most one such vertex when the end is free);
     R2  (target fixed) a non-target v in U adjacent to cur with exactly 2 neighbours in S is interior with path
         neighbours = those two, so it must be cur's successor: forced move (two such vertices: dead);
     R3  S induces a connected subgraph (bitmask flood fill from cur).

   Modes (one result line per order unless stated):
     census LO HI        undirected Hamiltonian paths (sum over starts of directed paths, / 2); undirected
                         Hamiltonian cycles (rooted at the LARGEST vertex N: paths N,a,...,b over neighbour pairs a<b);
                         paths 1,3->8, 1,3->15, 1,3->10, 1->3 (start 1, end 3), 1->8 (start 1, end 8)
     pairs LO HI KIND    pair census for n in [LO,HI]. KIND B = weak B (both start 1,3; even order ends 8, odd ends 10);
                         KIND G = Gerbicz nice pair normalisation (both start 1; odd order ends 3, even order ends 8);
                         KIND F = free (all directed Hamiltonian paths, no prefix/end condition).
                         Compatibility = equal one-based positional parity of every vertex <= n.
     pairdump n MAX [KIND]   print up to MAX compatible pairs as lines "  S: ..." and "  L: ..." (session-compatible
                         format, so the session's replay script can be fed by THIS program), then a count line
     list N A B END [CAP]  print every (or the first CAP) Hamiltonian path with prefix A (and B if B > 0) ending at
                         END (END = 0: any end)
     cycles N            print every undirected Hamiltonian cycle as a vertex list starting N,a,...,b with a < b
     terminal LO HI CAP  per N: number of paths 1,3 -> 8 (N even) or 1,3 -> 15 (N odd), stopping after CAP (CAP<0: all);
                         also the opposite closing; complete-search flag
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned long long u64;

static int N;
static u64 adjm[64];
static int path[64];
static int target;              /* 0 = free end */
static long long leaves, nodes, stop_at;
static void (*leaf_cb)(void);

static int issq(int x) { int r = 0; while ((r + 1) * (r + 1) <= x) r++; return r * r == x; }

static void setup(int n) {
    N = n;
    for (int v = 0; v < 64; v++) adjm[v] = 0;
    for (int a = 1; a <= n; a++)
        for (int b = a + 1; b <= n; b++)
            if (issq(a + b)) { adjm[a] |= 1ULL << b; adjm[b] |= 1ULL << a; }
}

static int flood_ok(u64 S, int src) {
    u64 reach = 1ULL << src, front = reach;
    while (front) {
        u64 nxt = 0, f = front;
        while (f) { int v = __builtin_ctzll(f); f &= f - 1; nxt |= adjm[v]; }
        nxt &= S & ~reach;
        reach |= nxt;
        front = nxt;
    }
    return reach == S;
}

static void ext(u64 U, int len) {
    nodes++;
    if (stop_at >= 0 && leaves >= stop_at) return;
    int cur = path[len - 1];
    if (!U) {
        if (!target || cur == target) { leaves++; if (leaf_cb) leaf_cb(); }
        return;
    }
    if (target && cur == target) return;
    u64 S = U | (1ULL << cur);
    int forced = 0, ones = 0;
    for (u64 t = U; t; t &= t - 1) {
        int v = __builtin_ctzll(t);
        int c = __builtin_popcountll(adjm[v] & S);
        if (c == 0) return;
        if (c == 1) {
            if (target) { if (v != target) return; }
            else if (++ones > 1) return;
        }
        if (target && c == 2 && v != target && ((adjm[cur] >> v) & 1ULL)) {
            if (forced) return;
            forced = v;
        }
    }
    if (!flood_ok(S, cur)) return;
    if (forced) { path[len] = forced; ext(U & ~(1ULL << forced), len + 1); return; }
    for (u64 c = adjm[cur] & U; c; c &= c - 1) {
        int w = __builtin_ctzll(c);
        if (target && w == target && U != (1ULL << w)) continue;
        path[len] = w;
        ext(U & ~(1ULL << w), len + 1);
    }
}

/* run a search on Q_n with a given prefix and end; returns leaves */
static long long run(int n, const int *pre, int npre, int end, long long cap, void (*cb)(void)) {
    if (n < 1 || n > 62) { fprintf(stderr, "order out of range\n"); exit(2); }
    setup(n);
    if (npre > n) return 0;
    u64 U = (1ULL << (n + 1)) - 2;
    for (int i = 0; i < npre; i++) {
        if (pre[i] < 1 || pre[i] > n || !((U >> pre[i]) & 1ULL)) return 0;
        if (i && !((adjm[pre[i - 1]] >> pre[i]) & 1ULL)) return 0;
        path[i] = pre[i];
        U &= ~(1ULL << pre[i]);
    }
    if (end && (end > n)) return 0;
    target = end; leaves = 0; stop_at = cap; leaf_cb = cb;
    if (end && npre < n && !((U >> end) & 1ULL)) return 0;   /* end already used inside the prefix */
    ext(U, npre);
    return leaves;
}

/* ---------------- census ---------------- */
static void mode_census(int lo, int hi) {
    for (int n = lo; n <= hi; n++) {
        long long dir = 0;
        for (int s = 1; s <= n; s++) { int pre[1] = {s}; dir += run(n, pre, 1, 0, -1, NULL); }
        long long cyc = 0;
        setup(n);
        u64 nb = adjm[n];
        int list[16], k = 0;
        for (u64 t = nb; t; t &= t - 1) list[k++] = __builtin_ctzll(t);
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++) { int pre[2] = {n, list[i]}; cyc += run(n, pre, 2, list[j], -1, NULL); }
        int p13[2] = {1, 3}, p1[1] = {1};
        long long t8 = run(n, p13, 2, 8, -1, NULL), t15 = run(n, p13, 2, 15, -1, NULL), t10 = run(n, p13, 2, 10, -1, NULL);
        long long g3 = run(n, p1, 1, 3, -1, NULL), g8 = run(n, p1, 1, 8, -1, NULL);
        printf("N=%d paths=%lld cycles=%lld p13_8=%lld p13_15=%lld p13_10=%lld p1_3=%lld p1_8=%lld\n",
               n, dir / 2, cyc, t8, t15, t10, g3, g8);
        fflush(stdout);
    }
}

/* ---------------- pair machinery ---------------- */
static int mask_n;
static u64 *mbuf; static long long mcnt, mcap;
static unsigned char *pbuf; static long long pcap; static int keep_paths;
static void push_mask(void) {
    u64 m = 0;
    for (int j = 0; j < N; j += 2) if (path[j] <= mask_n) m |= 1ULL << path[j];   /* one-based odd positions */
    if (mcnt == mcap) {
        mcap = mcap ? 2 * mcap : 4096;
        mbuf = realloc(mbuf, mcap * sizeof(u64));
        if (!mbuf) { fprintf(stderr, "out of memory\n"); exit(2); }
    }
    if (keep_paths) {
        if (mcnt >= pcap) {
            pcap = pcap ? 2 * pcap : 4096;
            pbuf = realloc(pbuf, pcap * 64);
            if (!pbuf) { fprintf(stderr, "out of memory\n"); exit(2); }
        }
        for (int j = 0; j < N; j++) pbuf[mcnt * 64 + j] = (unsigned char)path[j];
    }
    mbuf[mcnt++] = m;
}
static int cmpu(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return (x > y) - (x < y); }

static int end_for(char kind, int order) {
    if (kind == 'B') return (order % 2 == 0) ? 8 : 10;
    if (kind == 'G') return (order % 2 == 0) ? 8 : 3;
    return 0;
}
/* collect all paths of the given order for the pair kind; returns number of paths; masks in mbuf[0..mcnt) */
static long long collect(char kind, int order, int maskn, int keep) {
    mcnt = 0; mask_n = maskn; keep_paths = keep;
    long long tot = 0;
    if (kind == 'B') { int pre[2] = {1, 3}; tot = run(order, pre, 2, end_for('B', order), -1, push_mask); }
    else if (kind == 'G') { int pre[1] = {1}; tot = run(order, pre, 1, end_for('G', order), -1, push_mask); }
    else { for (int s = 1; s <= order; s++) { int pre[1] = {s}; tot += run(order, pre, 1, 0, -1, push_mask); } }
    return tot;
}

static long long count_eq(const u64 *arr, long long n, u64 key) {   /* arr sorted */
    long long lo = 0, hi = n;
    while (lo < hi) { long long mid = (lo + hi) / 2; if (arr[mid] < key) lo = mid + 1; else hi = mid; }
    long long a = lo; hi = n;
    while (lo < hi) { long long mid = (lo + hi) / 2; if (arr[mid] <= key) lo = mid + 1; else hi = mid; }
    return lo - a;
}

static void mode_pairs(int lo, int hi, char kind) {
    for (int n = lo; n <= hi; n++) {
        long long nL = collect(kind, n + 1, n, 0);
        u64 *Lm = malloc((mcnt ? mcnt : 1) * sizeof(u64)); memcpy(Lm, mbuf, mcnt * sizeof(u64)); long long Lc = mcnt;
        qsort(Lm, Lc, sizeof(u64), cmpu);
        long long nS = collect(kind, n, n, 0);
        u64 *Sm = malloc((mcnt ? mcnt : 1) * sizeof(u64)); memcpy(Sm, mbuf, mcnt * sizeof(u64)); long long Sc = mcnt;
        qsort(Sm, Sc, sizeof(u64), cmpu);
        long long Sok = 0, Lok = 0, pairs = 0, Ldist = 0;
        for (long long i = 0; i < Sc; i++) { long long c = count_eq(Lm, Lc, Sm[i]); if (c) { Sok++; pairs += c; } }
        for (long long i = 0; i < Lc; i++) { if (count_eq(Sm, Sc, Lm[i])) Lok++; if (!i || Lm[i] != Lm[i - 1]) Ldist++; }
        printf("kind=%c n=%d S_paths=%lld L_paths=%lld L_distinct_masks=%lld S_compatible=%lld L_compatible=%lld pairs=%lld exists=%s\n",
               kind, n, nS, nL, Ldist, Sok, Lok, pairs, pairs ? "YES" : "no");
        fflush(stdout);
        free(Lm); free(Sm);
    }
}

typedef struct { u64 m; long long idx; } MI;
static int cmpmi(const void *a, const void *b) { const MI *x = a, *y = b; return (x->m > y->m) - (x->m < y->m); }

static void mode_pairdump(int n, long long maxp, char kind) {
    long long nL = collect(kind, n + 1, n, 1);
    long long Lc = mcnt;
    MI *L = malloc((Lc ? Lc : 1) * sizeof(MI));
    for (long long i = 0; i < Lc; i++) { L[i].m = mbuf[i]; L[i].idx = i; }
    unsigned char *Lp = malloc((Lc ? Lc : 1) * 64); memcpy(Lp, pbuf, Lc * 64);
    qsort(L, Lc, sizeof(MI), cmpmi);
    long long nS = collect(kind, n, n, 1);
    long long Sc = mcnt, printed = 0, total = 0;
    for (long long i = 0; i < Sc; i++) {
        u64 key = mbuf[i];
        long long lo = 0, hi = Lc;
        while (lo < hi) { long long mid = (lo + hi) / 2; if (L[mid].m < key) lo = mid + 1; else hi = mid; }
        for (long long k = lo; k < Lc && L[k].m == key; k++) {
            total++;
            if (printed < maxp) {
                printed++;
                printf("  S: "); for (int j = 0; j < n; j++) printf("%d%s", pbuf[i * 64 + j], j + 1 < n ? "," : "\n");
                printf("  L: "); for (int j = 0; j < n + 1; j++) printf("%d%s", Lp[L[k].idx * 64 + j], j + 1 < n + 1 ? "," : "\n");
            }
        }
    }
    printf("kind=%c n=%d S_paths=%lld L_paths=%lld compatible (Shorter,Longer) path pairs = %lld\n", kind, n, nS, nL, total);
}

/* ---------------- list / cycles / terminal ---------------- */
static void print_path(void) { for (int j = 0; j < N; j++) printf("%d%s", path[j], j + 1 < N ? "," : "\n"); }

static void mode_list(int n, int a, int b, int end, long long cap) {
    int pre[2] = {a, b};
    long long c = run(n, pre, b > 0 ? 2 : 1, end, cap, print_path);
    printf("total=%lld%s\n", c, (cap >= 0 && c >= cap) ? "+" : "");
}

static void mode_cycles(int n) {
    setup(n);
    int list[16], k = 0;
    for (u64 t = adjm[n]; t; t &= t - 1) list[k++] = __builtin_ctzll(t);
    long long tot = 0;
    for (int i = 0; i < k; i++)
        for (int j = i + 1; j < k; j++) { int pre[2] = {n, list[i]}; tot += run(n, pre, 2, list[j], -1, print_path); }
    printf("total=%lld\n", tot);
}

static void mode_terminal(int lo, int hi, long long cap) {
    for (int n = lo; n <= hi; n++) {
        int p13[2] = {1, 3};
        int e = (n % 2 == 0) ? 8 : 15, o = (n % 2 == 0) ? 15 : 8;
        nodes = 0;
        long long a = run(n, p13, 2, e, cap, NULL);
        long long na = nodes;
        int complete_a = (cap < 0) || (a < cap);
        nodes = 0;
        long long b = run(n, p13, 2, o, cap, NULL);
        int complete_b = (cap < 0) || (b < cap);
        printf("N=%d terminal(%d-1-3)=%lld%s other(%d-1-3)=%lld%s nodes=%lld\n", n, e, a, complete_a ? "" : "+", o, b,
               complete_b ? "" : "+", na);
        fflush(stdout);
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: see header\n"); return 1; }
    const char *m = argv[1];
    if (!strcmp(m, "census")) mode_census(atoi(argv[2]), atoi(argv[3]));
    else if (!strcmp(m, "pairs")) mode_pairs(atoi(argv[2]), atoi(argv[3]), argc > 4 ? argv[4][0] : 'B');
    else if (!strcmp(m, "pairdump")) mode_pairdump(atoi(argv[2]), atoll(argv[3]), argc > 4 ? argv[4][0] : 'B');
    else if (!strcmp(m, "list")) mode_list(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), argc > 6 ? atoll(argv[6]) : -1);
    else if (!strcmp(m, "cycles")) mode_cycles(atoi(argv[2]));
    else if (!strcmp(m, "terminal")) mode_terminal(atoi(argv[2]), atoi(argv[3]), atoll(argv[4]));
    else { fprintf(stderr, "unknown mode\n"); return 1; }
    return 0;
}
