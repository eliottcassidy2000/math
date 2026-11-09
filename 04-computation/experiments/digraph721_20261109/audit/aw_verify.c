/* aw_verify.c - AUDIT D1 (THM-5055): independent verifier for explicit graphs.
 *
 * Written from scratch by the auditor; shares no code with wcheck.c / dom_even.c, and uses different algorithms:
 *   - d(G) by inclusion-exclusion over the set U of undominated vertices:
 *         d(G) = sum_{U subset V} (-1)^|U| 2^(n - |N[U]|)
 *     (sets S with S cap N[U] = empty leave all of U undominated).
 *   - von Bromssen parity: a strong generating set of Aut(G) from a stabiliser chain whose base is chosen by
 *     individualisation + equitable refinement (nauty-style ordered partitions, refinement with a splitter queue and
 *     a trace hash for pruning); every candidate image in the base point's cell is tested by a complete
 *     individualisation-refinement search, and each automorphism found is checked explicitly.  Since the sign is a
 *     homomorphism and the transversals generate Aut(G), G is even iff all found automorphisms are even.
 *   - flags from degrees.
 * Modes:
 *   aw_verify            read "tag set n k | u-v u-v ..." lines on stdin, print one result line per graph + summary
 *   aw_verify selftest N  compare against brute force over all n! permutations and direct subset enumeration on
 *                         N random graphs with 1..8 vertices (plus fixed known cases)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 32
static int n;
static uint32_t adj[MAXN];

/* ---------------- domination count by inclusion-exclusion ---------------- */
static uint32_t *NU = NULL; static size_t NUcap = 0;
static long long dcount_ie(void) {
    size_t N = (size_t)1 << n;
    if (N > NUcap) { free(NU); NU = malloc(N * sizeof(uint32_t)); NUcap = N; }
    uint32_t cl[MAXN];
    for (int v = 0; v < n; v++) cl[v] = adj[v] | (1u << v);
    long long s = 1LL << n;   /* U = empty: 2^n */
    NU[0] = 0;
    for (size_t U = 1; U < N; U++) {
        int b = __builtin_ctzll(U);
        NU[U] = NU[U & (U - 1)] | cl[b];
        long long t = 1LL << (n - __builtin_popcount(NU[U]));
        if (__builtin_popcountll(U) & 1) s -= t; else s += t;
    }
    return s;
}
static long long dcount_direct(void) {   /* only for selftest */
    uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1); long long c = 0;
    for (uint32_t S = 0; S <= full; S++) {
        uint32_t dm = 0; for (int v = 0; v < n; v++) if (S >> v & 1) dm |= adj[v] | (1u << v);
        if (dm == full) c++;
        if (S == full) break;
    }
    return c;
}

/* ---------------- ordered partitions ---------------- */
typedef struct {
    int lab[MAXN];      /* vertices in order */
    int pos[MAXN];      /* position of each vertex */
    int start[MAXN];    /* start[i] = 1 iff a cell starts at position i */
    int ncells;
} Part;

static int cell_end(const Part *P, int s) { int e = s + 1; while (e < n && !P->start[e]) e++; return e; }

static uint64_t mixh(uint64_t h, uint64_t x) { h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); h *= 0xff51afd7ed558ccdULL; h ^= h >> 29; return h; }

/* refine to the coarsest equitable refinement; returns trace hash.  Label-independent: processing order depends only
   on cell positions and counts. */
static uint64_t refine(Part *P) {
    uint64_t trace = 1469598103934665603ULL;
    int queue[4 * MAXN * MAXN], qh = 0, qt = 0, inq[MAXN];
    memset(inq, 0, sizeof inq);
    for (int i = 0; i < n; i++) if (P->start[i]) { queue[qt++] = i; inq[i] = 1; }
    while (qh < qt) {
        int ws = queue[qh++]; inq[ws] = 0;
        if (!P->start[ws]) continue;   /* cannot happen: starts are never removed */
        int we = cell_end(P, ws);
        uint32_t wm = 0; for (int i = ws; i < we; i++) wm |= 1u << P->lab[i];
        int cnt[MAXN]; for (int v = 0; v < n; v++) cnt[v] = __builtin_popcount(adj[v] & wm);
        /* scan cells in position order */
        int s = 0;
        while (s < n) {
            int e = cell_end(P, s);
            if (e - s > 1) {
                int mn = 1 << 30, mx = -1;
                for (int i = s; i < e; i++) { int c = cnt[P->lab[i]]; if (c < mn) mn = c; if (c > mx) mx = c; }
                if (mn != mx) {
                    /* stable counting sort of the cell by count ascending */
                    int tmp[MAXN], m = 0;
                    int vals[MAXN], nv = 0;
                    for (int c = mn; c <= mx; c++) {
                        int any = 0;
                        for (int i = s; i < e; i++) if (cnt[P->lab[i]] == c) { tmp[m++] = P->lab[i]; any++; }
                        if (any) vals[nv++] = c;
                    }
                    for (int i = 0; i < m; i++) { P->lab[s + i] = tmp[i]; P->pos[tmp[i]] = s + i; }
                    /* new cell starts */
                    int p = s;
                    trace = mixh(trace, (uint64_t)ws * 131 + s); trace = mixh(trace, nv);
                    for (int t = 0; t < nv; t++) {
                        int sz = 0; for (int i = s; i < e; i++) if (cnt[P->lab[i]] == vals[t]) sz++;
                        trace = mixh(trace, (uint64_t)vals[t] * 1000 + sz);
                        if (t > 0) { P->start[p] = 1; P->ncells++; }
                        if (!inq[p]) { queue[qt++] = p; inq[p] = 1; }
                        p += sz;
                    }
                }
            }
            s = e;
        }
    }
    return trace;
}

static void unit_part(Part *P) {
    for (int i = 0; i < n; i++) { P->lab[i] = i; P->pos[i] = i; P->start[i] = (i == 0); }
    P->ncells = (n > 0);
}

/* individualise vertex v (which lies in the cell starting at s): move v to position s, split */
static void individualise(Part *P, int s, int v) {
    int pv = P->pos[v];
    int u = P->lab[s];
    P->lab[s] = v; P->pos[v] = s; P->lab[pv] = u; P->pos[u] = pv;
    if (s + 1 < n && !P->start[s + 1]) { P->start[s + 1] = 1; P->ncells++; }
}

static int same_shape(const Part *A, const Part *B) {
    if (A->ncells != B->ncells) return 0;
    for (int i = 0; i < n; i++) if (A->start[i] != B->start[i]) return 0;
    return 1;
}

static int is_aut(const int *g) {
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++)
        if (((adj[u] >> v) & 1) != ((adj[g[u]] >> g[v]) & 1)) return 0;
    /* bijectivity */
    uint32_t seen = 0; for (int u = 0; u < n; u++) seen |= 1u << g[u];
    return __builtin_popcount(seen) == n;
}

static long nsearch_nodes;
/* complete search for an automorphism g with g(cell_i of L) = cell_i of R; L, R refined and same shape */
static int search(const Part *L, const Part *R, int *g) {
    nsearch_nodes++;
    if (L->ncells == n) {
        for (int i = 0; i < n; i++) g[L->lab[i]] = R->lab[i];
        return is_aut(g);
    }
    int s = 0; while (s < n && cell_end(L, s) - s == 1) s = cell_end(L, s);
    int e = cell_end(L, s);
    Part L2 = *L; int x = L->lab[s];
    individualise(&L2, s, x); uint64_t tl = refine(&L2);
    for (int i = s; i < e; i++) {
        Part R2 = *R; int y = R->lab[i];
        individualise(&R2, s, y); uint64_t tr = refine(&R2);
        if (tl != tr || !same_shape(&L2, &R2)) continue;
        if (search(&L2, &R2, g)) return 1;
    }
    return 0;
}

static int sign_odd(const int *g) {
    int par = 0;
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++)
        if (((adj[u] >> v) & 1) && g[u] > g[v]) par ^= 1;
    return par;
}

static long ngens;
/* returns 1 if ODD, 0 if EVEN */
static int parity_odd(void) {
    Part P; unit_part(&P); refine(&P);
    ngens = 0;
    while (P.ncells < n) {
        int s = 0; while (s < n && cell_end(&P, s) - s == 1) s = cell_end(&P, s);
        int e = cell_end(&P, s);
        int b = P.lab[s];
        Part PL = P; individualise(&PL, s, b); uint64_t tl = refine(&PL);
        for (int i = s; i < e; i++) {
            int c = P.lab[i]; if (c == b) continue;
            Part PR = P; individualise(&PR, s, c); uint64_t tr = refine(&PR);
            if (tl != tr || !same_shape(&PL, &PR)) continue;
            int g[MAXN];
            if (search(&PL, &PR, g)) {
                ngens++;
                if (g[b] != c) { fprintf(stderr, "internal: generator does not map base point\n"); exit(9); }
                if (sign_odd(g)) return 1;
            }
        }
        P = PL;
    }
    return 0;
}

/* brute force parity over all permutations (selftest only, n <= 8) */
static int bf_perm[MAXN], bf_used[MAXN], bf_odd;
static void bf_rec(int i) {
    if (bf_odd) return;
    if (i == n) { if (sign_odd(bf_perm)) bf_odd = 1; return; }
    for (int c = 0; c < n; c++) if (!bf_used[c]) {
        int ok = 1;
        for (int j = 0; j < i; j++) if (((adj[i] >> j) & 1) != ((adj[c] >> bf_perm[j]) & 1)) { ok = 0; break; }
        if (!ok) continue;
        bf_used[c] = 1; bf_perm[i] = c; bf_rec(i + 1); bf_used[c] = 0;
    }
}
static int parity_odd_bf(void) { memset(bf_used, 0, sizeof bf_used); bf_odd = 0; bf_rec(0); return bf_odd; }

static void flags(int *U, int *pend, int *pnu, int *r3) {
    int deg[MAXN]; for (int v = 0; v < n; v++) deg[v] = __builtin_popcount(adj[v]);
    *U = 0; *pend = 0; *pnu = 0; *r3 = 0;
    for (int v = 0; v < n; v++) if (deg[v] == n - 1) *U = 1;
    for (int v = 0; v < n; v++) {
        if (deg[v] == 1) { (*pend)++; int q = __builtin_ctz(adj[v]); if (deg[q] != n - 1) (*pnu)++; }
        if (deg[v] == 2) *r3 = 1;
    }
}

static uint64_t rng = 88172645463325252ULL;
static uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

static void set_edges(int nn, const int (*ed)[2], int m) {
    n = nn; memset(adj, 0, sizeof adj);
    for (int i = 0; i < m; i++) { adj[ed[i][0]] |= 1u << ed[i][1]; adj[ed[i][1]] |= 1u << ed[i][0]; }
}

static int selftest(int N) {
    int bad = 0;
    /* known cases: K2 odd, P3 even, P4 odd, C4 even, K33 odd, K22 even (= C4), K4-e odd, paw odd, E5 even */
    struct { const char *name; int nn; int m; int ed[16][2]; int odd; long long d; } kc[] = {
        {"K2", 2, 1, {{0,1}}, 1, 3},
        {"P3", 3, 2, {{0,1},{1,2}}, 0, 5},
        {"K3", 3, 3, {{0,1},{1,2},{0,2}}, 1, 7},
        {"P4", 4, 3, {{0,1},{1,2},{2,3}}, 1, 9},
        {"C4", 4, 4, {{0,1},{1,2},{2,3},{3,0}}, 0, 11},
        {"K4-e", 4, 5, {{0,1},{0,2},{0,3},{1,2},{1,3}}, 1, 13},
        {"K33", 6, 9, {{0,3},{0,4},{0,5},{1,3},{1,4},{1,5},{2,3},{2,4},{2,5}}, 1, -1},
        {"E5", 5, 0, {{0,0}}, 0, 1},
        {"P3uK1", 4, 2, {{0,1},{1,2}}, 0, 5},
    };
    for (size_t t = 0; t < sizeof kc / sizeof kc[0]; t++) {
        set_edges(kc[t].nn, (const int (*)[2])kc[t].ed, kc[t].m);
        int o = parity_odd(); long long d = dcount_ie();
        if (o != kc[t].odd || (kc[t].d > 0 && d != kc[t].d)) { bad++; printf("KNOWN FAIL %s odd=%d d=%lld\n", kc[t].name, o, d); }
    }
    for (int it = 0; it < N; it++) {
        n = 1 + (int)(xr() % 8); memset(adj, 0, sizeof adj);
        int p = (int)(xr() % 100);
        for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if ((int)(xr() % 100) < p) { adj[u] |= 1u << v; adj[v] |= 1u << u; }
        /* sometimes add twins / isolated structure: random relabel of a symmetric pattern */
        if (it % 3 == 0 && n >= 4) { /* make vertex 1 a twin of vertex 0 (adjacent or not) */
            uint32_t m0 = adj[0] & ~(1u << 1) & ~1u;
            for (int v = 0; v < n; v++) { adj[v] &= ~(1u << 1); }
            adj[1] = m0; for (int v = 0; v < n; v++) if (m0 >> v & 1) adj[v] |= 1u << 1;
            if (xr() & 1) { adj[0] |= 2; adj[1] |= 1; } else { adj[0] &= ~2u; adj[1] &= ~1u; }
        }
        int o1 = parity_odd(), o2 = parity_odd_bf();
        long long d1 = dcount_ie(), d2 = dcount_direct();
        if (o1 != o2 || d1 != d2) { bad++; if (bad < 20) printf("FAIL n=%d odd(sc)=%d odd(bf)=%d d(ie)=%lld d(direct)=%lld\n", n, o1, o2, d1, d2); }
    }
    printf("selftest: %d random graphs + known cases, failures %d\n", N, bad);
    return bad != 0;
}

#ifndef AW_NO_MAIN
int main(int argc, char **argv) {
    if (argc > 1 && !strcmp(argv[1], "selftest")) return selftest(argc > 2 ? atoi(argv[2]) : 20000);
    int bfmode = (argc > 1 && !strcmp(argv[1], "bf"));   /* also brute-force parity (n <= 11) and direct d */
    static char line[1 << 17];
    long ok = 0, bad = 0, lines = 0;
    while (fgets(line, sizeof line, stdin)) {
        char tag[128], set[16]; long long k; int pos = 0;
        if (sscanf(line, "%127s %15s %d %lld |%n", tag, set, &n, &k, &pos) < 4 || pos == 0) continue;
        if (n < 1 || n > 26) { printf("SKIP %s n=%d\n", tag, n); continue; }
        memset(adj, 0, sizeof adj);
        char *p = line + pos; int u, v, m, badedge = 0;
        while (sscanf(p, " %d-%d%n", &u, &v, &m) == 2) {
            if (u < 0 || v < 0 || u >= n || v >= n || u == v) badedge = 1;
            else { adj[u] |= 1u << v; adj[v] |= 1u << u; }
            p += m;
        }
        lines++;
        long long d = dcount_ie();
        int odd = parity_odd();
        int U, pend, pnu, r3; flags(&U, &pend, &pnu, &r3);
        int fl = 1;
        if (!strcmp(set, "W")) fl = !U;
        else if (!strcmp(set, "P")) fl = pend > 0;
        else if (!strcmp(set, "S")) fl = pnu > 0;
        else if (!strcmp(set, "S2")) fl = pnu >= 2;
        else if (!strcmp(set, "S3")) fl = r3;
        else if (!strcmp(set, "G")) fl = (pnu >= 2) && r3;
        else if (!strcmp(set, "E")) fl = 1;
        else fl = -1;
        if (bfmode && n <= 11) {
            int o2 = parity_odd_bf(); long long d2 = dcount_direct();
            printf("BF %s n=%d odd(sc)=%d odd(bf)=%d d(ie)=%lld d(direct)=%lld %s\n", tag, n, odd, o2, d, d2,
                   (o2 == odd && d2 == d) ? "AGREE" : "DISAGREE");
            if (o2 != odd || d2 != d) bad++;
        }
        int good = (d == k) && !odd && fl == 1 && !badedge;
        printf("%s %s %d %lld d=%lld %s U=%d pend=%d pnu=%d r3=%d flags=%d gens=%ld %s\n", tag, set, n, k, d,
               odd ? "ODD" : "EVEN", U, pend, pnu, r3, fl, ngens, good ? "OK" : "BAD");
        if (good) ok++; else bad++;
    }
    printf("SUMMARY lines=%ld ok=%ld bad=%ld\n", lines, ok, bad);
    return bad != 0;
}
#endif
