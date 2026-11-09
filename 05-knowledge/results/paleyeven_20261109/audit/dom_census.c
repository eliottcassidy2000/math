/* dom_census.c - audit P1 (independent implementation): connected graphs with at most K dominating sets,
 * their Royle parity, |Aut| for even ones, and I(G,2).
 *
 * Differences from the session's dom_even.c (written without reusing it):
 *  - isomorphism: exact CANONICAL FORMS by individualisation-refinement (equitable refinement, first non-singleton
 *    target cell, full search tree with no pruning, canonical form = lexicographically least relabelled adjacency
 *    string over all leaves). Two graphs are isomorphic iff their forms agree.
 *  - parity: the set of leaves attaining the least form is a single Aut(G)-orbit on which Aut(G) acts regularly, so
 *    {l o l0^-1} enumerates ALL of Aut(G); the sign of every automorphism is computed (early exit at the first odd
 *    one). The session instead used stabiliser-chain transversals.
 *  - d of a child G + v_N from d(G + v_N) = A(N) + d(G) - C(N) with A(N) = #{S : V \ dom(S) subset of N},
 *    C(N) = #{dominating S : S cap N empty} (my own zeta transforms); this formula is tested against direct counts
 *    in mode "test", and every stored class has its d recounted directly.
 * Exhaustiveness: for connected G and a new vertex v with a neighbour in G, d(G + v) > d(G) (S -> S + v injects,
 * V(G) is extra). Every connected graph has a connected-prefix vertex order (BFS), so level-by-level extension of
 * all classes with d <= K by every nonempty neighbourhood reaches every connected graph with d <= K.
 *
 * usage: dom_census census K maxn            > classes.txt   (stderr: per-level summary)
 *        dom_census trees maxn               > trees.txt     (all trees up to maxn vertices, by leaf addition)
 *        dom_census test                                     (formula and canonical-form self-tests)
 * class line: n m d parity naut I2 edges     (parity E/O; naut = |Aut| for E, 0 for O; edges as hex pairs uv)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned __int128 u128;
typedef unsigned long long u64;
#define MAXN 16

/* ------------------------------------------------------------------ canonical forms */
static int N;
static uint16_t ADJ[MAXN];

typedef struct { int8_t lab[MAXN]; uint32_t starts; } part;   /* bit s of starts: a cell starts at position s */

static inline int cell_end(uint32_t starts, int s) {
    uint32_t rest = starts & ~((2u << s) - 1u);
    return rest ? __builtin_ctz(rest) : N;
}

static void refine(part *p, uint32_t queue) {
    while (queue) {
        int s = __builtin_ctz(queue); queue &= queue - 1;
        int e = cell_end(p->starts, s);
        uint32_t W = 0; for (int i = s; i < e; i++) W |= 1u << p->lab[i];
        int c = 0;
        while (c < N) {
            int ce = cell_end(p->starts, c);
            if (ce - c > 1) {
                int cnt[MAXN], mn = 99, mx = -1;
                for (int i = c; i < ce; i++) {
                    int k = __builtin_popcount((uint32_t)ADJ[(int)p->lab[i]] & W);
                    cnt[i] = k; if (k < mn) mn = k; if (k > mx) mx = k;
                }
                if (mn != mx) {
                    for (int i = c + 1; i < ce; i++) {
                        int kv = cnt[i]; int8_t lv = p->lab[i]; int j = i - 1;
                        while (j >= c && cnt[j] > kv) { cnt[j + 1] = cnt[j]; p->lab[j + 1] = p->lab[j]; j--; }
                        cnt[j + 1] = kv; p->lab[j + 1] = lv;
                    }
                    for (int i = c + 1; i < ce; i++) if (cnt[i] != cnt[i - 1]) p->starts |= 1u << i;
                    for (int i = c; i < ce; i++) if (p->starts >> i & 1) queue |= 1u << i;
                }
            }
            c = ce;
        }
    }
}

static int mode;                 /* 0 = find least form, 1 = enumerate automorphisms */
static u128 best; static int have_best;
static int8_t reflab[MAXN]; static int have_ref, found_odd; static long naut; static int early_exit;
static long leaves;

static u128 leaf_form(const int8_t *lab) {
    u128 f = 0;
    for (int i = 0; i < N; i++) { int u = lab[i]; for (int j = i + 1; j < N; j++) f = (f << 1) | (u128)((ADJ[u] >> lab[j]) & 1); }
    return f;
}

static void search(part *p, uint32_t queue) {
    refine(p, queue);
    if (__builtin_popcount(p->starts) == N) {
        u128 f = leaf_form(p->lab); leaves++;
        if (mode == 0) { if (!have_best || f < best) { best = f; have_best = 1; } }
        else if (f == best) {
            if (!have_ref) { memcpy(reflab, p->lab, N); have_ref = 1; }
            int g[MAXN];
            for (int i = 0; i < N; i++) g[(int)reflab[i]] = p->lab[i];
            naut++;
            int par = 0;
            for (int u = 0; u < N; u++) {
                uint32_t r = (uint32_t)ADJ[u] & ~((2u << u) - 1u);
                while (r) { int w = __builtin_ctz(r); r &= r - 1; if (g[u] > g[w]) par ^= 1; }
            }
            if (par) found_odd = 1;
        }
        return;
    }
    int s = 0;
    for (;;) { int e = cell_end(p->starts, s); if (e - s > 1) break; s = e; }
    int e = cell_end(p->starts, s);
    for (int i = s; i < e; i++) {
        part q = *p;
        int8_t v = q.lab[i]; q.lab[i] = q.lab[s]; q.lab[s] = v;
        q.starts |= 1u << (s + 1);
        search(&q, (1u << s) | (1u << (s + 1)));
        if (mode == 1 && found_odd && early_exit) return;
    }
}

static u128 canon(int n, const uint16_t *adj) {
    N = n; memcpy(ADJ, adj, n * sizeof(uint16_t));
    part p; for (int i = 0; i < n; i++) p.lab[i] = (int8_t)i; p.starts = 1;
    mode = 0; have_best = 0; leaves = 0;
    search(&p, 1);
    return best;
}

/* returns 1 if some automorphism reverses an odd number of edges (reference orientation u -> v for u < v) */
static int royle_odd(int n, const uint16_t *adj, long *nautp) {
    canon(n, adj);
    mode = 1; have_ref = 0; found_odd = 0; naut = 0; early_exit = 1;
    part p; for (int i = 0; i < n; i++) p.lab[i] = (int8_t)i; p.starts = 1;
    search(&p, 1);
    mode = 0;
    *nautp = found_odd ? 0 : naut;
    return found_odd;
}

/* ------------------------------------------------------------------ statistics */
static uint16_t DM[1 << MAXN];
static u64 count_dom(int n, const uint16_t *adj) {
    uint32_t full = (1u << n) - 1; u64 c = 0; DM[0] = 0;
    for (uint32_t S = 1; S <= full; S++) {
        int v = __builtin_ctz(S);
        DM[S] = DM[S & (S - 1)] | adj[v] | (uint16_t)(1u << v);
        if (DM[S] == full) c++;
    }
    return c + (n == 0);
}

static u64 indep2(int n, const uint16_t *adj) {
    uint32_t full = (1u << n) - 1; u64 t = 0;
    for (uint32_t S = 0; S <= full; S++) {
        int ok = 1; uint32_t s = S;
        while (s) { int v = __builtin_ctz(s); s &= s - 1; if (adj[v] & S) { ok = 0; break; } }
        if (ok) t += 1ULL << __builtin_popcount(S);
    }
    return t;
}

/* child counts: dch[N] = d(G + v_N) for all N, via A(N) + d(G) - C(N) */
static u64 *ZA, *ZC;
static void child_counts(int n, const uint16_t *adj, u64 dG) {
    uint32_t full = (1u << n) - 1, NS = 1u << n;
    memset(ZA, 0, NS * sizeof(u64)); memset(ZC, 0, NS * sizeof(u64));
    DM[0] = 0; ZA[full] += 1;          /* S = empty: undominated set = V */
    for (uint32_t S = 1; S <= full; S++) {
        int v = __builtin_ctz(S);
        DM[S] = DM[S & (S - 1)] | adj[v] | (uint16_t)(1u << v);
        ZA[full & ~(uint32_t)DM[S]] += 1;
        if (DM[S] == full) ZC[S] += 1;
    }
    for (int b = 0; b < n; b++) for (uint32_t S = 0; S < NS; S++) if (S >> b & 1) { ZA[S] += ZA[S ^ (1u << b)]; ZC[S] += ZC[S ^ (1u << b)]; }
    (void)dG;
}
static inline u64 child_d(uint32_t Nn, int n, u64 dG) { uint32_t full = (1u << n) - 1; return ZA[Nn] + dG - ZC[full & ~Nn]; }

/* ------------------------------------------------------------------ hash set of forms */
typedef struct { u128 *key; int64_t *val; size_t cap, cnt; } hset;
static void hs_init(hset *h, size_t cap) { h->cap = cap; h->cnt = 0; h->key = calloc(cap, sizeof(u128)); h->val = malloc(cap * sizeof(int64_t)); for (size_t i = 0; i < cap; i++) h->val[i] = -1; }
static void hs_free(hset *h) { free(h->key); free(h->val); }
static inline size_t hmix(u128 k) { u64 a = (u64)k, b = (u64)(k >> 64); a ^= b * 0x9E3779B97F4A7C15ULL; a ^= a >> 31; a *= 0xBF58476D1CE4E5B9ULL; a ^= a >> 29; return (size_t)a; }
static int64_t hs_find_or_add(hset *h, u128 k, int64_t v) {   /* returns existing value or -1 after adding */
    if (2 * (h->cnt + 1) > h->cap) {
        hset g; hs_init(&g, h->cap * 2);
        for (size_t i = 0; i < h->cap; i++) if (h->val[i] >= 0) { size_t j = hmix(h->key[i]) & (g.cap - 1); while (g.val[j] >= 0) j = (j + 1) & (g.cap - 1); g.key[j] = h->key[i]; g.val[j] = h->val[i]; g.cnt++; }
        hs_free(h); *h = g;
    }
    size_t j = hmix(k) & (h->cap - 1);
    while (h->val[j] >= 0) { if (h->key[j] == k) return h->val[j]; j = (j + 1) & (h->cap - 1); }
    h->key[j] = k; h->val[j] = v; h->cnt++;
    return -1;
}

/* ------------------------------------------------------------------ levels */
typedef struct { uint16_t adj[MAXN]; u64 d; } rec;
typedef struct { rec *r; size_t cnt, cap; } lvl;
static void lv_push(lvl *L, const uint16_t *adj, int n, u64 d) {
    if (L->cnt == L->cap) { L->cap = L->cap ? 2 * L->cap : 1024; L->r = realloc(L->r, L->cap * sizeof(rec)); }
    memset(L->r[L->cnt].adj, 0, sizeof L->r[L->cnt].adj); memcpy(L->r[L->cnt].adj, adj, n * sizeof(uint16_t)); L->r[L->cnt].d = d; L->cnt++;
}

static void print_class(int n, const uint16_t *adj, u64 d, int odd, long na) {
    int m = 0; for (int v = 0; v < n; v++) m += __builtin_popcount(adj[v]); m /= 2;
    printf("%d %d %llu %c %ld %llu ", n, m, d, odd ? 'O' : 'E', na, indep2(n, adj));
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if (adj[u] >> v & 1) printf("%x%x", u, v);
    printf("\n");
}

static int run_census(u64 K, int maxn, int trees) {
    ZA = malloc(sizeof(u64) << MAXN); ZC = malloc(sizeof(u64) << MAXN);
    lvl cur = {0}; uint16_t a1[MAXN] = {0};
    lv_push(&cur, a1, 1, 1);
    for (int n = 1; n <= maxn; n++) {
        long ne = 0, no = 0;
        for (size_t i = 0; i < cur.cnt; i++) {
            rec *g = &cur.r[i];
            u64 dd = count_dom(n, g->adj);
            if (!trees && dd != g->d) { fprintf(stderr, "d mismatch at n=%d: formula %llu direct %llu\n", n, g->d, dd); return 2; }
            if (!(dd & 1)) { fprintf(stderr, "even d?\n"); return 3; }
            long na; int odd = royle_odd(n, g->adj, &na);
            print_class(n, g->adj, dd, odd, na);
            if (odd) no++; else ne++;
        }
        fprintf(stderr, "n=%d classes=%zu even=%ld odd=%ld\n", n, cur.cnt, ne, no);
        fflush(stdout); fflush(stderr);
        if (n == maxn) break;
        lvl nxt = {0}; hset hs; hs_init(&hs, 1 << 16);
        long cand = 0;
        for (size_t i = 0; i < cur.cnt; i++) {
            rec *g = &cur.r[i];
            uint32_t full = (1u << n) - 1;
            if (!trees) child_counts(n, g->adj, g->d);
            for (uint32_t Nn = 1; Nn <= full; Nn++) {
                if (trees && __builtin_popcount(Nn) != 1) continue;
                u64 dn = trees ? 0 : child_d(Nn, n, g->d);
                if (!trees && dn > K) continue;
                uint16_t h[MAXN]; memcpy(h, g->adj, n * sizeof(uint16_t)); h[n] = (uint16_t)Nn;
                for (int v = 0; v < n; v++) if (Nn >> v & 1) h[v] |= (uint16_t)(1u << n);
                cand++;
                u128 f = canon(n + 1, h);
                if (hs_find_or_add(&hs, f, (int64_t)nxt.cnt) < 0) lv_push(&nxt, h, n + 1, dn);
            }
        }
        fprintf(stderr, "   level %d: %ld candidate children canonicalised\n", n + 1, cand);
        hs_free(&hs); free(cur.r); cur = nxt;
        if (cur.cnt == 0) { fprintf(stderr, "no classes at n=%d\n", n + 1); break; }
    }
    return 0;
}

/* ------------------------------------------------------------------ self-tests */
static u64 rng = 0x1234567ULL;
static u64 xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static int run_tests(void) {
    ZA = malloc(sizeof(u64) << MAXN); ZC = malloc(sizeof(u64) << MAXN);
    long bad = 0, tests = 0;
    /* (1) child-count formula vs direct count, random graphs n <= 10, all N */
    for (int t = 0; t < 3000; t++) {
        int n = 1 + xr() % 10; uint16_t a[MAXN] = {0};
        int pm = xr() % 100;
        for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if ((int)(xr() % 100) < pm) { a[u] |= 1u << v; a[v] |= 1u << u; }
        u64 dG = count_dom(n, a);
        child_counts(n, a, dG);
        for (uint32_t Nn = 1; Nn < (1u << n); Nn++) {
            uint16_t h[MAXN]; memcpy(h, a, sizeof h); h[n] = (uint16_t)Nn;
            for (int v = 0; v < n; v++) if (Nn >> v & 1) h[v] |= (uint16_t)(1u << n);
            u64 direct = count_dom(n + 1, h), viaf = child_d(Nn, n, dG);
            tests++; if (direct != viaf) bad++;
            if (direct <= dG) { fprintf(stderr, "monotonicity fails?\n"); bad++; }
        }
    }
    printf("child formula: %ld tests, %ld mismatches (also checked d(G+v) > d(G))\n", tests, bad);
    fflush(stdout);
    /* (2) canonical form invariance under random relabelling, n <= 9 (no automorphism pruning: the search tree of
     *     K_n has n! leaves, so the test stays at n <= 9 like the census, whose densest graph is K_9) */
    long cb = 0;
    for (int t = 0; t < 20000; t++) {
        int n = 1 + xr() % 9; uint16_t a[MAXN] = {0}, b[MAXN] = {0};
        int pm = xr() % 100;
        for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if ((int)(xr() % 100) < pm) { a[u] |= 1u << v; a[v] |= 1u << u; }
        int pi[MAXN]; for (int i = 0; i < n; i++) pi[i] = i;
        for (int i = n - 1; i > 0; i--) { int j = xr() % (i + 1); int tt = pi[i]; pi[i] = pi[j]; pi[j] = tt; }
        for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (a[u] >> v & 1) b[pi[u]] |= 1u << pi[v];
        u128 fa = canon(n, a), fb = canon(n, b);
        long na, nb; int oa = royle_odd(n, a, &na), ob = royle_odd(n, b, &nb);
        if (fa != fb || oa != ob || na != nb) cb++;
    }
    printf("canonical form / parity / |Aut| relabelling invariance: 20000 random graphs, %ld failures\n", cb);
    return (bad || cb) ? 1 : 0;
}

int main(int argc, char **argv) {
    if (argc >= 2 && !strcmp(argv[1], "test")) return run_tests();
    if (argc >= 4 && !strcmp(argv[1], "census")) return run_census(strtoull(argv[2], 0, 10), atoi(argv[3]) > MAXN ? MAXN : atoi(argv[3]), 0);
    if (argc >= 3 && !strcmp(argv[1], "trees")) return run_census(0, atoi(argv[2]) > MAXN ? MAXN : atoi(argv[2]), 1);
    fprintf(stderr, "usage: dom_census census K maxn | trees maxn | test\n");
    return 1;
}
