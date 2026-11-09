/* dom_even.c - all connected graphs with at most K dominating sets, with their Royle parity (paleyeven lane,
 * mac-mini-2026-11-09).
 *
 * d(G) = number of dominating sets (S such that every vertex is in S or adjacent to S); always odd (Brouwer 2009).
 * Monotonicity: if v has at least one neighbour in G then d(G + v) >= d(G)  (S -> S u {v} is injective), with
 * equality impossible (V(G) itself dominates G + v and does not contain v).  Every connected graph has a vertex order
 * whose prefixes are connected, so ALL connected graphs with d <= K are reached by adding, level by level, a vertex
 * with a nonempty neighbourhood to the connected graphs with d <= K of the previous level.
 * Child counts: d(G + v_N) = A(N) + d(G) - C(N), where A(N) = #{S : V \ dom(S) subset of N} and
 * C(N) = #{dominating S of G disjoint from N} (two subset-zeta transforms per parent).
 * Isomorphism: 1-WL colour refinement with canonical colour names gives an invariant hash; exact backtracking
 * isomorphism tests inside a hash bucket.
 * Royle parity: reference orientation u -> v for u < v; an automorphism g reverses {u,v} iff g(u) > g(v); G is ODD iff
 * some automorphism reverses an odd number of edges.  We compute a full transversal of every level of the stabiliser
 * chain (base 0,1,...,n-1): for each base point b and each candidate image c we search for an automorphism fixing the
 * earlier base points with b -> c.  These transversals generate Aut(G), and the sign is a homomorphism, so G is even
 * iff every transversal element is even.
 * Output: one line per isomorphism class: "n m d parity I2 edges" (parity E/O; I2 = independence polynomial at 2).
 * usage: dom_even K maxn > out
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 16
typedef uint16_t row;
typedef struct { int n; row adj[MAXN]; uint64_t d; uint64_t h; } graph;

static uint64_t K;

/* ---------- WL invariant ---------- */
static uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }

static void wl_colors(const graph *g, int *col, uint64_t *hash_out) {
    int n = g->n;
    uint64_t sig[MAXN];
    for (int v = 0; v < n; v++) col[v] = __builtin_popcount(g->adj[v]);
    uint64_t H = mix(n * 1000003ULL);
    for (int round = 0; round < n; round++) {
        for (int v = 0; v < n; v++) {
            /* signature: old colour and sorted neighbour colours, hashed order-independently */
            uint64_t s = mix((uint64_t)col[v] * 0x9E3779B97F4A7C15ULL + 7);
            uint64_t acc = 0;
            row r = g->adj[v];
            while (r) { int u = __builtin_ctz(r); r &= r - 1; acc += mix((uint64_t)col[u] + 0x1234567ULL); }
            sig[v] = mix(s ^ (acc * 0x2545F4914F6CDD1DULL));
        }
        /* canonical new colours: rank of signature among distinct sorted signatures */
        uint64_t ds[MAXN]; int nd = 0;
        for (int v = 0; v < n; v++) { int f = 0; for (int i = 0; i < nd; i++) if (ds[i] == sig[v]) { f = 1; break; } if (!f) ds[nd++] = sig[v]; }
        for (int i = 1; i < nd; i++) { uint64_t t = ds[i]; int j = i - 1; while (j >= 0 && ds[j] > t) { ds[j + 1] = ds[j]; j--; } ds[j + 1] = t; }
        int newc[MAXN], changed = 0;
        for (int v = 0; v < n; v++) for (int i = 0; i < nd; i++) if (ds[i] == sig[v]) { newc[v] = i; break; }
        /* graph-level hash: sorted multiset of signatures */
        uint64_t ms[MAXN]; for (int v = 0; v < n; v++) ms[v] = sig[v];
        for (int i = 1; i < n; i++) { uint64_t t = ms[i]; int j = i - 1; while (j >= 0 && ms[j] > t) { ms[j + 1] = ms[j]; j--; } ms[j + 1] = t; }
        for (int v = 0; v < n; v++) H = mix(H ^ (ms[v] + 0x632BE59BD9B4E019ULL * (round + 1)));
        int oldnd = 0; { int seen[MAXN * 4] = {0}; (void)seen; }
        /* count distinct old colours */
        int distinct_old = 0; { int c2[MAXN]; for (int v = 0; v < n; v++) { int f = 0; for (int i = 0; i < distinct_old; i++) if (c2[i] == col[v]) { f = 1; break; } if (!f) c2[distinct_old++] = col[v]; } }
        (void)oldnd;
        for (int v = 0; v < n; v++) { if (newc[v] != col[v]) changed = 1; col[v] = newc[v]; }
        if (nd == distinct_old && round > 0) break;
        (void)changed;
    }
    *hash_out = H;
}

/* ---------- exact isomorphism / automorphism search ---------- */
static int iso_n; static const graph *GA, *GB; static int ca[MAXN], cb[MAXN];
static int mapAB[MAXN], usedB;

static int iso_rec(int i) {
    if (i == iso_n) return 1;
    for (int c = 0; c < iso_n; c++) {
        if (usedB >> c & 1) continue;
        if (cb[c] != ca[i]) continue;
        int ok = 1;
        for (int j = 0; j < i; j++) {
            int ea = GA->adj[i] >> j & 1, eb = GB->adj[c] >> mapAB[j] & 1;
            if (ea != eb) { ok = 0; break; }
        }
        if (!ok) continue;
        mapAB[i] = c; usedB |= 1 << c;
        if (iso_rec(i + 1)) return 1;
        usedB &= ~(1 << c);
    }
    return 0;
}

static int isomorphic(const graph *a, const graph *b) {
    if (a->n != b->n || a->d != b->d || a->h != b->h) return 0;
    uint64_t h1, h2;
    wl_colors(a, ca, &h1); wl_colors(b, cb, &h2);
    iso_n = a->n; GA = a; GB = b; usedB = 0;
    return iso_rec(0);
}

/* automorphism search with fixed prefix: g(b_j) = b_j for j < i, g(i) = c */
static int aut_n; static const graph *GG; static int cc[MAXN]; static int amap[MAXN]; static int aused;
static int aut_rec(int v) {
    if (v == aut_n) return 1;
    if (amap[v] >= 0) {   /* preassigned */
        int c = amap[v];
        for (int j = 0; j < v; j++) if ((GG->adj[v] >> j & 1) != (GG->adj[c] >> amap[j] & 1)) return 0;
        return aut_rec(v + 1);
    }
    for (int c = 0; c < aut_n; c++) {
        if (aused >> c & 1) continue;
        if (cc[c] != cc[v]) continue;
        int ok = 1;
        for (int j = 0; j < v; j++) if ((GG->adj[v] >> j & 1) != (GG->adj[c] >> amap[j] & 1)) { ok = 0; break; }
        if (!ok) continue;
        amap[v] = c; aused |= 1 << c;
        /* also check consistency with preassigned later vertices */
        if (aut_rec(v + 1)) return 1;
        amap[v] = -1; aused &= ~(1 << c);
    }
    return 0;
}

static int sign_of(const graph *g, const int *p) {
    int par = 0;
    for (int u = 0; u < g->n; u++) { row r = g->adj[u] & ~((1u << (u + 1)) - 1); while (r) { int v = __builtin_ctz(r); r &= r - 1; if (p[u] > p[v]) par ^= 1; } }
    return par;
}

/* returns 1 if ODD (some automorphism has odd sign), 0 if EVEN */
static int royle_odd(const graph *g) {
    int n = g->n; uint64_t hh;
    wl_colors(g, cc, &hh);
    GG = g; aut_n = n;
    for (int i = 0; i < n; i++) {           /* level i: fix 0..i-1, move i */
        for (int c = 0; c < n; c++) {
            if (c == i) continue;
            if (cc[c] != cc[i]) continue;
            if (c < i) continue;            /* c < i is fixed by the stabiliser of 0..i-1 */
            for (int v = 0; v < n; v++) amap[v] = -1;
            aused = 0;
            for (int j = 0; j < i; j++) { amap[j] = j; aused |= 1 << j; }
            amap[i] = c; aused |= 1 << c;
            /* the recursion assigns vertices in order; preassigned ones are checked */
            if (aut_rec(0)) { if (sign_of(g, amap)) return 1; }
        }
    }
    return 0;
}

/* ---------- domination ---------- */
static uint64_t count_dom(const graph *g) {
    int n = g->n; uint32_t full = (1u << n) - 1; uint64_t c = 0;
    row cl[MAXN]; for (int v = 0; v < n; v++) cl[v] = g->adj[v] | (1u << v);
    for (uint32_t S = 0; S <= full; S++) {
        uint32_t dm = 0, s = S; while (s) { int v = __builtin_ctz(s); s &= s - 1; dm |= cl[v]; }
        if (dm == full) c++;
    }
    return c;
}

static uint64_t indep2(const graph *g) {   /* sum over independent sets of 2^|S| */
    int n = g->n; uint32_t full = (1u << n) - 1; uint64_t t = 0;
    for (uint32_t S = 0; S <= full; S++) {
        int ok = 1; uint32_t s = S; while (s) { int v = __builtin_ctz(s); s &= s - 1; if (g->adj[v] & S) { ok = 0; break; } }
        if (ok) t += 1ULL << __builtin_popcount(S);
    }
    return t;
}

/* ---------- level storage with hash buckets ---------- */
typedef struct { graph *g; size_t cnt, cap; int64_t *tab; size_t tcap; int64_t *next; } level;

static void lv_init(level *L) { L->cnt = 0; L->cap = 1024; L->g = malloc(L->cap * sizeof(graph)); L->tcap = 1 << 16; L->tab = malloc(L->tcap * sizeof(int64_t)); for (size_t i = 0; i < L->tcap; i++) L->tab[i] = -1; L->next = malloc(L->cap * sizeof(int64_t)); }

static void lv_rehash(level *L) {
    size_t nc = L->tcap * 4; int64_t *nt = malloc(nc * sizeof(int64_t)); for (size_t i = 0; i < nc; i++) nt[i] = -1;
    for (size_t i = 0; i < L->cnt; i++) { size_t b = L->g[i].h & (nc - 1); L->next[i] = nt[b]; nt[b] = i; }
    free(L->tab); L->tab = nt; L->tcap = nc;
}

static int lv_insert(level *L, const graph *g) {
    size_t b = g->h & (L->tcap - 1);
    for (int64_t i = L->tab[b]; i >= 0; i = L->next[i]) if (isomorphic(&L->g[i], g)) return 0;
    if (L->cnt == L->cap) { L->cap *= 2; L->g = realloc(L->g, L->cap * sizeof(graph)); L->next = realloc(L->next, L->cap * sizeof(int64_t)); }
    L->g[L->cnt] = *g; L->next[L->cnt] = L->tab[b]; L->tab[b] = L->cnt; L->cnt++;
    if (L->cnt > L->tcap) lv_rehash(L);
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: dom_even K maxn\n"); return 1; }
    K = strtoull(argv[1], 0, 10); int maxn = atoi(argv[2]); if (maxn > MAXN) maxn = MAXN;
    level cur, nxt; lv_init(&cur);
    graph g0; memset(&g0, 0, sizeof g0); g0.n = 1; g0.d = 1; { int c[MAXN]; wl_colors(&g0, c, &g0.h); }
    lv_insert(&cur, &g0);
    for (int n = 1; n <= maxn; n++) {
        /* report level n */
        long ne = 0, no = 0;
        for (size_t i = 0; i < cur.cnt; i++) {
            graph *g = &cur.g[i];
            int odd = royle_odd(g);
            if (count_dom(g) != g->d) { fprintf(stderr, "dom mismatch\n"); return 2; }
            if ((g->d & 1) == 0) { fprintf(stderr, "EVEN dominating count?!\n"); return 3; }
            int m = 0; for (int v = 0; v < n; v++) m += __builtin_popcount(g->adj[v]); m /= 2;
            printf("%d %d %llu %c %llu ", n, m, (unsigned long long)g->d, odd ? 'O' : 'E', (unsigned long long)indep2(g));
            for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if (g->adj[u] >> v & 1) printf("%x%x", u, v);
            printf("\n");
            if (odd) no++; else ne++;
        }
        fprintf(stderr, "n=%d classes=%zu even=%ld odd=%ld\n", n, cur.cnt, ne, no);
        fflush(stdout);
        if (n == maxn) break;
        lv_init(&nxt);
        for (size_t i = 0; i < cur.cnt; i++) {
            graph *g = &cur.g[i];
            uint32_t full = (1u << n) - 1, Ssz = 1u << n;
            uint64_t *f = calloc(Ssz, sizeof(uint64_t)), *cd = calloc(Ssz, sizeof(uint64_t));
            row cl[MAXN]; for (int v = 0; v < n; v++) cl[v] = g->adj[v] | (1u << v);
            for (uint32_t S = 0; S <= full; S++) {
                uint32_t dm = 0, s = S; while (s) { int v = __builtin_ctz(s); s &= s - 1; dm |= cl[v]; }
                f[full & ~dm]++;
                if (dm == full) cd[S]++;
            }
            for (int b = 0; b < n; b++) for (uint32_t S = 0; S <= full; S++) if (S >> b & 1) { f[S] += f[S ^ (1u << b)]; cd[S] += cd[S ^ (1u << b)]; }
            for (uint32_t N = 1; N <= full; N++) {
                uint64_t dn = f[N] + g->d - cd[full & ~N];
                if (dn > K) continue;
                graph h = *g; h.n = n + 1; h.adj[n] = (row)N;
                for (int v = 0; v < n; v++) if (N >> v & 1) h.adj[v] |= 1u << n;
                h.d = dn; { int c[MAXN]; wl_colors(&h, c, &h.h); }
                lv_insert(&nxt, &h);
            }
            free(f); free(cd);
        }
        free(cur.g); free(cur.tab); free(cur.next);
        cur = nxt;
        if (cur.cnt == 0) { fprintf(stderr, "no graphs at n=%d\n", n + 1); break; }
    }
    return 0;
}
