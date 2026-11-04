/* uh_hpath_values_20261104.c (kind-pasteur-2026-11-04-S10): the generator of uh_noncritical_run_20261104.c with the
 * per-class checks replaced by an exact Hamiltonian-path count; prints the H values of every uniquely Hamiltonian
 * class of orders 4..NMAX (all asymmetric by THM-4868 (E)). Usage: uh_hpath_values_20261104 NMAX 0 0 */
/* uh_noncritical_run_20261104.c -- the non-critical run of a uniquely Hamiltonian tournament (exactly one directed
 * Hamiltonian cycle), checked exhaustively by one-vertex extension.
 *
 * GENERATION (complete; each step proved in the note):
 *   (M) a strong tournament of order >= 4 has a vertex v with T - v strong (Moon 1966, pancyclicity);
 *   (H) T - v strong  =>  every Hamiltonian cycle of T - v extends to one of T, so T - v is again
 *       uniquely Hamiltonian;
 *   (I) the in-neighbourhood I of v is one cyclic interval of the cycle C' of T - v, and
 *       hc(T) = sum_{w in O, u in I} hp_{T-v}(w, u)   (hp = Hamiltonian paths w -> u, capped at 2).
 *   Classes are stored in cycle coordinates (vertex i -> i+1 mod n is the cycle); the canonical form
 *   is the least adjacency over the n rotations (an isomorphism preserves the unique cycle).
 *
 * CHECKS per class (order n >= 4):
 *   round  : every non-critical v beats exactly the next d+(v) vertices of C;
 *   run    : the non-critical set NC is ONE cyclic run, 2 <= |NC| <= n-1, its first vertex has
 *            in-degree 1 and its last vertex has out-degree 1;
 *   crit   : out-degree-1 x => x+1 critical; in-degree-1 y => y-1 critical;
 *   cut    : the arc i -> i+1 is a cut arc (T minus it is not strong) iff not (i, i+1 both in NC);
 *            hence #cut arcs = n + 1 - |NC|;
 *   deg1   : #out1, #in1 in {1, 2};
 *   asym   : no nontrivial rotation is an automorphism;
 *   cut1   : a one-arc cut exists (interval form, plus brute force over all subsets for n <= BF).
 * TABLES: (#out1,#in1) counts, N(n, r) = #classes with |NC| = r, #cut-arc distribution.
 *
 * Usage: uh_noncritical_run_20261104 NMAX [XC_upto=12] [BF_upto=12]     (OpenMP over parent classes)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXN 22
typedef uint32_t row_t;
typedef struct { int n; row_t adj[MAXN]; } tour_t;
static int popc(uint32_t x) { return __builtin_popcount(x); }
static uint32_t fullmask(int n) { return (n >= 32) ? 0xffffffffu : ((1u << n) - 1u); }

static void rotate(const tour_t *t, int r, row_t *out) {
    int n = t->n;
    for (int i = 0; i < n; i++) {
        row_t src = t->adj[(i + r) % n], dst = 0;
        for (int j = 0; j < n; j++) if (src >> ((j + r) % n) & 1u) dst |= 1u << j;
        out[i] = dst;
    }
}
static int canon(const tour_t *t, tour_t *c) {   /* returns number of rotations fixing t (incl. identity) */
    int n = t->n, fix = 1; row_t best[MAXN], cur[MAXN];
    rotate(t, 0, best);
    for (int r = 1; r < n; r++) {
        rotate(t, r, cur);
        if (!memcmp(cur, t->adj, sizeof(row_t) * n)) fix++;
        if (memcmp(cur, best, sizeof(row_t) * n) < 0) memcpy(best, cur, sizeof(row_t) * n);
    }
    c->n = n; memset(c->adj, 0, sizeof c->adj); memcpy(c->adj, best, sizeof(row_t) * n);
    return fix;
}

/* ---------- hash set ---------- */
typedef struct { tour_t *v; size_t cnt, cap; int64_t *tab; size_t tcap; } set_t;
static uint64_t hsh(const tour_t *t) {
    uint64_t h = 1469598103934665603ULL ^ (uint64_t)t->n;
    for (int i = 0; i < t->n; i++) { h ^= t->adj[i]; h *= 1099511628211ULL; h ^= h >> 29; }
    return h;
}
static void set_init(set_t *s) { s->cnt = 0; s->cap = 1024; s->v = malloc(s->cap * sizeof(tour_t));
    s->tcap = 1u << 12; s->tab = malloc(s->tcap * sizeof(int64_t)); for (size_t i = 0; i < s->tcap; i++) s->tab[i] = -1; }
static void set_grow(set_t *s) {
    size_t nc = s->tcap * 2; int64_t *nt = malloc(nc * sizeof(int64_t));
    for (size_t i = 0; i < nc; i++) nt[i] = -1;
    for (size_t k = 0; k < s->cnt; k++) { size_t p = hsh(&s->v[k]) & (nc - 1); while (nt[p] >= 0) p = (p + 1) & (nc - 1); nt[p] = (int64_t)k; }
    free(s->tab); s->tab = nt; s->tcap = nc;
}
static int set_add(set_t *s, const tour_t *t) {
    if (2 * (s->cnt + 1) > s->tcap) set_grow(s);
    size_t p = hsh(t) & (s->tcap - 1);
    while (s->tab[p] >= 0) { const tour_t *u = &s->v[s->tab[p]];
        if (u->n == t->n && !memcmp(u->adj, t->adj, sizeof(row_t) * t->n)) return 0; p = (p + 1) & (s->tcap - 1); }
    if (s->cnt == s->cap) { s->cap *= 2; s->v = realloc(s->v, s->cap * sizeof(tour_t)); }
    s->v[s->cnt] = *t; s->tab[p] = (int64_t)s->cnt; s->cnt++; return 1;
}

/* ---------- Hamiltonian path endpoint counts capped at 2 ---------- */
static void hp_matrix(const tour_t *t, uint32_t *ge1, uint32_t *ge2, uint32_t *b1, uint32_t *b2) {
    int n = t->n; size_t M = (size_t)1 << n; uint32_t full = fullmask(n);
    memset(b1, 0, M * (size_t)n * sizeof(uint32_t)); memset(b2, 0, M * (size_t)n * sizeof(uint32_t));
    for (int s = 0; s < n; s++) b1[((size_t)1 << s) * n + s] = 1u << s;
    for (size_t mask = 1; mask < M; mask++) {
        for (int e = 0; e < n; e++) {
            uint32_t a1 = b1[mask * n + e]; if (!a1) continue;
            uint32_t a2 = b2[mask * n + e];
            uint32_t nxt = t->adj[e] & ~(uint32_t)mask & full;
            while (nxt) { int w = __builtin_ctz(nxt); nxt &= nxt - 1;
                size_t idx = (mask | ((size_t)1 << w)) * n + w;
                uint32_t c1 = b1[idx], c2 = b2[idx];
                b2[idx] = c2 | a2 | (c1 & a1);
                b1[idx] = c1 | a1; }
        }
    }
    for (int e = 0; e < n; e++) { ge1[e] = b1[(M - 1) * n + e]; ge2[e] = b2[(M - 1) * n + e]; }
}
/* independent exact Hamiltonian cycle count (64-bit) */
static uint64_t hc_exact(const tour_t *t) {
    int n = t->n; size_t M = (size_t)1 << n; uint64_t *dp = calloc(M * n, sizeof(uint64_t));
    dp[1 * n + 0] = 1;
    for (size_t mask = 1; mask < M; mask += 2) for (int e = 0; e < n; e++) { uint64_t c = dp[mask * n + e]; if (!c) continue;
        uint32_t nxt = t->adj[e] & ~(uint32_t)mask & (uint32_t)(M - 1);
        while (nxt) { int w = __builtin_ctz(nxt); nxt &= nxt - 1; dp[(mask | ((size_t)1 << w)) * n + w] += c; } }
    uint64_t tot = 0; for (int e = 1; e < n; e++) if (t->adj[e] & 1u) tot += dp[(M - 1) * n + e];
    free(dp); return tot;
}

/* strong connectivity of the digraph adj restricted to vertex set S (bitmask) */
static int strong_on(const row_t *adj, int n, uint32_t S) {
    if (popc(S) <= 1) return 1;
    int s0 = __builtin_ctz(S);
    uint32_t seen = 1u << s0, frontier = seen;
    while (frontier) { uint32_t nf = 0; uint32_t f = frontier;
        while (f) { int x = __builtin_ctz(f); f &= f - 1; nf |= adj[x] & S; }
        nf &= ~seen; seen |= nf; frontier = nf; }
    if (seen != S) return 0;
    /* backward */
    row_t radj[MAXN]; memset(radj, 0, sizeof radj);
    for (int x = 0; x < n; x++) if (S >> x & 1u) { uint32_t o = adj[x] & S; while (o) { int y = __builtin_ctz(o); o &= o - 1; radj[y] |= 1u << x; } }
    seen = 1u << s0; frontier = seen;
    while (frontier) { uint32_t nf = 0; uint32_t f = frontier;
        while (f) { int x = __builtin_ctz(f); f &= f - 1; nf |= radj[x]; }
        nf &= ~seen; seen |= nf; frontier = nf; }
    return seen == S;
}
static int back_arcs(const tour_t *t, uint32_t A) {
    int c = 0; for (int x = 0; x < t->n; x++) if (!(A >> x & 1u)) c += popc(t->adj[x] & A); return c;
}

typedef struct {
    long cls, asym_bad, round_bad, run_bad, crit_bad, cut_bad, deg1_bad, nocut1, xc_bad, bf_bad;
    long oi[3][3];           /* (#out1, #in1) */
    long Nr[MAXN + 1];       /* |NC| = r */
    long Kc[MAXN + 1];       /* #cut arcs = k */
} stats_t;

static void check_class(const tour_t *t, stats_t *st, int XC, int BF) {
    int n = t->n; uint32_t F = fullmask(n);
    st->cls++;
    { tour_t c; if (canon(t, &c) != 1) st->asym_bad++; }
    int od[MAXN], id[MAXN]; int out1 = 0, in1 = 0;
    for (int i = 0; i < n; i++) { od[i] = popc(t->adj[i]); id[i] = n - 1 - od[i]; if (od[i] == 1) out1++; if (id[i] == 1) in1++; }
    if (out1 < 1 || out1 > 2 || in1 < 1 || in1 > 2) st->deg1_bad++;
    else st->oi[out1][in1]++;
    uint32_t NC = 0;
    for (int v = 0; v < n; v++) if (strong_on(t->adj, n, F & ~(1u << v))) NC |= 1u << v;
    int r = popc(NC); st->Nr[r]++;
    /* round */
    for (int v = 0; v < n; v++) if (NC >> v & 1u) {
        uint32_t want = 0; for (int s = 1; s <= od[v]; s++) want |= 1u << ((v + s) % n);
        if (t->adj[v] != want) { st->round_bad++; break; }
    }
    /* run: exactly one cyclic run, 2 <= r <= n-1, first in-degree 1, last out-degree 1 */
    {
        int starts = 0, ok = 1, first = -1, last = -1;
        for (int i = 0; i < n; i++) {
            int inNC = NC >> i & 1u, prevNC = NC >> ((i + n - 1) % n) & 1u, nextNC = NC >> ((i + 1) % n) & 1u;
            if (inNC && !prevNC) { starts++; first = i; }
            if (inNC && !nextNC) last = i;
        }
        if (starts != 1 || r < 2 || r > n - 1) ok = 0;
        if (ok && (id[first] != 1 || od[last] != 1)) ok = 0;
        if (!ok) st->run_bad++;
    }
    /* crit */
    for (int i = 0; i < n; i++) {
        if (od[i] == 1 && (NC >> ((i + 1) % n) & 1u)) { st->crit_bad++; break; }
        if (id[i] == 1 && (NC >> ((i + n - 1) % n) & 1u)) { st->crit_bad++; break; }
    }
    /* cut arcs */
    {
        int k = 0, bad = 0;
        for (int i = 0; i < n; i++) {
            row_t a2[MAXN]; memcpy(a2, t->adj, sizeof(row_t) * n);
            a2[i] &= ~(1u << ((i + 1) % n));
            int cut = !strong_on(a2, n, F);
            k += cut;
            int pred = !((NC >> i & 1u) && (NC >> ((i + 1) % n) & 1u));
            if (cut != pred) bad = 1;
        }
        if (bad || k != n + 1 - r) st->cut_bad++;
        st->Kc[k]++;
    }
    /* one-arc cut, interval form and brute force */
    {
        int found = 0;
        for (int s = 0; s < n && !found; s++) for (int len = 1; len < n && !found; len++) {
            uint32_t A = 0; for (int q = 0; q < len; q++) A |= 1u << ((s + q) % n);
            if (back_arcs(t, A) == 1) found = 1;
        }
        if (!found) st->nocut1++;
        if (n <= BF) { int bf = 0; for (uint32_t A = 1; A < F && !bf; A++) if (back_arcs(t, A) == 1) bf = 1; if (bf != found) st->bf_bad++; }
    }
    if (n <= XC && hc_exact(t) != 1) st->xc_bad++;
}

static uint64_t fib(int k) { if (k < 0) return (k % 2 == 0) ? (uint64_t)(-(int64_t)fib(-k)) : fib(-k); uint64_t a = 0, b = 1; for (int i = 0; i < k; i++) { uint64_t c = a + b; a = b; b = c; } return a; }

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 12;
    int XC = argc > 2 ? atoi(argv[2]) : 12;
    int BF = argc > 3 ? atoi(argv[3]) : 12;
    if (NMAX > 20) { fprintf(stderr, "NMAX too large\n"); return 1; }
    set_t lev[MAXN + 1];
    set_init(&lev[3]);
    { tour_t c3 = { .n = 3 }; c3.adj[0] = 2; c3.adj[1] = 4; c3.adj[2] = 1; tour_t c; canon(&c3, &c); set_add(&lev[3], &c); }
    printf("# uh_noncritical_run_20261104 NMAX=%d XC<=%d BF<=%d\n", NMAX, XC, BF);
    for (int n = 4; n <= NMAX; n++) {
        set_t *S = &lev[n]; set_init(S);
        set_t *P = &lev[n - 1]; int np = n - 1; size_t M = (size_t)1 << np;
        long nparent = (long)P->cnt;
        #pragma omp parallel
        {
            uint32_t *b1 = malloc(M * np * sizeof(uint32_t)), *b2 = malloc(M * np * sizeof(uint32_t));
            #pragma omp for schedule(dynamic, 16)
            for (long k = 0; k < nparent; k++) {
                const tour_t *tp = &P->v[k];
                uint32_t ge1[MAXN], ge2[MAXN];
                hp_matrix(tp, ge1, ge2, b1, b2);
                tour_t local[MAXN * MAXN]; int nl = 0;
                for (int p = 0; p < np; p++) for (int L = 1; L < np; L++) {
                    uint32_t I = 0; for (int q = 0; q < L; q++) I |= 1u << ((p + q) % np);
                    uint32_t O = fullmask(np) & ~I;
                    int tot = 0;
                    for (int u = 0; u < np && tot < 2; u++) if (I >> u & 1u) tot += popc(ge1[u] & O) + popc(ge2[u] & O);
                    if (tot != 1) continue;
                    int q = (p + L - 1) % np;
                    tour_t t = { .n = n }; memset(t.adj, 0, sizeof t.adj);
                    int nid[MAXN]; for (int x = 0; x < np; x++) nid[x] = (x <= q) ? x : x + 1;
                    int vid = q + 1;
                    for (int x = 0; x < np; x++) { row_t rr = 0; for (int y = 0; y < np; y++) if (tp->adj[x] >> y & 1u) rr |= 1u << nid[y];
                        if (I >> x & 1u) rr |= 1u << vid; t.adj[nid[x]] = rr; }
                    { row_t rr = 0; for (int y = 0; y < np; y++) if (O >> y & 1u) rr |= 1u << nid[y]; t.adj[vid] = rr; }
                    canon(&t, &local[nl++]);
                }
                #pragma omp critical(setadd)
                { for (int j = 0; j < nl; j++) set_add(S, &local[j]); }
            }
            free(b1); free(b2);
        }
        free(P->v); free(P->tab); P->v = NULL; P->tab = NULL;   /* parent level no longer needed */

        {
            long ncls = (long)S->cnt;
            uint64_t *Hs = malloc(sizeof(uint64_t) * ncls);
            #pragma omp parallel for schedule(dynamic, 16)
            for (long k = 0; k < ncls; k++) {
                const tour_t *t = &S->v[k]; int nn = t->n; size_t M = (size_t)1 << nn;
                uint64_t *dp = calloc(M * nn, sizeof(uint64_t));
                for (int s = 0; s < nn; s++) dp[((size_t)1 << s) * nn + s] = 1;
                for (size_t mask = 1; mask < M; mask++) for (int e = 0; e < nn; e++) { uint64_t c = dp[mask * nn + e]; if (!c) continue;
                    uint32_t nx = t->adj[e] & ~(uint32_t)mask & (uint32_t)(M - 1);
                    while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dp[(mask | ((size_t)1 << w)) * nn + w] += c; } }
                uint64_t tot = 0; for (int e = 0; e < nn; e++) tot += dp[(M - 1) * nn + e];
                Hs[k] = tot; free(dp);
            }
            printf("UH n=%d classes=%ld H:", n, ncls);
            for (long k = 0; k < ncls; k++) printf(" %llu", (unsigned long long)Hs[k]);
            printf("\n");
            free(Hs);
        }
        fflush(stdout);
    }
    return 0;
}
