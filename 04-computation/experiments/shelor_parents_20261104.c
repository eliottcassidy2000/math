/* shelor_parents_20261104.c -- THM-4868 (G) and (H) checked on every uniquely Hamiltonian class.
 * Classes are generated exactly as in uh_noncritical_run_20261104.c (interval extension; Garey-checked).
 * For each class U (cycle coordinates) with run start y = y_U:
 *   (G)  adm(U) = {p : e_p(U) = 1} equals {y, y-} plus y-- iff y- is a sole entry;
 *   (H)  the child U^p (new vertex beating all but p, p -> new) has run length r(U), r(U)+1, r(U)+2
 *        for p = y, y-, y-- respectively.
 * Usage: shelor_parents_20261104 NMAX      (checks parents of orders 4..NMAX) */
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

static uint32_t ncset(const uint32_t *adj, int n) {
    uint32_t F = (1u << n) - 1, NC = 0;
    for (int v = 0; v < n; v++) if (strong_on(adj, n, F & ~(1u << v))) NC |= 1u << v;
    return NC;
}
static void endcounts(const uint32_t *adj, int n, uint64_t *e) {  /* e[p] = Ham paths ending at p (capped 1e18) */
    size_t M = (size_t)1 << n; uint64_t *dp = calloc(M * n, sizeof(uint64_t));
    for (int s = 0; s < n; s++) dp[((size_t)1 << s) * n + s] = 1;
    for (size_t m = 1; m < M; m++) for (int a = 0; a < n; a++) { uint64_t c = dp[m * n + a]; if (!c) continue;
        uint32_t nx = adj[a] & ~(uint32_t)m & (uint32_t)(M - 1);
        while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dp[(m | ((size_t)1 << w)) * n + w] += c; } }
    for (int p = 0; p < n; p++) e[p] = dp[(M - 1) * n + p];
    free(dp);
}

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 12;
    if (NMAX > 19) { fprintf(stderr, "NMAX too large\n"); return 1; }
    set_t lev[MAXN + 1];
    set_init(&lev[3]);
    { tour_t c3 = { .n = 3 }; c3.adj[0] = 2; c3.adj[1] = 4; c3.adj[2] = 1; tour_t c; canon(&c3, &c); set_add(&lev[3], &c); }
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
        free(P->v); free(P->tab); P->v = NULL; P->tab = NULL;
        /* (G) and (H) on every class of order n (as a parent) */
        long cls = (long)S->cnt, gbad = 0, hbad = 0, E2 = 0, E3 = 0;
        #pragma omp parallel for schedule(dynamic, 16) reduction(+:gbad,hbad,E2,E3)
        for (long k = 0; k < cls; k++) {
            const uint32_t *adj = S->v[k].adj;
            uint32_t NC = ncset(adj, n); int r = popc(NC), y = -1;
            for (int v = 0; v < n; v++) if ((NC >> v & 1u) && !(NC >> ((v + n - 1) % n) & 1u)) y = v;
            int ym = (y + n - 1) % n, ymm = (y + n - 2) % n;
            int ym_entry = (popc(adj[ym]) == n - 2);
            uint64_t e[MAXN]; endcounts(adj, n, e);
            uint32_t want = (1u << y) | (1u << ym) | (ym_entry ? (1u << ymm) : 0u), adm = 0;
            for (int p = 0; p < n; p++) if (e[p] == 1) adm |= 1u << p;
            if (adm != want) gbad++;
            if (popc(adm) == 2) E2++; else if (popc(adm) == 3) E3++;
            for (int p = 0; p < n; p++) if (adm >> p & 1u) {
                uint32_t ch[MAXN]; memcpy(ch, adj, sizeof(uint32_t) * n);
                ch[p] |= 1u << n; ch[n] = ((1u << n) - 1) & ~(1u << p);
                int d = popc(ncset(ch, n + 1)) - r, dw = (p == y) ? 0 : (p == ym) ? 1 : 2;
                if (d != dw) hbad++;
            }
        }
        printf("parent order %2d: classes=%ld  |adm|=2:%ld |adm|=3:%ld  G_bad=%ld  H_bad=%ld\n", n, cls, E2, E3, gbad, hbad);
        fflush(stdout);
    }
    return 0;
}
