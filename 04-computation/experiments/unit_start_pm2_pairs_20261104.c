/*
 * unit_start_pm2_pairs_20261104.c  (opus-2026-11-04-S76)
 *
 * Generates every strong tournament with at most K Hamiltonian cycles, up to isomorphism, order by order, and tabulates
 *   t_k(n) = number of classes with exactly k Hamiltonian cycles,
 *   u_k(n) = number of those with no vertex of out-degree 1 ("sole-exit-free"; klein's HYP-9459),
 * split by arc-strong connectivity (lambda = 1, or lambda >= 2) and vertex connectivity (2-strong or a cut vertex).
 *
 * Completeness of the generation.
 * - Every strong T of order n >= 4 has a vertex v with T - v strong (Moon 1966).
 * - Every Hamiltonian cycle C' of T - v extends to T at each in->out switch of v along C' (at least one switch, since v
 *   has in- and out-neighbours). Different (C', switch) pairs give different cycles of T (shortcut v to recover C').
 * - Hence hc(T - v) <= hc(T), and every strong T with hc <= K is T' + v for a strong T' with hc <= K.
 * For a parent T' of order m and a new vertex with out-set O and in-set I = complement:
 *     hc(T' + v) = sum over w in O, u in I of P_T'(w -> u)   (Hamiltonian paths of T' from w to u),
 * which a Gray-code walk over all O evaluates in O(m) per subset.
 *
 * For every sole-exit-free class with lambda = 1 the program also records the minimal one-arc-leaving sets D (exactly
 * one arc u -> l leaves D) and checks hc(T) = hc(T/D) * e_u(T[D]), where e_u counts Hamiltonian paths of T[D] ending at u.
 *
 * Usage: fewcycle NMAX K [DUMPN]   (NMAX <= 15, K <= MAXK; classes of order <= DUMPN with no sole exit are printed)
 *   fewcycle 14 5      : the census of HYP-9461 (orders 3-14)
 *   fewcycle 11 40     : built with -DMAXK=60, the minimum number of Hamiltonian cycles of 2-strong tournaments
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define MAXN 16
#ifndef MAXK
#define MAXK 8   /* compile with -DMAXK=60 for the minimum-cycle table of 2-strong tournaments */
#endif
typedef unsigned long long u64;
typedef long long i64;
typedef unsigned __int128 u128;

/* ---------------- canonical form (same routine as uh_one_arc_cut_20261104.c) ---------------- */
typedef struct {
    int n; const uint32_t *out;
    int cellof[MAXN], cellverts[MAXN][MAXN], cellsize[MAXN];
    int used[MAXN], pos[MAXN];
    u128 best; int have;
} CCtx;

static void canon_rec(CCtx *C, int q, u128 code) {
    int n = C->n;
    if (q == n) { if (!C->have || code < C->best) { C->best = code; C->have = 1; } return; }
    int c = C->cellof[q];
    int nb = (q + 1) * q / 2, tot = n * (n - 1) / 2;
    for (int k = 0; k < C->cellsize[c]; k++) {
        int v = C->cellverts[c][k];
        if (C->used[v]) continue;
        u128 nc = code;
        for (int p = 0; p < q; p++) nc = (nc << 1) | ((C->out[C->pos[p]] >> v) & 1u);
        if (C->have && nc > (C->best >> (tot - nb))) continue;
        C->used[v] = 1; C->pos[q] = v;
        canon_rec(C, q + 1, nc);
        C->used[v] = 0;
    }
}

static u128 canon(int n, const uint32_t *out) {
    CCtx C; memset(&C, 0, sizeof C);
    int col[MAXN], ncol; long key[MAXN];
    for (int v = 0; v < n; v++) {
        int s = __builtin_popcount(out[v]);
        uint32_t inv = (~out[v]) & ((1u << n) - 1) & ~(1u << v);
        int c3 = 0;
        for (int u = 0; u < n; u++) if ((out[v] >> u) & 1) c3 += __builtin_popcount(out[u] & inv);
        key[v] = (long)s * 100000 + c3;
    }
    {
        long tmp[MAXN]; memcpy(tmp, key, sizeof(long) * n);
        for (int i = 1; i < n; i++) { long x = tmp[i]; int j = i - 1; while (j >= 0 && tmp[j] > x) { tmp[j + 1] = tmp[j]; j--; } tmp[j + 1] = x; }
        long d[MAXN]; int nd = 0;
        for (int i = 0; i < n; i++) if (i == 0 || tmp[i] != tmp[i - 1]) d[nd++] = tmp[i];
        for (int v = 0; v < n; v++) for (int i = 0; i < nd; i++) if (d[i] == key[v]) { col[v] = i; break; }
        ncol = nd;
    }
    for (;;) {
        int sig[MAXN][MAXN + 1];
        for (int v = 0; v < n; v++) {
            sig[v][0] = col[v];
            for (int c = 0; c < ncol; c++) sig[v][c + 1] = 0;
            for (int u = 0; u < n; u++) if ((out[v] >> u) & 1) sig[v][col[u] + 1]++;
        }
        int ord[MAXN]; for (int i = 0; i < n; i++) ord[i] = i;
        for (int i = 1; i < n; i++) {
            int x = ord[i], j = i - 1;
            while (j >= 0 && memcmp(sig[ord[j]], sig[x], sizeof(int) * (ncol + 1)) > 0) { ord[j + 1] = ord[j]; j--; }
            ord[j + 1] = x;
        }
        int ncol2 = 0, newcol[MAXN];
        for (int i = 0; i < n; i++) {
            if (i > 0 && memcmp(sig[ord[i]], sig[ord[i - 1]], sizeof(int) * (ncol + 1)) != 0) ncol2++;
            newcol[ord[i]] = ncol2;
        }
        ncol2++;
        memcpy(col, newcol, sizeof(int) * n);
        if (ncol2 == ncol) break;
        ncol = ncol2;
    }
    for (int c = 0; c < ncol; c++) C.cellsize[c] = 0;
    for (int v = 0; v < n; v++) { int c = col[v]; C.cellverts[c][C.cellsize[c]++] = v; }
    int p = 0;
    for (int c = 0; c < ncol; c++) for (int k = 0; k < C.cellsize[c]; k++) C.cellof[p++] = c;
    C.n = n; C.out = out;
    canon_rec(&C, 0, 0);
    return C.best;
}

static void decode(int n, u128 code, uint32_t *out) {
    int tot = n * (n - 1) / 2, b = tot - 1;
    for (int i = 0; i < n; i++) out[i] = 0;
    for (int q = 1; q < n; q++) for (int p = 0; p < q; p++) {
        int bit = (int)((code >> b) & 1); b--;
        if (bit) out[p] |= 1u << q; else out[q] |= 1u << p;
    }
}

/* ---------------- Hamiltonian path / cycle counts ---------------- */
/* all-pairs Hamiltonian path counts P[s*m + t] of the subtournament on the vertex list vs[0..m-1] */
static void all_pairs_sub(int m, const int *vs, const uint32_t *out, i64 *P, i64 *buf) {
    uint32_t loc[MAXN];
    for (int i = 0; i < m; i++) {
        loc[i] = 0;
        for (int j = 0; j < m; j++) if ((out[vs[i]] >> vs[j]) & 1) loc[i] |= 1u << j;
    }
    int full = (1 << m) - 1;
    for (int s = 0; s < m; s++) {
        i64 *dp = buf;
        memset(dp, 0, sizeof(i64) * (size_t)(full + 1) * m);
        dp[(size_t)(1 << s) * m + s] = 1;
        for (int S = 1; S <= full; S++) {
            if (!((S >> s) & 1)) continue;
            for (int v = 0; v < m; v++) {
                i64 c = dp[(size_t)S * m + v]; if (!c) continue;
                uint32_t nb = loc[v] & ~(uint32_t)S;
                while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << w)) * m + w] += c; }
            }
        }
        for (int t = 0; t < m; t++) P[s * m + t] = dp[(size_t)full * m + t];
    }
}

static i64 ham_cycles(int n, const uint32_t *out, i64 *buf) {
    int full = (1 << n) - 1;
    i64 *dp = buf;
    memset(dp, 0, sizeof(i64) * (size_t)(full + 1) * n);
    dp[(size_t)1 * n + 0] = 1;
    for (int S = 1; S <= full; S++) {
        if (!(S & 1)) continue;
        for (int v = 0; v < n; v++) {
            i64 c = dp[(size_t)S * n + v]; if (!c) continue;
            uint32_t nb = out[v] & ~(uint32_t)S;
            while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << w)) * n + w] += c; }
        }
    }
    i64 hc = 0;
    for (int v = 1; v < n; v++) if (out[v] & 1u) hc += dp[(size_t)full * n + v];
    return hc;
}

/* ---------------- connectivity ---------------- */
static int is_strong_mask(int n, const uint32_t *out, uint32_t V) {
    if (!V) return 1;
    int s = __builtin_ctz(V);
    uint32_t seen = 1u << s, fr = seen;
    while (fr) { int x = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = out[x] & V & ~seen; seen |= nb; fr |= nb; }
    if (seen != V) return 0;
    seen = 1u << s; fr = seen;
    while (fr) {
        int x = __builtin_ctz(fr); fr &= fr - 1;
        uint32_t nb = 0;
        for (uint32_t r = V & ~seen; r; r &= r - 1) { int y = __builtin_ctz(r); if ((out[y] >> x) & 1) nb |= 1u << y; }
        seen |= nb; fr |= nb;
    }
    return seen == V;
}

/* ---------------- hash set ---------------- */
typedef struct { u128 *key; int *hc; char *used; size_t cap, cnt; } HSet;
static void hs_init(HSet *h, size_t cap) {
    h->cap = cap; h->cnt = 0; h->key = malloc(sizeof(u128) * cap); h->hc = malloc(sizeof(int) * cap); h->used = calloc(cap, 1);
    if (!h->key || !h->hc || !h->used) { fprintf(stderr, "alloc\n"); exit(1); }
}
static int hs_insert(HSet *h, u128 k, int hc) {
    u64 x = (u64)k * 0x9E3779B97F4A7C15ULL ^ (u64)(k >> 64) * 0xC2B2AE3D27D4EB4FULL;
    size_t i = (size_t)(x >> 11) % h->cap;
    while (h->used[i]) { if (h->key[i] == k) return 0; i++; if (i == h->cap) i = 0; }
    h->used[i] = 1; h->key[i] = k; h->hc[i] = hc; h->cnt++; return 1;
}

typedef struct { u128 code; int hc; } Rec;

static u64 fib(int k) { u64 a = 0, b = 1; if (k < 0) return 0; for (int i = 0; i < k; i++) { u64 t = a + b; a = b; b = t; } return a; }


/*
 * All tournament classes of orders 3..NMAX by one-vertex extension (every out-set, canonical form, dedupe).
 * For every strong class and every root v: rooted canonical code (root individualised) and the profile
 * (h, s, q, F) of adaptive_metagraph_reformulation_20261103/ADAPTIVE.md.  Reports the number of rooted strong classes,
 * those with unit root start (s = 1), and the unordered pairs with s = s' = 1, F = F', |h - h'| = 2 (the +-2 port pairs:
 * Session189 wrappers H(D_m) = h + 2^m F then differ by exactly 2 at every clock m).
 */
static int strong_mask(int n, const uint32_t *out) {
    uint32_t full = (1u << n) - 1, seen = 1, fr = 1;
    while (fr) { int x = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = out[x] & ~seen; seen |= nb; fr |= nb; }
    if (seen != full) return 0;
    seen = 1; fr = 1;
    while (fr) { int x = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = 0;
        for (int y = 0; y < n; y++) if (!((seen >> y) & 1) && ((out[y] >> x) & 1)) nb |= 1u << y; seen |= nb; fr |= nb; }
    return seen == full;
}
/* rooted canonical code: relabel so that the root is vertex 0 and canonise the rest with the root fixed first */
static u128 rooted_code(int n, const uint32_t *out, int root) {
    /* move root to position 0, then canonise with a colouring that individualises it: append a new vertex trick is
       avoided; instead brute-force over the canonical routine with key forcing: we encode the root by adding a pendant
       structure-free marker: compute canon of the tournament where the root is replaced first and prefix the code by
       the root's out-row in the canonical order.  Simplest exact method: minimum over all labelings with root first,
       found by the refinement search with the root as its own first cell. */
    CCtx C; memset(&C, 0, sizeof C);
    int col[MAXN], ncol; long key[MAXN];
    for (int v = 0; v < n; v++) {
        int s = __builtin_popcount(out[v]);
        uint32_t inv = (~out[v]) & ((1u << n) - 1) & ~(1u << v);
        int c3 = 0;
        for (int u = 0; u < n; u++) if ((out[v] >> u) & 1) c3 += __builtin_popcount(out[u] & inv);
        key[v] = (v == root) ? -1 : (long)s * 100000 + c3;
    }
    {
        long tmp[MAXN]; memcpy(tmp, key, sizeof(long) * n);
        for (int i = 1; i < n; i++) { long x = tmp[i]; int j = i - 1; while (j >= 0 && tmp[j] > x) { tmp[j + 1] = tmp[j]; j--; } tmp[j + 1] = x; }
        long d[MAXN]; int nd = 0;
        for (int i = 0; i < n; i++) if (i == 0 || tmp[i] != tmp[i - 1]) d[nd++] = tmp[i];
        for (int v = 0; v < n; v++) for (int i = 0; i < nd; i++) if (d[i] == key[v]) { col[v] = i; break; }
        ncol = nd;
    }
    for (;;) {
        int sig[MAXN][MAXN + 1];
        for (int v = 0; v < n; v++) { sig[v][0] = col[v]; for (int c = 0; c < ncol; c++) sig[v][c + 1] = 0;
            for (int u = 0; u < n; u++) if ((out[v] >> u) & 1) sig[v][col[u] + 1]++; }
        int ord[MAXN]; for (int i = 0; i < n; i++) ord[i] = i;
        for (int i = 1; i < n; i++) { int x = ord[i], j = i - 1;
            while (j >= 0 && memcmp(sig[ord[j]], sig[x], sizeof(int) * (ncol + 1)) > 0) { ord[j + 1] = ord[j]; j--; } ord[j + 1] = x; }
        int ncol2 = 0, newcol[MAXN];
        for (int i = 0; i < n; i++) { if (i > 0 && memcmp(sig[ord[i]], sig[ord[i - 1]], sizeof(int) * (ncol + 1)) != 0) ncol2++; newcol[ord[i]] = ncol2; }
        ncol2++; memcpy(col, newcol, sizeof(int) * n);
        if (ncol2 == ncol) break; ncol = ncol2;
    }
    for (int c = 0; c < ncol; c++) C.cellsize[c] = 0;
    for (int v = 0; v < n; v++) { int c = col[v]; C.cellverts[c][C.cellsize[c]++] = v; }
    int p = 0; for (int c = 0; c < ncol; c++) for (int k = 0; k < C.cellsize[c]; k++) C.cellof[p++] = c;
    C.n = n; C.out = out; canon_rec(&C, 0, 0);
    return C.best;
}
typedef struct { long long h, s, q, F; } Prof;
static Prof rprof(int n, const uint32_t *out, int v, long long *dpe, long long *Hs) {
    int full = (1 << n) - 1;
    for (size_t i = 0; i < (size_t)(full + 1) * n; i++) dpe[i] = 0;
    for (int x = 0; x < n; x++) dpe[(size_t)(1 << x) * n + x] = 1;
    for (int S = 1; S <= full; S++) for (int x = 0; x < n; x++) { long long c = dpe[(size_t)S * n + x]; if (!c) continue;
        uint32_t nb = out[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dpe[(size_t)(S | (1 << y)) * n + y] += c; } }
    Hs[0] = 1; for (int S = 1; S <= full; S++) { long long t = 0; for (int x = 0; x < n; x++) t += dpe[(size_t)S * n + x]; Hs[S] = t; }
    Prof pr; pr.h = Hs[full]; pr.q = Hs[full & ~(1 << v)];
    long long F = 0; for (int S = 1; S <= full; S++) if ((S >> v) & 1) F += dpe[(size_t)S * n + v] * Hs[full & ~S];
    pr.F = F;
    uint32_t rev[MAXN]; for (int a = 0; a < n; a++) { rev[a] = 0; for (int b = 0; b < n; b++) if (b != a && !((out[a] >> b) & 1)) rev[a] |= 1u << b; }
    for (size_t i = 0; i < (size_t)(full + 1) * n; i++) dpe[i] = 0;
    for (int x = 0; x < n; x++) dpe[(size_t)(1 << x) * n + x] = 1;
    for (int S = 1; S <= full; S++) for (int x = 0; x < n; x++) { long long c = dpe[(size_t)S * n + x]; if (!c) continue;
        uint32_t nb = rev[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dpe[(size_t)(S | (1 << y)) * n + y] += c; } }
    pr.s = dpe[(size_t)full * n + v];
    return pr;
}
typedef struct { u128 code; Prof p; } RRec;
static int cmp_unit(const void *a, const void *b) {
    const RRec *x = a, *y = b;
    if (x->p.F != y->p.F) return x->p.F < y->p.F ? -1 : 1;
    if (x->p.h != y->p.h) return x->p.h < y->p.h ? -1 : 1;
    return x->code < y->code ? -1 : (x->code > y->code);
}

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 9;
    printf("# unit_start_pm2_pairs_20261104 NMAX=%d\n", NMAX);
    u128 *cur = malloc(sizeof(u128)); size_t ncur = 1;
    { uint32_t o[3] = {2u, 4u, 1u}; uint32_t t[3] = {6u, 4u, 0u}; (void)o; cur[0] = canon(3, t); }
    /* order 3: both classes (transitive and cyclic) */
    { u128 *c3 = malloc(sizeof(u128) * 2); uint32_t tt[3] = {6u, 4u, 0u}, cy[3] = {2u, 4u, 1u}; c3[0] = canon(3, tt); c3[1] = canon(3, cy); free(cur); cur = c3; ncur = 2; }
    for (int n = 3; n <= NMAX; n++) {
        long long *dpe = malloc(sizeof(long long) * ((size_t)1 << n) * n), *Hs = malloc(sizeof(long long) * ((size_t)1 << n));
        RRec *units = malloc(sizeof(RRec) * (ncur * n + 1)); size_t nu = 0; size_t nstrong = 0, nrooted = 0;
        HSet rs; hs_init(&rs, (ncur * n + 1024) * 4);
        for (size_t i = 0; i < ncur; i++) {
            uint32_t out[MAXN]; decode(n, cur[i], out);
            if (!strong_mask(n, out)) continue;
            nstrong++;
            for (int v = 0; v < n; v++) {
                u128 rc = rooted_code(n, out, v);
                if (!hs_insert(&rs, rc, 0)) continue;
                nrooted++;
                Prof pr = rprof(n, out, v, dpe, Hs);
                if (pr.s == 1) { units[nu].code = rc; units[nu].p = pr; nu++; }
            }
        }
        qsort(units, nu, sizeof(RRec), cmp_unit);
        long pairs = 0; long long exF = -1, exh = -1;
        for (size_t a = 0; a < nu; a++) for (size_t b = a + 1; b < nu && units[b].p.F == units[a].p.F; b++)
            if (units[b].p.h - units[a].p.h == 2) { pairs++; if (exF < 0) { exF = units[a].p.F; exh = units[a].p.h; } }
        printf("n=%d classes=%zu strong=%zu rooted strong=%zu unit-start rooted=%zu | unit +-2 pairs (s=s'=1, F=F', h'=h+2): %ld",
               n, ncur, nstrong, nrooted, nu, pairs);
        if (pairs) printf("  first: F=%lld h=%lld,%lld", exF, exh, exh + 2);
        printf("\n"); fflush(stdout);
        free(units); free(dpe); free(Hs); free(rs.key); free(rs.hc); free(rs.used);
        if (n == NMAX) break;
        /* extend all classes */
        HSet hs; hs_init(&hs, ncur * (1u << n) / 2 + 1024);
        u128 *nxt = malloc(sizeof(u128) * (ncur * (1u << n) / 4 + 1024)); size_t nn = 0, capn = ncur * (1u << n) / 4 + 1024;
        for (size_t i = 0; i < ncur; i++) {
            uint32_t out[MAXN]; decode(n, cur[i], out);
            for (uint32_t O = 0; O < (1u << n); O++) {
                uint32_t o2[MAXN]; for (int v = 0; v < n; v++) o2[v] = out[v] | (((O >> v) & 1) ? 0u : (1u << n)); o2[n] = O;
                u128 c = canon(n + 1, o2);
                if (hs_insert(&hs, c, 0)) { if (nn == capn) { capn *= 2; nxt = realloc(nxt, sizeof(u128) * capn); } nxt[nn++] = c; }
            }
        }
        free(hs.key); free(hs.hc); free(hs.used); free(cur); cur = nxt; ncur = nn;
    }
    return 0;
}
