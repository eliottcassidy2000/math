/* classes_audit2.c -- independent audit (audit2, 2026-11-04) of THM-4868 (G) and (H), class level.
 *
 * Written from scratch for the audit; it does not reuse the author's programs or outputs.
 *
 * Generation of the uniquely Hamiltonian (UH) classes, order by order, by TWO methods whose
 * completeness does not depend on the statements being audited:
 *   order 4..LAB : exhaustive over all 2^C(n,2) labelled tournaments (exact DP count hc = 1);
 *   method A     : "Shelor child" U^p for every class U of order m and EVERY vertex p, kept when an
 *                  exact DP Hamiltonian-cycle count of the (m+1)-vertex tournament is 1;
 *   method B     : "Moon extension": U + v with an ARBITRARY out-set O (2^m - 2 choices), kept when
 *                  hc = sum_{w in O, u not in O} hp_U(w,u) = 1, hp from an independent start-fixed DP.
 *                  Complete by Moon 1966 (a strong tournament of order >= 4 has a strong
 *                  subtournament of order one less; T - v strong and hc(T) = 1 force hc(T - v) = 1,
 *                  since each Hamiltonian cycle of T - v extends to one of T and distinct cycles give
 *                  distinct extensions).
 * Classes are identified by the minimum over the n rotations of the adjacency word read along the
 * unique Hamiltonian cycle (a complete invariant for UH tournaments: isomorphisms preserve the
 * unique cycle).  At orders <= 9 every class is also given a brute-force canonical form over all n!
 * permutations, as a sanity check of the rotation invariant.
 *
 * Checks on every class of order 4..NCHECK (children of order NCHECK + 1 must be generated):
 *   hc = 1 by exact DP; asymmetry (the n rotation words are distinct);
 *   (C-input) NC one run from the unique non-critical sole entry y to a sole exit x;
 *   (G) adm = {p : e_p = 1} = {y, y-} + (y-- iff y- sole entry);
 *   (L3) hc(T^p) (exact DP, from method A) = e_p(T) for every p;
 *   (H) r(T^p) - r(T) = 0/1/2 for p = y, y-, y--; the four facts used in the proof;
 *   bijection: every class of order n+1 is hit by exactly one admissible pair (T, p), and
 *              T^p - y_child = T with y_child- = p;
 *   decomposition (n >= 5): with U = T - y, p0 = y-, for every p not in {y, p0}:
 *              e_p(T) = e_p(U) - hp_U(p0,p) + W(p0,p),  W = sum_S E_S(p0) E_{V-S}(p),
 *              and e_y(T) = e_{p0}(U) = 1 = e_{p0}(T);  min W, min (W - hp) recorded;
 *   cut arcs: k = #arcs e of C with T - e not strong;  k + r = n + 1;  T(n,k) table.
 *
 * Usage: classes_audit2 NMAX NCHECK [LAB] [VERIFYB]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define MAXV 16
typedef uint64_t u64;
typedef unsigned __int128 u128;
typedef struct { int n; uint32_t out[MAXV], in[MAXV]; } Tour;

static void finish(Tour *T) {
    uint32_t full = (1u << T->n) - 1;
    for (int v = 0; v < T->n; v++) T->in[v] = full & ~T->out[v] & ~(1u << v);
}
static inline int arc(const Tour *T, int u, int v) { return (T->out[u] >> v) & 1; }

/* ---------- counting ---------- */
static u64 hc_exact(const Tour *T, u64 *f) {      /* f: >= 2^(n-1) * n */
    int n = T->n; size_t N = (size_t)1 << (n - 1);
    memset(f, 0, sizeof(u64) * N * n);
    for (int v = 1; v < n; v++) if (arc(T, 0, v)) f[((size_t)1 << (v - 1)) * n + v] = 1;
    for (size_t S = 1; S < N; S++)
        for (int v = 1; v < n; v++) {
            if (!((S >> (v - 1)) & 1)) continue;
            u64 val = f[S * n + v]; if (!val) continue;
            uint32_t cand = (T->out[v] >> 1) & ~(uint32_t)S;
            while (cand) { int b = __builtin_ctz(cand); cand &= cand - 1; f[(S | ((size_t)1 << b)) * n + b + 1] += val; }
        }
    u64 tot = 0;
    for (int v = 1; v < n; v++) if (arc(T, v, 0)) tot += f[(N - 1) * n + v];
    return tot;
}
/* f[S*n+v] = number of paths covering exactly S and ending at v (any start); single vertex = 1 */
static void end_table(const Tour *T, u64 *f) {
    int n = T->n; size_t N = (size_t)1 << n;
    memset(f, 0, sizeof(u64) * N * n);
    for (int v = 0; v < n; v++) f[((size_t)1 << v) * n + v] = 1;
    for (size_t S = 1; S < N; S++)
        for (int v = 0; v < n; v++) {
            if (!((S >> v) & 1)) continue;
            u64 val = f[S * n + v]; if (!val) continue;
            uint32_t cand = T->out[v] & ~(uint32_t)S;
            while (cand) { int w = __builtin_ctz(cand); cand &= cand - 1; f[(S | ((size_t)1 << w)) * n + w] += val; }
        }
}
/* g[S*n+v] = number of paths starting at s, covering exactly S, ending at v */
static void start_table(const Tour *T, int s, u64 *g) {
    int n = T->n; size_t N = (size_t)1 << n;
    memset(g, 0, sizeof(u64) * N * n);
    g[((size_t)1 << s) * n + s] = 1;
    for (size_t S = 1; S < N; S++) {
        if (!((S >> s) & 1)) continue;
        for (int v = 0; v < n; v++) {
            if (!((S >> v) & 1)) continue;
            u64 val = g[S * n + v]; if (!val) continue;
            uint32_t cand = T->out[v] & ~(uint32_t)S;
            while (cand) { int w = __builtin_ctz(cand); cand &= cand - 1; g[(S | ((size_t)1 << w)) * n + w] += val; }
        }
    }
}

/* ---------- DFS for the cycle ---------- */
typedef struct { const Tour *T; int cnt, cap; int path[MAXV]; int cyc[MAXV]; } Dfs;
static int dfs_rec(Dfs *c, int cur, uint32_t vis, int depth) {
    const Tour *T = c->T;
    if (depth == T->n) {
        if (T->out[cur] & 1u) { if (c->cnt == 0) memcpy(c->cyc, c->path, sizeof(int) * T->n); c->cnt++; }
        return c->cnt >= c->cap;
    }
    uint32_t cand = T->out[cur] & ~vis;
    while (cand) {
        int w = __builtin_ctz(cand); cand &= cand - 1;
        c->path[depth] = w;
        if (dfs_rec(c, w, vis | (1u << w), depth + 1)) return 1;
    }
    return 0;
}
static int hc_capped(const Tour *T, int cap, int *cyc) {
    Dfs c; c.T = T; c.cnt = 0; c.cap = cap; c.path[0] = 0;
    dfs_rec(&c, 0, 1u, 1);
    if (c.cnt > 0 && cyc) memcpy(cyc, c.cyc, sizeof(int) * T->n);
    return c.cnt;
}

/* ---------- canonical form ---------- */
static u128 rot_word(const Tour *T, const int *cyc, int s) {
    int n = T->n; u128 c = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++)
        c = (c << 1) | (u128)arc(T, cyc[(s + i) % n], cyc[(s + j) % n]);
    return c;
}
static u128 canon_cycle(const Tour *T, const int *cyc, int *asym) {
    int n = T->n; u128 w[MAXV], best = 0;
    for (int s = 0; s < n; s++) { w[s] = rot_word(T, cyc, s); if (s == 0 || w[s] < best) best = w[s]; }
    if (asym) { *asym = 1; for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (w[a] == w[b]) *asym = 0; }
    return best;
}
static void decode(u128 code, int n, Tour *T) {
    T->n = n; memset(T->out, 0, sizeof T->out);
    int b = n * (n - 1) / 2 - 1;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        int bit = (int)((code >> b) & 1); b--;
        if (bit) T->out[i] |= 1u << j; else T->out[j] |= 1u << i;
    }
    finish(T);
}
/* brute-force canonical form over all permutations (small n only) */
static u128 canon_brute(const Tour *T) {
    int n = T->n, perm[MAXV], c[MAXV];
    for (int i = 0; i < n; i++) { perm[i] = i; c[i] = 0; }
    u128 best = 0; int first = 1;
    /* Heap's algorithm */
    int i = 0;
    for (;;) {
        u128 w = 0;
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) w = (w << 1) | (u128)arc(T, perm[a], perm[b]);
        if (first || w < best) { best = w; first = 0; }
        while (i < n && c[i] >= i) { c[i] = 0; i++; }
        if (i >= n) break;
        if (i % 2 == 0) { int t = perm[0]; perm[0] = perm[i]; perm[i] = t; }
        else { int t = perm[c[i]]; perm[c[i]] = perm[i]; perm[i] = t; }
        c[i]++; i = 0;
    }
    return best;
}

/* ---------- hash set of u128 ---------- */
typedef struct { u128 *k; int *val; size_t cap, cnt; } HSet;
static u64 hmix(u128 k) {
    u64 a = (u64)k, b = (u64)(k >> 64);
    u64 h = a * 0x9E3779B97F4A7C15ULL ^ (b + 0x632BE59BD9B4E019ULL) * 0xC2B2AE3D27D4EB4FULL;
    h ^= h >> 31; h *= 0x94D049BB133111EBULL; h ^= h >> 29; return h;
}
static void hs_init(HSet *h, size_t expect) {
    size_t cap = 64; while (cap < expect * 2 + 64) cap <<= 1;
    h->cap = cap; h->cnt = 0; h->k = calloc(cap, sizeof(u128)); h->val = calloc(cap, sizeof(int));
}
static void hs_free(HSet *h) { free(h->k); free(h->val); }
static long hs_find(const HSet *h, u128 key) {
    size_t m = h->cap - 1, i = hmix(key) & m;
    while (h->k[i]) { if (h->k[i] == key) return (long)i; i = (i + 1) & m; }
    return -1;
}
static int hs_insert(HSet *h, u128 key) {
    size_t m = h->cap - 1, i = hmix(key) & m;
    while (h->k[i]) { if (h->k[i] == key) return 0; i = (i + 1) & m; }
    h->k[i] = key; h->cnt++; return 1;
}
static int cmp128(const void *a, const void *b) {
    u128 x = *(const u128 *)a, y = *(const u128 *)b; return x < y ? -1 : x > y ? 1 : 0;
}
static u128 *hs_sorted(const HSet *h) {
    u128 *arr = malloc(sizeof(u128) * (h->cnt + 1)); size_t c = 0;
    for (size_t i = 0; i < h->cap; i++) if (h->k[i]) arr[c++] = h->k[i];
    qsort(arr, c, sizeof(u128), cmp128);
    return arr;
}

/* ---------- structure ---------- */
static int strong_sub(const Tour *T, uint32_t M) {
    if (!M) return 1;
    int s = __builtin_ctz(M);
    uint32_t seen = 1u << s, fr = seen;
    while (fr) { int u = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = T->out[u] & M & ~seen; seen |= nb; fr |= nb; }
    if (seen != M) return 0;
    seen = 1u << s; fr = seen;
    while (fr) { int u = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = T->in[u] & M & ~seen; seen |= nb; fr |= nb; }
    return seen == M;
}
static uint32_t nc_mask(const Tour *T) {
    uint32_t full = (1u << T->n) - 1, m = 0;
    for (int v = 0; v < T->n; v++) if (strong_sub(T, full & ~(1u << v))) m |= 1u << v;
    return m;
}
static void child(const Tour *T, int p, Tour *C) {
    int n = T->n; C->n = n + 1;
    for (int v = 0; v < n; v++) C->out[v] = T->out[v];
    C->out[p] |= 1u << n;
    C->out[n] = ((1u << n) - 1) & ~(1u << p);
    finish(C);
}
typedef struct { int y, x, r; uint32_t NC; int succ[MAXV], pred[MAXV]; } Run;
static int analyse_run(const Tour *T, const int *cyc, Run *R) {
    int n = T->n;
    for (int i = 0; i < n; i++) { R->succ[cyc[i]] = cyc[(i + 1) % n]; R->pred[cyc[i]] = cyc[(i + n - 1) % n]; }
    R->NC = nc_mask(T); R->r = __builtin_popcount(R->NC);
    uint32_t full = (1u << n) - 1;
    if (R->NC == 0 || R->NC == full) return 0;
    int starts = 0, y = -1;
    for (int v = 0; v < n; v++) if (((R->NC >> v) & 1) && !((R->NC >> R->pred[v]) & 1)) { starts++; y = v; }
    if (starts != 1) return 0;
    int x = y; while ((R->NC >> R->succ[x]) & 1) x = R->succ[x];
    R->y = y; R->x = x;
    if (__builtin_popcount(T->in[y]) != 1 || __builtin_popcount(T->out[x]) != 1) return 0;
    for (int v = 0; v < n; v++)
        if (v != y && ((R->NC >> v) & 1) && __builtin_popcount(T->in[v]) == 1) return 0;
    return 1;
}
/* remove vertex z, relabel the rest in increasing order; map[old] = new or -1 */
static void delete_vertex(const Tour *T, int z, Tour *U, int *map) {
    int n = T->n, k = 0;
    for (int v = 0; v < n; v++) map[v] = (v == z) ? -1 : k++;
    U->n = n - 1; memset(U->out, 0, sizeof U->out);
    for (int u = 0; u < n; u++) if (u != z)
        for (int v = 0; v < n; v++) if (v != z && arc(T, u, v)) U->out[map[u]] |= 1u << map[v];
    finish(U);
}

static long long fib(int k) { long long a = 0, b = 1; for (int i = 0; i < k; i++) { long long t = a + b; a = b; b = t; } return a; }

/* ---------- generation ---------- */
static u128 *cls[MAXV + 2]; static size_t ncls[MAXV + 2];
static u64 *hcA[MAXV + 2];   /* hcA[m][i*m + p] = exact hc(U_i^p), U_i the i-th class of order m */

static void gen_methodA(int m, HSet *out) {
    size_t cnt = ncls[m];
    u128 *cand = calloc(cnt * m, sizeof(u128));
    hcA[m] = calloc(cnt * m, sizeof(u64));
    size_t fsz = ((size_t)1 << m) * (m + 1);
#pragma omp parallel
    {
        u64 *f = malloc(sizeof(u64) * fsz);
#pragma omp for schedule(dynamic, 1)
        for (long i = 0; i < (long)cnt; i++) {
            Tour U; decode(cls[m][i], m, &U);
            for (int p = 0; p < m; p++) {
                Tour C; child(&U, p, &C);
                u64 h = hc_exact(&C, f);
                hcA[m][(size_t)i * m + p] = h;
                if (h == 1) {
                    int cyc[MAXV]; int hc2 = hc_capped(&C, 2, cyc);
                    if (hc2 != 1) { fprintf(stderr, "DFS/DP disagree (A)\n"); exit(2); }
                    cand[(size_t)i * m + p] = canon_cycle(&C, cyc, NULL);
                }
            }
        }
        free(f);
    }
    for (size_t i = 0; i < cnt * m; i++) if (cand[i]) hs_insert(out, cand[i]);
    free(cand);
}

static void gen_methodB(int m, HSet *out, int verify, long long *hits_out, long long *verify_bad) {
    size_t cnt = ncls[m];
    int nth = omp_get_max_threads();
    u128 **buf = calloc(nth, sizeof(u128 *)); size_t *bn = calloc(nth, sizeof(size_t)), *bc = calloc(nth, sizeof(size_t));
    long long vb = 0;
    size_t gsz = ((size_t)1 << m) * m, fsz = ((size_t)1 << m) * (m + 1);
#pragma omp parallel reduction(+:vb)
    {
        int t = omp_get_thread_num();
        u64 *g = malloc(sizeof(u64) * gsz), *f = malloc(sizeof(u64) * fsz);
        u64 hp[MAXV][MAXV];
        bc[t] = 1024; buf[t] = malloc(sizeof(u128) * bc[t]); bn[t] = 0;
#pragma omp for schedule(dynamic, 1)
        for (long i = 0; i < (long)cnt; i++) {
            Tour U; decode(cls[m][i], m, &U);
            size_t full = ((size_t)1 << m) - 1;
            for (int s = 0; s < m; s++) { start_table(&U, s, g); for (int v = 0; v < m; v++) hp[s][v] = g[full * m + v]; }
            for (uint32_t O = 1; O < (uint32_t)full; O++) {
                u64 h = 0;
                for (int w = 0; w < m && h < 2; w++) if ((O >> w) & 1)
                    for (int u = 0; u < m; u++) if (!((O >> u) & 1)) h += hp[w][u];
                if (h != 1) continue;
                Tour T; T.n = m + 1; memcpy(T.out, U.out, sizeof U.out);
                T.out[m] = O;
                for (int u = 0; u < m; u++) if (!((O >> u) & 1)) T.out[u] |= 1u << m;
                finish(&T);
                if (verify && hc_exact(&T, f) != 1) vb++;
                int cyc[MAXV]; if (hc_capped(&T, 2, cyc) != 1) vb++;
                if (bn[t] == bc[t]) { bc[t] *= 2; buf[t] = realloc(buf[t], sizeof(u128) * bc[t]); }
                buf[t][bn[t]++] = canon_cycle(&T, cyc, NULL);
            }
        }
        free(g); free(f);
    }
    long long hits = 0;
    for (int t = 0; t < nth; t++) { for (size_t j = 0; j < bn[t]; j++) hs_insert(out, buf[t][j]); hits += bn[t]; free(buf[t]); }
    free(buf); free(bn); free(bc);
    *hits_out = hits; *verify_bad = vb;
}

static void exhaustive(int n, HSet *out) {
    int pi[64], pj[64], M = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { pi[M] = i; pj[M] = j; M++; }
    long long total = 1LL << M;
    int nth = omp_get_max_threads();
    u128 **buf = calloc(nth, sizeof(u128 *)); size_t *bn = calloc(nth, sizeof(size_t)), *bc = calloc(nth, sizeof(size_t));
#pragma omp parallel
    {
        int t = omp_get_thread_num();
        bc[t] = 1024; buf[t] = malloc(sizeof(u128) * bc[t]); bn[t] = 0;
        u64 *f = malloc(sizeof(u64) * ((size_t)1 << (n - 1)) * n);
#pragma omp for schedule(dynamic, 4096)
        for (long long code = 0; code < total; code++) {
            Tour T; T.n = n; memset(T.out, 0, sizeof T.out);
            for (int b = 0; b < M; b++) { if ((code >> b) & 1) T.out[pi[b]] |= 1u << pj[b]; else T.out[pj[b]] |= 1u << pi[b]; }
            finish(&T);
            if (hc_exact(&T, f) != 1) continue;
            int cyc[MAXV]; if (hc_capped(&T, 2, cyc) != 1) { fprintf(stderr, "DFS/DP disagree (E)\n"); exit(2); }
            if (bn[t] == bc[t]) { bc[t] *= 2; buf[t] = realloc(buf[t], sizeof(u128) * bc[t]); }
            buf[t][bn[t]++] = canon_cycle(&T, cyc, NULL);
        }
        free(f);
    }
    for (int t = 0; t < nth; t++) { for (size_t j = 0; j < bn[t]; j++) hs_insert(out, buf[t][j]); free(buf[t]); }
    free(buf); free(bn); free(bc);
}

static int sets_equal(const HSet *a, const HSet *b) {
    if (a->cnt != b->cnt) return 0;
    for (size_t i = 0; i < a->cap; i++) if (a->k[i] && hs_find(b, a->k[i]) < 0) return 0;
    return 1;
}

/* ---------- checks ---------- */
typedef struct {
    long long bad_hc, bad_asym, bad_cin, bad_G, bad_L3, bad_H, bad_facts, bad_canon, bad_decomp, bad_cut, bad_pathsplit;
    long long adm2, adm3, admo, delta[3];
    long long Nr[MAXV + 2], N2r[MAXV + 2], Tk[MAXV + 2];
    long long CNr[MAXV + 3], CN2r[MAXV + 3], CTk[MAXV + 3], bad_childcut;   /* rows of the children, measured directly */
    u64 minW, minWminusHp, min_ep_U_at_pminus, min_hp_at_pminus;
} Chk;

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s NMAX NCHECK [LAB=7] [VERIFYB=12]\n", argv[0]); return 1; }
    int NMAX = atoi(argv[1]), NCHECK = atoi(argv[2]);
    int LAB = argc > 3 ? atoi(argv[3]) : 7, VERIFYB = argc > 4 ? atoi(argv[4]) : 12;
    if (NCHECK > NMAX) NCHECK = NMAX;   /* NCHECK == NMAX: "top-order mode" for the last order, see below */
    double t0 = omp_get_wtime();
    printf("# classes_audit2 NMAX=%d NCHECK=%d LAB=%d VERIFYB=%d threads=%d\n", NMAX, NCHECK, LAB, VERIFYB, omp_get_max_threads());

    /* order 4 exhaustively */
    HSet E4; hs_init(&E4, 16); exhaustive(4, &E4);
    cls[4] = hs_sorted(&E4); ncls[4] = E4.cnt;
    printf("order  4: exhaustive labelled classes = %zu (F(2)=%lld)\n", ncls[4], fib(2));
    hs_free(&E4);

    for (int m = 4; m < NMAX; m++) {
        int n1 = m + 1;
        double ta = omp_get_wtime();
        HSet A, B; hs_init(&A, (size_t)fib(2 * n1 - 6) + 16); hs_init(&B, (size_t)fib(2 * n1 - 6) + 16);
        gen_methodA(m, &A);
        double tb = omp_get_wtime();
        long long hitsB = 0, vbad = 0;
        gen_methodB(m, &B, n1 <= VERIFYB, &hitsB, &vbad);
        double tc = omp_get_wtime();
        int eqAB = sets_equal(&A, &B);
        printf("order %2d: methodA=%zu methodB=%zu (B hits %lld, B hc-verify bad %lld)  A==B:%s  F(2n-6)=%lld %s  [A %.1fs, B %.1fs]\n",
               n1, A.cnt, B.cnt, hitsB, vbad, eqAB ? "yes" : "NO", fib(2 * n1 - 6),
               (A.cnt == (size_t)fib(2 * n1 - 6)) ? "OK" : "MISMATCH", tb - ta, tc - tb);
        if (n1 <= LAB) {
            HSet E; hs_init(&E, (size_t)fib(2 * n1 - 6) + 16); exhaustive(n1, &E);
            printf("          exhaustive labelled classes = %zu, equal to method A set: %s\n", E.cnt, sets_equal(&A, &E) ? "yes" : "NO");
            hs_free(&E);
        }
        cls[n1] = hs_sorted(&A); ncls[n1] = A.cnt;
        hs_free(&A); hs_free(&B);
        fflush(stdout);
    }
    /* brute-force canonical forms as a sanity check of the rotation invariant */
    for (int n = 4; n <= 9 && n <= NMAX; n++) {
        HSet Bf; hs_init(&Bf, ncls[n] + 16);
        u128 *bfc = malloc(sizeof(u128) * ncls[n]);
#pragma omp parallel for schedule(dynamic, 1)
        for (long i = 0; i < (long)ncls[n]; i++) { Tour T; decode(cls[n][i], n, &T); bfc[i] = canon_brute(&T); }
        for (size_t i = 0; i < ncls[n]; i++) hs_insert(&Bf, bfc[i]);
        printf("order %2d: brute-force (all n! relabellings) distinct canonical forms = %zu of %zu classes\n", n, Bf.cnt, ncls[n]);
        free(bfc); hs_free(&Bf);
    }

    /* checks */
    printf("\n# checks of (G), (H), Lemma 3, decomposition, bijection, cut arcs\n");
    for (int n = 4; n <= NCHECK; n++) {
        Chk g; memset(&g, 0, sizeof g);
        g.minW = ~0ULL; g.minWminusHp = ~0ULL; g.min_ep_U_at_pminus = ~0ULL; g.min_hp_at_pminus = ~0ULL;
        size_t cnt = ncls[n];
        /* top-order mode (n == NMAX): order n+1 is not generated here.  Lemma 3 is then checked as
         * min(hc(T^p), 2) = min(e_p, 2) by a capped DFS, and the bijection is certified by collecting the
         * children's canonical words: they must be pairwise distinct and number |UH(n+1)| (established by an
         * independent generation, e.g. classes_audit2 16 ...).  */
        int top = (n == NMAX);
        int nth_top = omp_get_max_threads();
        u128 **tb = calloc(nth_top, sizeof(u128 *)); size_t *tbn = calloc(nth_top, sizeof(size_t)), *tbc = calloc(nth_top, sizeof(size_t));
        int *hitcount = calloc(ncls[n + 1] + 1, sizeof(int));
        HSet idx; hs_init(&idx, ncls[n + 1] + 16);
        for (size_t j = 0; j < ncls[n + 1]; j++) { hs_insert(&idx, cls[n + 1][j]); }
        /* store index in val */
        for (size_t j = 0; j < ncls[n + 1]; j++) { long pos = hs_find(&idx, cls[n + 1][j]); idx.val[pos] = (int)j; }
        size_t fsz = ((size_t)1 << (n + 1)) * (n + 1);
#pragma omp parallel
        {
            Chk l; memset(&l, 0, sizeof l);
            l.minW = ~0ULL; l.minWminusHp = ~0ULL; l.min_ep_U_at_pminus = ~0ULL; l.min_hp_at_pminus = ~0ULL;
            u64 *f = malloc(sizeof(u64) * fsz), *fu = malloc(sizeof(u64) * fsz), *gu = malloc(sizeof(u64) * fsz);
            int tid = omp_get_thread_num();
            if (top) { tbc[tid] = 4096; tb[tid] = malloc(sizeof(u128) * tbc[tid]); tbn[tid] = 0; }
#pragma omp for schedule(dynamic, 1)
            for (long i = 0; i < (long)cnt; i++) {
                Tour T; decode(cls[n][i], n, &T);
                int cyc[MAXV]; for (int k = 0; k < n; k++) cyc[k] = k;   /* representative is labelled along its cycle */
                { int cyc_ok = 1; for (int k = 0; k < n; k++) if (!arc(&T, k, (k + 1) % n)) cyc_ok = 0;
                  if (!cyc_ok) { l.bad_hc++; continue; } }
                if (hc_exact(&T, f) != 1) { l.bad_hc++; continue; }
                int asym; canon_cycle(&T, cyc, &asym); if (!asym) l.bad_asym++;
                Run R; if (!analyse_run(&T, cyc, &R)) { l.bad_cin++; continue; }
                int y = R.y, ym = R.pred[y], ymm = R.pred[ym];
                int ym_sole = __builtin_popcount(T.in[ym]) == 1;
                int nsole = 0; for (int v = 0; v < n; v++) if (__builtin_popcount(T.in[v]) == 1) nsole++;
                l.Nr[R.r]++; if (nsole == 2) l.N2r[R.r]++;
                /* e_p(T) */
                end_table(&T, f);
                size_t FULL = ((size_t)1 << n) - 1;
                u64 e[MAXV]; for (int p = 0; p < n; p++) e[p] = f[FULL * n + p];
                uint32_t adm = 0; for (int p = 0; p < n; p++) if (e[p] == 1) adm |= 1u << p;
                uint32_t pred_adm = (1u << y) | (1u << ym) | (ym_sole ? (1u << ymm) : 0u);
                if (adm != pred_adm) l.bad_G++;
                int a = __builtin_popcount(adm); if (a == 2) l.adm2++; else if (a == 3) l.adm3++; else l.admo++;
                /* Lemma 3 from method A's exact counts */
                if (!top) { for (int p = 0; p < n; p++) if (hcA[n][(size_t)i * n + p] != e[p]) l.bad_L3++; }
                else {      /* capped form; admissible p are checked (hc = 1) in the children loop */
                    for (int p = 0; p < n; p++) if (!((adm >> p) & 1)) {
                        Tour C; child(&T, p, &C);
                        if (hc_capped(&C, 2, NULL) != (e[p] >= 2 ? 2 : (int)e[p])) l.bad_L3++;
                    }
                }
                /* children */
                for (int p = 0; p < n; p++) {
                    if (!((adm >> p) & 1)) continue;
                    Tour C; child(&T, p, &C);
                    int ccyc[MAXV]; if (hc_capped(&C, 2, ccyc) != 1) { l.bad_H++; continue; }
                    Run RC; if (!analyse_run(&C, ccyc, &RC)) { l.bad_H++; continue; }
                    int expect = (p == y) ? 0 : (p == ym) ? 1 : 2;
                    if (RC.r - R.r != expect) l.bad_H++; else l.delta[expect]++;
                    int fok = 1;
                    if ((RC.NC >> p) & 1) fok = 0;
                    if (R.x == p) fok = 0;
                    if (__builtin_popcount(C.out[R.x]) != 1) fok = 0;
                    if ((RC.NC >> R.succ[R.x]) & 1) fok = 0;
                    if (expect == 2 && !((RC.NC >> ym) & 1)) fok = 0;
                    for (int v = 0; v < n; v++) if (v != p && ((R.NC >> v) & 1) && !((RC.NC >> v) & 1)) fok = 0;
                    if (!fok) l.bad_facts++;
                    /* the child's own run, sole entries and cut arcs, measured directly */
                    l.CNr[RC.r]++;
                    { int cs = 0; for (int v = 0; v <= n; v++) if (__builtin_popcount(C.in[v]) == 1) cs++; if (cs == 2) l.CN2r[RC.r]++; }
                    { int ck = 0;
                      for (int v = 0; v <= n; v++) {
                          Tour D = C; int w = RC.succ[v];
                          D.out[v] &= ~(1u << w); D.in[w] &= ~(1u << v);
                          if (!strong_sub(&D, (1u << (n + 1)) - 1)) ck++;
                      }
                      if (ck + RC.r != n + 2) l.bad_childcut++;
                      l.CTk[ck]++; }
                    /* canonicity: child's run start is the new vertex n with predecessor p; child - y = T */
                    if (RC.y != n || RC.pred[n] != p) l.bad_canon++;
                    u128 cc = canon_cycle(&C, ccyc, NULL);
                    if (top) {
                        if (tbn[tid] == tbc[tid]) { tbc[tid] *= 2; tb[tid] = realloc(tb[tid], sizeof(u128) * tbc[tid]); }
                        tb[tid][tbn[tid]++] = cc;
                    } else {
                        long pos = hs_find(&idx, cc);
                        if (pos < 0) l.bad_canon++;
                        else {
#pragma omp atomic
                            hitcount[idx.val[pos]]++;
                        }
                    }
                }
                /* decomposition identity (n >= 5): T = U^{p0}, U = T - y */
                if (n >= 5) {
                    Tour U; int map[MAXV]; delete_vertex(&T, y, &U, map);
                    int m = n - 1, p0 = map[ym];
                    /* check T = U^{p0}: y's in-neighbourhood is exactly {ym} */
                    if (T.in[y] != (1u << ym)) l.bad_decomp++;
                    end_table(&U, fu);
                    size_t FU = ((size_t)1 << m) - 1;
                    start_table(&U, p0, gu);
                    if (e[y] != fu[FU * m + p0] || e[y] != 1 || e[ym] != fu[FU * m + p0]) l.bad_decomp++;
                    for (int p = 0; p < n; p++) {
                        if (p == y || p == ym) continue;
                        int q = map[p];
                        u64 eU = fu[FU * m + q], hp = gu[FU * m + q];
                        u64 W = 0;
                        for (size_t S = 1; S < FU; S++) {
                            if (!((S >> p0) & 1) || ((S >> q) & 1)) continue;
                            W += fu[S * m + p0] * fu[(FU & ~S) * m + q];
                        }
                        if (e[p] != eU - hp + W) l.bad_decomp++;
                        if (W < l.minW) l.minW = W;
                        if (W < hp) l.bad_pathsplit++; else if (W - hp < l.minWminusHp) l.minWminusHp = W - hp;
                        if (p == ymm && !ym_sole) {   /* the case p = p0-, p0 not a sole entry */
                            if (eU < l.min_ep_U_at_pminus) l.min_ep_U_at_pminus = eU;
                            if (hp < l.min_hp_at_pminus) l.min_hp_at_pminus = hp;
                        }
                    }
                }
                /* cut arcs on C */
                int k = 0;
                for (int v = 0; v < n; v++) {
                    Tour D = T; int w = R.succ[v];
                    D.out[v] &= ~(1u << w); D.in[w] &= ~(1u << v);
                    if (!strong_sub(&D, (1u << n) - 1)) k++;
                }
                if (k + R.r != n + 1) l.bad_cut++;
                l.Tk[k]++;
            }
            free(f); free(fu); free(gu);
#pragma omp critical
            {
                g.bad_hc += l.bad_hc; g.bad_asym += l.bad_asym; g.bad_cin += l.bad_cin; g.bad_G += l.bad_G; g.bad_L3 += l.bad_L3;
                g.bad_H += l.bad_H; g.bad_facts += l.bad_facts; g.bad_canon += l.bad_canon; g.bad_decomp += l.bad_decomp;
                g.bad_cut += l.bad_cut; g.bad_pathsplit += l.bad_pathsplit;
                g.adm2 += l.adm2; g.adm3 += l.adm3; g.admo += l.admo;
                for (int q = 0; q < 3; q++) g.delta[q] += l.delta[q];
                for (int q = 0; q < MAXV + 2; q++) { g.Nr[q] += l.Nr[q]; g.N2r[q] += l.N2r[q]; g.Tk[q] += l.Tk[q]; }
                for (int q = 0; q < MAXV + 3; q++) { g.CNr[q] += l.CNr[q]; g.CN2r[q] += l.CN2r[q]; g.CTk[q] += l.CTk[q]; }
                g.bad_childcut += l.bad_childcut;
                if (l.minW < g.minW) g.minW = l.minW;
                if (l.minWminusHp < g.minWminusHp) g.minWminusHp = l.minWminusHp;
                if (l.min_ep_U_at_pminus < g.min_ep_U_at_pminus) g.min_ep_U_at_pminus = l.min_ep_U_at_pminus;
                if (l.min_hp_at_pminus < g.min_hp_at_pminus) g.min_hp_at_pminus = l.min_hp_at_pminus;
            }
        }
        long long hit0 = 0, hit1 = 0, hitmore = 0;
        for (size_t j = 0; j < ncls[n + 1]; j++) { if (hitcount[j] == 0) hit0++; else if (hitcount[j] == 1) hit1++; else hitmore++; }
        printf("n=%2d classes=%zu  bad: hc=%lld asym=%lld C-input=%lld G=%lld L3=%lld H=%lld facts=%lld canon=%lld decomp=%lld pathsplit=%lld cut=%lld\n",
               n, cnt, g.bad_hc, g.bad_asym, g.bad_cin, g.bad_G, g.bad_L3, g.bad_H, g.bad_facts, g.bad_canon, g.bad_decomp, g.bad_pathsplit, g.bad_cut);
        printf("      |adm|=2:%lld |adm|=3:%lld other:%lld   child deltas 0/1/2: %lld %lld %lld\n", g.adm2, g.adm3, g.admo, g.delta[0], g.delta[1], g.delta[2]);
        if (!top)
            printf("      bijection onto order %d (%zu classes): hit once %lld, never %lld, more than once %lld\n", n + 1, ncls[n + 1], hit1, hit0, hitmore);
        else {
            long long tot = 0; for (int t = 0; t < nth_top; t++) tot += (long long)tbn[t];
            HSet ch; hs_init(&ch, (size_t)tot + 16);
            for (int t = 0; t < nth_top; t++) for (size_t j = 0; j < tbn[t]; j++) hs_insert(&ch, tb[t][j]);
            printf("      top-order mode: children of order %d: %lld admissible pairs, %zu distinct children, %lld repeated; F(2(n+1)-6) = %lld\n",
                   n + 1, tot, ch.cnt, tot - (long long)ch.cnt, fib(2 * (n + 1) - 6));
            hs_free(&ch);
        }
        for (int t = 0; t < nth_top; t++) free(tb[t]);
        free(tb); free(tbn); free(tbc);
        if (n >= 5) {
            printf("      decomposition: min W = %llu, min (W - hp) = %llu", (unsigned long long)g.minW, (unsigned long long)g.minWminusHp);
            if (g.min_ep_U_at_pminus != ~0ULL)
                printf(", case p=p0-, p0 not sole entry: min e_p(U) = %llu, min hp_U(p0,p) = %llu", (unsigned long long)g.min_ep_U_at_pminus, (unsigned long long)g.min_hp_at_pminus);
            printf("\n");
        }
        printf("      N(n,r) r=2..n-1:"); for (int r = 2; r <= n - 1; r++) printf(" %lld", g.Nr[r]);
        printf("\n      N2(n,r) r=2..n-1:"); for (int r = 2; r <= n - 1; r++) printf(" %lld", g.N2r[r]);
        printf("\n      T(n,k) k=2..n-1:"); for (int k = 2; k <= n - 1; k++) printf(" %lld", g.Tk[k]);
        printf("\n      children (all of order %d, by the bijection above), measured directly; child cut+run != n+2: %lld\n", n + 1, g.bad_childcut);
        printf("      childN(%d,r) r=2..%d:", n + 1, n); for (int r = 2; r <= n; r++) printf(" %lld", g.CNr[r]);
        printf("\n      childN2(%d,r) r=2..%d:", n + 1, n); for (int r = 2; r <= n; r++) printf(" %lld", g.CN2r[r]);
        printf("\n      childT(%d,k) k=2..%d:", n + 1, n); for (int k = 2; k <= n; k++) printf(" %lld", g.CTk[k]);
        printf("\n");
        fflush(stdout);
        free(hitcount); hs_free(&idx);
    }
    printf("# total time %.1fs\n", omp_get_wtime() - t0);
    return 0;
}
