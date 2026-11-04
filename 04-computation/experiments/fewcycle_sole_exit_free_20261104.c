/*
 * fewcycle_sole_exit_free_20261104.c  (opus-2026-11-04-S75)
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
 * Usage: fewcycle NMAX K [DUMPN]   (NMAX <= 15, K <= 8; classes of order <= DUMPN with no sole exit are printed)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define MAXN 16
#define MAXK 8
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
typedef struct { u128 *key; signed char *hc; char *used; size_t cap, cnt; } HSet;
static void hs_init(HSet *h, size_t cap) {
    h->cap = cap; h->cnt = 0; h->key = malloc(sizeof(u128) * cap); h->hc = malloc(cap); h->used = calloc(cap, 1);
    if (!h->key || !h->hc || !h->used) { fprintf(stderr, "alloc\n"); exit(1); }
}
static int hs_insert(HSet *h, u128 k, int hc) {
    u64 x = (u64)k * 0x9E3779B97F4A7C15ULL ^ (u64)(k >> 64) * 0xC2B2AE3D27D4EB4FULL;
    size_t i = (size_t)(x >> 11) % h->cap;
    while (h->used[i]) { if (h->key[i] == k) return 0; i++; if (i == h->cap) i = 0; }
    h->used[i] = 1; h->key[i] = k; h->hc[i] = (signed char)hc; h->cnt++; return 1;
}

typedef struct { u128 code; int hc; } Rec;

static u64 fib(int k) { u64 a = 0, b = 1; if (k < 0) return 0; for (int i = 0; i < k; i++) { u64 t = a + b; a = b; b = t; } return a; }

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 10;
    int K = argc > 2 ? atoi(argv[2]) : 4;
    int DUMPN = argc > 3 ? atoi(argv[3]) : 0;
    if (NMAX > 15 || K > MAXK) { fprintf(stderr, "NMAX <= 15, K <= %d\n", MAXK); return 1; }
    printf("# fewcycle_sole_exit_free_20261104 NMAX=%d K=%d DUMPN=%d threads=%d\n", NMAX, K, DUMPN, omp_get_max_threads());

    /* level 3: the cyclic triangle */
    Rec *cur = malloc(sizeof(Rec)); size_t ncur = 1;
    { uint32_t o[3] = {2u, 4u, 1u}; cur[0].code = canon(3, o); cur[0].hc = 1; }

    for (int n = 3; n <= NMAX; n++) {
        /* ---------- statistics of the classes of order n ---------- */
        u64 t[MAXK + 1] = {0}, u[MAXK + 1] = {0}, uL1[MAXK + 1] = {0}, uL2cut[MAXK + 1] = {0}, u2s[MAXK + 1] = {0};
        u64 tL2cut[MAXK + 1] = {0}, t2s[MAXK + 1] = {0}, prodbad = 0, mincount_hist[8] = {0};
        u64 coretab[MAXK + 1][MAXN + 1][MAXK + 1]; memset(coretab, 0, sizeof coretab); /* [hc][|D|][e_u] for u classes, lambda = 1 */
        u64 entry_hist[MAXK + 1][4]; memset(entry_hist, 0, sizeof entry_hist);       /* sole entries among u classes */
        int dumped = 0; u64 slhist[MAXK + 1] = {0};
        #pragma omp parallel
        {
            i64 *buf = malloc(sizeof(i64) * ((size_t)1 << n) * n);
            i64 P[MAXN * MAXN];
            u64 lt[MAXK + 1] = {0}, lu[MAXK + 1] = {0}, luL1[MAXK + 1] = {0}, luL2cut[MAXK + 1] = {0}, lu2s[MAXK + 1] = {0};
            u64 ltL2cut[MAXK + 1] = {0}, lt2s[MAXK + 1] = {0}, lprodbad = 0, lmin[8] = {0};
            u64 lcore[MAXK + 1][MAXN + 1][MAXK + 1]; memset(lcore, 0, sizeof lcore);
            u64 lentry[MAXK + 1][4]; memset(lentry, 0, sizeof lentry);
            u64 lsl[MAXK + 1] = {0};
            #pragma omp for schedule(dynamic, 64)
            for (size_t i = 0; i < ncur; i++) {
                uint32_t out[MAXN]; decode(n, cur[i].code, out);
                int hc = cur[i].hc;
                uint32_t full = (1u << n) - 1;
                int soleexit = 0, soleentry = 0;
                for (int v = 0; v < n; v++) {
                    if (__builtin_popcount(out[v]) == 1) soleexit++;
                    if (n - 1 - __builtin_popcount(out[v]) == 1) soleentry++;
                }
                /* arc-strong connectivity: lambda = 1 iff some X has exactly one leaving arc */
                int lambda1 = 0;
                for (uint32_t X = 1; X < full && !lambda1; X++) {
                    int leave = 0;
                    for (uint32_t r = X; r && leave < 2; r &= r - 1) leave += __builtin_popcount(out[__builtin_ctz(r)] & ~X & full);
                    if (leave == 1) lambda1 = 1;
                }
                int cutv = 0;
                for (int v = 0; v < n && !cutv; v++) if (!is_strong_mask(n, out, full & ~(1u << v))) cutv = 1;
                lt[hc]++;
                if (!lambda1) { if (cutv) ltL2cut[hc]++; else lt2s[hc]++; }
                if (soleexit) continue;
                lu[hc]++;
                lentry[hc][soleentry < 3 ? soleentry : 3]++;
                if (!lambda1) { if (cutv) luL2cut[hc]++; else lu2s[hc]++; continue; }
                luL1[hc]++;
                /* minimal one-arc-leaving sets */
                uint32_t mins[64]; int nm = 0;
                for (uint32_t X = 1; X < full; X++) {
                    int leave = 0;
                    for (uint32_t r = X; r && leave < 2; r &= r - 1) leave += __builtin_popcount(out[__builtin_ctz(r)] & ~X & full);
                    if (leave != 1) continue;
                    /* minimal: no proper non-empty subset Y of X is one-arc-leaving */
                    int minimal = 1;
                    for (uint32_t Y = (X - 1) & X; Y && minimal; Y = (Y - 1) & X) {
                        int l2 = 0;
                        for (uint32_t r = Y; r && l2 < 2; r &= r - 1) l2 += __builtin_popcount(out[__builtin_ctz(r)] & ~Y & full);
                        if (l2 == 1) minimal = 0;
                    }
                    if (minimal && nm < 64) mins[nm++] = X;
                }
                lmin[nm < 7 ? nm : 7]++;
                for (int a = 0; a < nm; a++) {
                    uint32_t D = mins[a];
                    int vs[MAXN], m = 0, uu = -1, ll = -1;
                    for (int v = 0; v < n; v++) if ((D >> v) & 1) vs[m++] = v;
                    for (int k2 = 0; k2 < m; k2++) { uint32_t o = out[vs[k2]] & ~D & full; if (o) { uu = vs[k2]; ll = __builtin_ctz(o); } }
                    i64 eu = 0;
                    if (m == 1) eu = 1;
                    else {
                        all_pairs_sub(m, vs, out, P, buf);
                        int iu = 0; for (int k2 = 0; k2 < m; k2++) if (vs[k2] == uu) iu = k2;
                        for (int s = 0; s < m; s++) eu += P[s * m + iu];
                    }
                    /* prefix: s_l = Hamiltonian paths of T[V - D] starting at l (1 if |V - D| = 1) */
                    int qv[MAXN], q = 0;
                    for (int v = 0; v < n; v++) if (!((D >> v) & 1)) qv[q++] = v;
                    i64 sl = 0;
                    if (q == 1) sl = 1;
                    else {
                        all_pairs_sub(q, qv, out, P, buf);
                        int il = 0; for (int k2 = 0; k2 < q; k2++) if (qv[k2] == ll) il = k2;
                        for (int t2 = 0; t2 < q; t2++) sl += P[il * q + t2];
                    }
                    if (sl * eu != hc) lprodbad++;
                    if (a == 0) { int si = sl <= MAXK ? (int)sl : MAXK; lsl[si]++; }
                    int eidx = eu <= MAXK ? (int)eu : MAXK;
                    if (a == 0) lcore[hc][m][eidx]++;
                }
                if (DUMPN && n <= DUMPN) {
                    #pragma omp critical
                    {
                        printf("  DUMP n=%d hc=%d lambda=%d cutvertex=%d soleentries=%d minimal-one-arc-sets=%d sizes:", n, hc, lambda1 ? 1 : 2, cutv, soleentry, nm);
                        for (int a = 0; a < nm; a++) printf(" %d", __builtin_popcount(mins[a]));
                        printf("  rows:");
                        for (int v = 0; v < n; v++) { printf(" "); for (int w = 0; w < n; w++) printf("%d", (out[v] >> w) & 1); }
                        printf("\n");
                        dumped++;
                    }
                }
            }
            #pragma omp critical
            {
                for (int k = 0; k <= K; k++) {
                    t[k] += lt[k]; u[k] += lu[k]; uL1[k] += luL1[k]; uL2cut[k] += luL2cut[k]; u2s[k] += lu2s[k];
                    tL2cut[k] += ltL2cut[k]; t2s[k] += lt2s[k];
                    for (int a = 0; a <= MAXN; a++) for (int b = 0; b <= MAXK; b++) coretab[k][a][b] += lcore[k][a][b];
                    for (int b = 0; b < 4; b++) entry_hist[k][b] += lentry[k][b];
                }
                prodbad += lprodbad;
                for (int b = 0; b <= MAXK; b++) slhist[b] += lsl[b];
                for (int b = 0; b < 8; b++) mincount_hist[b] += lmin[b];
            }
            free(buf);
        }
        printf("n=%2d classes(hc<=%d)=%zu   F(2n-6)=%llu F(2n-9)=%llu 2F(2n-11)=%llu\n", n, K, ncur,
               fib(2 * n - 6), fib(2 * n - 9), 2 * fib(2 * n - 11));
        for (int k = 1; k <= K; k++)
            printf("   k=%d  t=%llu  [2-strong %llu, 2-arc-strong with cut vertex %llu]   u=%llu  [lambda=1 %llu, 2-arc-strong cut vertex %llu, 2-strong %llu]  u by #sole entries (0,1,2,3+): %llu %llu %llu %llu\n",
                   k, t[k], t2s[k], tL2cut[k], u[k], uL1[k], uL2cut[k], u2s[k],
                   entry_hist[k][0], entry_hist[k][1], entry_hist[k][2], entry_hist[k][3]);
        printf("   product check hc(T) = s_l(T - D) e_u(D) over minimal one-arc sets: failures %llu;  s_l of the first minimal set (0..%d):", prodbad, MAXK);
        for (int b = 0; b <= MAXK; b++) printf(" %llu", slhist[b]);
        printf(";  #minimal sets per u-class (0..6,7+):");
        for (int b = 0; b < 8; b++) printf(" %llu", mincount_hist[b]);
        printf("\n");
        for (int k = 1; k <= K; k++) for (int a = 1; a <= n; a++) for (int b = 0; b <= MAXK; b++)
            if (coretab[k][a][b]) printf("   core[k=%d] first minimal set |D|=%d e_u(D)=%d : %llu\n", k, a, b, coretab[k][a][b]);
        fflush(stdout);
        if (n == NMAX) break;

        /* ---------- extend to order n + 1 ---------- */
        int m = n;
        size_t cap = 1024; while (cap < ncur * 8 + 1024) cap <<= 1;
        HSet hs; hs_init(&hs, cap * 4);
        Rec *nxt = malloc(sizeof(Rec) * cap); size_t nnxt = 0;
        #pragma omp parallel
        {
            i64 *buf = malloc(sizeof(i64) * ((size_t)1 << m) * m);
            i64 P[MAXN * MAXN];
            Rec *loc = malloc(sizeof(Rec) * 4096); size_t nloc = 0, caploc = 4096;
            #pragma omp for schedule(dynamic, 16)
            for (size_t i = 0; i < ncur; i++) {
                uint32_t out[MAXN]; decode(m, cur[i].code, out);
                int vs[MAXN]; for (int k2 = 0; k2 < m; k2++) vs[k2] = k2;
                all_pairs_sub(m, vs, out, P, buf);
                uint32_t full = (1u << m) - 1, O = 0;
                i64 f = 0;
                for (uint32_t g = 1; g < ((uint32_t)1 << m); g++) {
                    int x = __builtin_ctz(g);
                    if ((O >> x) & 1) {
                        uint32_t O2 = O & ~(1u << x);
                        for (uint32_t r = O2; r; r &= r - 1) f += P[__builtin_ctz(r) * m + x];
                        for (uint32_t r = full & ~O; r; r &= r - 1) f -= P[x * m + __builtin_ctz(r)];
                        O = O2;
                    } else {
                        for (uint32_t r = O; r; r &= r - 1) f -= P[__builtin_ctz(r) * m + x];
                        for (uint32_t r = full & ~O & ~(1u << x); r; r &= r - 1) f += P[x * m + __builtin_ctz(r)];
                        O |= 1u << x;
                    }
                    if (O == 0 || O == full) continue;
                    if (f < 1 || f > K) continue;
                    uint32_t o2[MAXN];
                    for (int v = 0; v < m; v++) o2[v] = out[v] | (((O >> v) & 1) ? 0u : (1u << m));
                    o2[m] = O;
                    u128 c = canon(m + 1, o2);
                    if (nloc == caploc) { caploc *= 2; loc = realloc(loc, sizeof(Rec) * caploc); }
                    loc[nloc].code = c; loc[nloc].hc = (int)f; nloc++;
                }
            }
            #pragma omp critical
            {
                for (size_t j = 0; j < nloc; j++) {
                    if (hs_insert(&hs, loc[j].code, loc[j].hc)) {
                        if (nnxt == cap) { cap *= 2; nxt = realloc(nxt, sizeof(Rec) * cap); }
                        nxt[nnxt++] = loc[j];
                    }
                }
            }
            free(loc); free(buf);
        }
        free(hs.key); free(hs.hc); free(hs.used);
        free(cur); cur = nxt; ncur = nnxt;
    }
    return 0;
}
