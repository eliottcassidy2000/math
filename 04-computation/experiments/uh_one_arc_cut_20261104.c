/*
 * uh_one_arc_cut_20261104.c  (opus-2026-11-04-S74)
 *
 * Generates every uniquely Hamiltonian tournament (exactly one directed Hamiltonian cycle) up to isomorphism, order by
 * order, and tests the one-arc-cut property:
 *     T has a vertex set A (empty != A != V) with exactly one arc from V - A to A,
 * i.e. T is a one-arc join J(A, a; B, b) (THM-4866 (B)); together with vertex-degree statistics.
 *
 * Generation (complete, by the extension lemma of the note, section 2):
 *   every uniquely Hamiltonian T of order n >= 4 has a vertex v with T - v strong (Moon 1966); then T - v is uniquely
 *   Hamiltonian, and the out-neighbourhood O of v is a cyclic interval of the Hamiltonian cycle C' of T - v (both O and
 *   I = complement non-empty).  So U(n) = { T' + v : T' in U(n-1), O a proper non-empty cyclic interval of C',
 *   hc(T' + v) = 1 }, deduplicated by a canonical form.  hc(T' + v) = sum over w in O, u in I of P_T'(w -> u), where
 *   P_T'(w -> u) counts Hamiltonian paths of T' from w to u.  The class counts are checked against Garey's F(2n - 6).
 *
 * Usage: uh_one_arc_cut NMAX   (NMAX <= 16)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define MAXN 16
typedef unsigned long long u64;
typedef unsigned __int128 u128;

/* ---------------- canonical form (refined ordered partition + pruned search), codes in 128 bits ---------------- */
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

/* ---------------- Hamiltonian cycles ---------------- */
/* number of Hamiltonian cycles; if exactly one, its vertex order is written to cyc (cyc[i] -> cyc[i+1]) */
static u64 ham_cycles(int n, const uint32_t *out, int *cyc, u64 *buf) {
    int full = (1 << n) - 1;
    u64 *dp = buf;
    memset(dp, 0, sizeof(u64) * (size_t)(full + 1) * n);
    dp[(size_t)1 * n + 0] = 1;
    for (int S = 1; S <= full; S++) {
        if (!(S & 1)) continue;
        for (int v = 0; v < n; v++) {
            u64 c = dp[(size_t)S * n + v]; if (!c) continue;
            uint32_t nb = out[v] & ~(uint32_t)S;
            while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << w)) * n + w] += c; }
        }
    }
    u64 hc = 0; int last = -1;
    for (int v = 1; v < n; v++) if (out[v] & 1u) { hc += dp[(size_t)full * n + v]; if (dp[(size_t)full * n + v]) last = v; }
    if (hc == 1 && cyc) {
        /* backtrack the unique path 0 -> ... -> last */
        int S = full, v = last, k = n - 1;
        while (1) {
            cyc[k] = v;
            if (S == 1) break;
            int S2 = S & ~(1 << v), found = -1;
            for (int u = 0; u < n; u++) if (((S2 >> u) & 1) && ((out[u] >> v) & 1) && dp[(size_t)S2 * n + u]) { found = u; break; }
            S = S2; v = found; k--;
        }
    }
    return hc;
}

/* all-pairs Hamiltonian path counts P[s*m + t] of a tournament of order m */
static void all_pairs(int m, const uint32_t *out, u64 *P, u64 *buf) {
    int full = (1 << m) - 1;
    for (int s = 0; s < m; s++) {
        u64 *dp = buf;
        memset(dp, 0, sizeof(u64) * (size_t)(full + 1) * m);
        dp[(size_t)(1 << s) * m + s] = 1;
        for (int S = 1; S <= full; S++) {
            if (!((S >> s) & 1)) continue;
            for (int v = 0; v < m; v++) {
                u64 c = dp[(size_t)S * m + v]; if (!c) continue;
                uint32_t nb = out[v] & ~(uint32_t)S;
                while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << w)) * m + w] += c; }
            }
        }
        for (int t = 0; t < m; t++) P[s * m + t] = dp[(size_t)full * m + t];
    }
}

/* ---------------- hash set of 128-bit codes ---------------- */
typedef struct { u128 *key; char *used; size_t cap, cnt; } HSet;
static void hs_init(HSet *h, size_t cap) { h->cap = cap; h->cnt = 0; h->key = malloc(sizeof(u128) * cap); h->used = calloc(cap, 1); }
static int hs_insert(HSet *h, u128 k) {
    u64 x = (u64)k * 0x9E3779B97F4A7C15ULL ^ (u64)(k >> 64) * 0xC2B2AE3D27D4EB4FULL;
    size_t i = (size_t)(x >> 11) % h->cap;
    while (h->used[i]) { if (h->key[i] == k) return 0; i++; if (i == h->cap) i = 0; }
    h->used[i] = 1; h->key[i] = k; h->cnt++; return 1;
}

/* ---------------- tests ---------------- */
/* number of one-arc cuts: pairs (i, k) of cycle arcs; A = {c_(i+1), ..., c_k} (cyclic), W = rest; no arc W -> A except c_i -> c_(i+1) */
static int one_arc_cuts(int n, const uint32_t *out, const int *cyc, int *ncutarcs) {
    int cnt = 0; char isarc[MAXN] = {0};
    for (int i = 0; i < n; i++) {
        uint32_t A = 0;
        for (int len = 1; len < n; len++) {
            int k = (i + len) % n;            /* A = c_(i+1) .. c_(i+len) */
            A |= 1u << cyc[k];
            uint32_t W = ((1u << n) - 1) & ~A;
            int back = 0;
            for (int y = 0; y < n && back < 2; y++) if ((W >> y) & 1) back += __builtin_popcount(out[y] & A);
            if (back == 1) { cnt++; isarc[i] = 1; }
        }
    }
    int a = 0; for (int i = 0; i < n; i++) a += isarc[i];
    *ncutarcs = a;
    return cnt;
}

typedef struct { u128 code; int cyc[MAXN]; } Rec;

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 13;
    if (NMAX > MAXN) NMAX = MAXN;
    long long fib[64]; fib[1] = fib[2] = 1; for (int i = 3; i < 64; i++) fib[i] = fib[i - 1] + fib[i - 2];
    printf("uh_one_arc_cut NMAX=%d threads=%d\n", NMAX, omp_get_max_threads());
    /* order 3: the cyclic triangle 0 -> 1 -> 2 -> 0 */
    Rec *cur = malloc(sizeof(Rec)); size_t ncur = 1;
    { uint32_t o[MAXN] = {0}; o[0] = 1u << 1; o[1] = 1u << 2; o[2] = 1u << 0; cur[0].code = canon(3, o);
      uint32_t d[MAXN]; decode(3, cur[0].code, d); u64 buf[8 * 3]; ham_cycles(3, d, cur[0].cyc, buf); }
    for (int n = 3; n <= NMAX; n++) {
        double t0 = omp_get_wtime();
        if (n > 3) {
            int m = n - 1;
            size_t cap = 16; while (cap < (size_t)(fib[2 * n - 6] * 4 + 64)) cap <<= 1;
            HSet hs; hs_init(&hs, cap);
            size_t ncap = (size_t)fib[2 * n - 6] * 2 + 64;
            Rec *nxt = malloc(sizeof(Rec) * ncap); size_t nn = 0;
            long long cand = 0, ext1 = 0, ext2 = 0;
#pragma omp parallel
            {
                u64 *buf = malloc(sizeof(u64) * ((size_t)1 << n) * n);
                u64 *P = malloc(sizeof(u64) * m * m);
#pragma omp for schedule(dynamic, 1) reduction(+:cand, ext1, ext2)
                for (long r = 0; r < (long)ncur; r++) {
                    uint32_t tp[MAXN]; decode(m, cur[r].code, tp);
                    int cyc[MAXN]; memcpy(cyc, cur[r].cyc, sizeof(int) * m);
                    all_pairs(m, tp, P, buf);
                    for (int st = 0; st < m; st++) for (int L = 1; L < m; L++) {
                        uint32_t O = 0;
                        for (int k = 0; k < L; k++) O |= 1u << cyc[(st + k) % m];
                        uint32_t I = ((1u << m) - 1) & ~O;
                        u64 hc = 0;
                        for (int w = 0; w < m && hc < 2; w++) if ((O >> w) & 1) for (int u = 0; u < m; u++) if ((I >> u) & 1) hc += P[w * m + u];
                        cand++;
                        if (hc != 1) continue;
                        if (L == 1 || L == m - 1) ext1++; else ext2++;
                        uint32_t t[MAXN];
                        for (int i = 0; i < m; i++) t[i] = tp[i] | (((I >> i) & 1) ? (1u << m) : 0);
                        t[m] = O;
                        u128 code = canon(n, t);
                        int isnew;
#pragma omp critical
                        { isnew = hs_insert(&hs, code); if (isnew) { nxt[nn].code = code; nn++; } }
                    }
                }
                free(buf); free(P);
            }
            /* recompute the unique cycle in canonical labelling */
#pragma omp parallel
            {
                u64 *buf = malloc(sizeof(u64) * ((size_t)1 << n) * n);
#pragma omp for schedule(dynamic, 64)
                for (long r = 0; r < (long)nn; r++) {
                    uint32_t t[MAXN]; decode(n, nxt[r].code, t);
                    u64 hc = ham_cycles(n, t, nxt[r].cyc, buf);
                    if (hc != 1) { fprintf(stderr, "generated class with hc=%llu\n", hc); exit(1); }
                }
                free(buf);
            }
            free(cur); free(hs.key); free(hs.used);
            cur = nxt; ncur = nn;
            printf("n=%d: candidates %lld; single-cycle extensions with |O| = 1 or |I| = 1: %lld, with |O|, |I| >= 2: %lld\n", n, cand, ext1, ext2);
        }
        /* tests */
        long long withcut = 0, minin1 = 1 << 30, minout1 = 1 << 30, noin1 = 0, noout1 = 0, neither = 0, mincuts = 1 << 30, mincutarcs = 1 << 30;
        long long cutarc_hist[MAXN + 1] = {0};
        for (size_t r = 0; r < ncur; r++) {
            uint32_t t[MAXN]; decode(n, cur[r].code, t);
            int ca, c = one_arc_cuts(n, t, cur[r].cyc, &ca);
            if (c > 0) withcut++;
            if (c < mincuts) mincuts = c;
            if (ca < mincutarcs) mincutarcs = ca;
            cutarc_hist[ca]++;
            int in1 = 0, out1 = 0;
            for (int v = 0; v < n; v++) {
                int od = __builtin_popcount(t[v]);
                if (od == 1) out1++;
                if (n - 1 - od == 1) in1++;
            }
            if (in1 < minin1) minin1 = in1;
            if (out1 < minout1) minout1 = out1;
            if (!in1) noin1++;
            if (!out1) noout1++;
            if (!in1 && !out1) neither++;
        }
        printf("n=%d: %zu uniquely Hamiltonian classes (Garey F(2n-6) = %lld) | with a one-arc cut: %lld | min one-arc cuts per class %lld, min cut arcs on the cycle %lld | min #vertices of in-degree 1: %lld, of out-degree 1: %lld; classes with no in-degree-1 vertex %lld, with no out-degree-1 vertex %lld, with neither %lld | [%.1f s]\n",
               n, ncur, n >= 4 ? fib[2 * n - 6] : 1, withcut, mincuts, mincutarcs, minin1, minout1, noin1, noout1, neither, omp_get_wtime() - t0);
        printf("   cut arcs on the cycle (count: classes):");
        for (int k = 0; k <= n; k++) if (cutarc_hist[k]) printf(" %d:%lld", k, cutarc_hist[k]);
        printf("\n");
        fflush(stdout);
    }
    return 0;
}
