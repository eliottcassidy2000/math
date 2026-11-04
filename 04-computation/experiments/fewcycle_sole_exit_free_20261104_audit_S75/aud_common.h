/*
 * aud_common.h  (independent audit of THM-4869, opus-2026-11-04-S75 audit)
 *
 * Shared routines written for the audit (no code taken from the session's programs):
 * - canonical form by individualization-refinement (equitable refinement at every node, all leaves, minimum code),
 *   with an optional initial vertex colouring (used for rooted tournaments and for (W, A) blocks);
 * - Hamiltonian cycle and path counts by subset DP (an endpoint recurrence from a fixed root);
 * - strong connectivity, 2-strong, arc-strong connectivity >= 2 (via single-arc deletions).
 */
#ifndef AUD_COMMON_H
#define AUD_COMMON_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXV 16
typedef unsigned long long u64;
typedef long long i64;
typedef unsigned __int128 u128;

/* ------------------------------------------------------------------ connectivity */
static int a_strong(int n, const uint32_t *out, uint32_t V) {
    (void)n;
    if (!V) return 1;
    int s = __builtin_ctz(V);
    uint32_t seen = 1u << s, fr = seen;
    while (fr) { int x = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = out[x] & V & ~seen; seen |= nb; fr |= nb; }
    if (seen != V) return 0;
    seen = 1u << s; fr = seen;
    while (fr) {
        int x = __builtin_ctz(fr); fr &= fr - 1;
        for (uint32_t r = V & ~seen; r; r &= r - 1) { int y = __builtin_ctz(r); if ((out[y] >> x) & 1u) { seen |= 1u << y; fr |= 1u << y; } }
    }
    return seen == V;
}
static int a_two_strong(int n, const uint32_t *out) {
    uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    if (!a_strong(n, out, full)) return 0;
    for (int v = 0; v < n; v++) if (!a_strong(n, out, full & ~(1u << v))) return 0;
    return 1;
}
/* arc-strong connectivity >= 2 (assumes T strong): no single arc whose deletion destroys strongness */
static int a_lambda2(int n, const uint32_t *out0) {
    uint32_t out[MAXV]; memcpy(out, out0, sizeof(uint32_t) * n);
    uint32_t full = (1u << n) - 1;
    for (int a = 0; a < n; a++) for (uint32_t r = out0[a]; r; r &= r - 1) {
        int b = __builtin_ctz(r);
        out[a] &= ~(1u << b);
        int st = a_strong(n, out, full);
        out[a] |= 1u << b;
        if (!st) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ path / cycle counts */
/* Hamiltonian cycles of T[V] (vertex mask V), via paths from the lowest vertex r of V.
 * buf must hold (1 << |V|) * |V| i64. */
static i64 a_hc_sub(const uint32_t *out, uint32_t V, i64 *buf) {
    int vs[MAXV], m = 0;
    for (uint32_t r = V; r; r &= r - 1) vs[m++] = __builtin_ctz(r);
    if (m < 3) return 0;
    uint32_t loc[MAXV];
    for (int i = 0; i < m; i++) { loc[i] = 0; for (int j = 0; j < m; j++) if ((out[vs[i]] >> vs[j]) & 1u) loc[i] |= 1u << j; }
    /* root = local 0; masks over local vertices 1..m-1 (bit j-1) */
    int k = m - 1; uint32_t full = (1u << k) - 1;
    memset(buf, 0, sizeof(i64) * (size_t)(full + 1) * k);
    for (int j = 1; j < m; j++) if ((loc[0] >> j) & 1u) buf[(size_t)(1u << (j - 1)) * k + (j - 1)] = 1;
    for (uint32_t S = 1; S <= full; S++) {
        for (uint32_t r = S; r; r &= r - 1) {
            int b = __builtin_ctz(r); i64 c = buf[(size_t)S * k + b]; if (!c) continue;
            uint32_t nb = (loc[b + 1] >> 1) & full & ~S;
            while (nb) { int b2 = __builtin_ctz(nb); nb &= nb - 1; buf[(size_t)(S | (1u << b2)) * k + b2] += c; }
        }
    }
    i64 h = 0;
    for (int j = 1; j < m; j++) if (loc[j] & 1u) h += buf[(size_t)full * k + (j - 1)];
    return h;
}

/* end counts e[x] (Hamiltonian paths of T[V] ending at x) and start counts s[x], for x in V (indexed by global vertex).
 * Uses a "paths ending at v over mask S" DP from all starts simultaneously for the end counts, and the same on the
 * converse for the start counts. buf: (1 << |V|) * |V| i64. */
static void a_end_counts(const uint32_t *out, uint32_t V, i64 *e, i64 *buf, int converse) {
    int vs[MAXV], m = 0;
    for (uint32_t r = V; r; r &= r - 1) vs[m++] = __builtin_ctz(r);
    uint32_t loc[MAXV];
    for (int i = 0; i < m; i++) {
        loc[i] = 0;
        for (int j = 0; j < m; j++) {
            int arc = converse ? ((out[vs[j]] >> vs[i]) & 1u) : ((out[vs[i]] >> vs[j]) & 1u);
            if (arc) loc[i] |= 1u << j;
        }
    }
    uint32_t full = (1u << m) - 1;
    memset(buf, 0, sizeof(i64) * (size_t)(full + 1) * m);
    for (int j = 0; j < m; j++) buf[(size_t)(1u << j) * m + j] = 1;
    for (uint32_t S = 1; S <= full; S++) {
        for (uint32_t r = S; r; r &= r - 1) {
            int b = __builtin_ctz(r); i64 c = buf[(size_t)S * m + b]; if (!c) continue;
            uint32_t nb = loc[b] & ~S;
            while (nb) { int b2 = __builtin_ctz(nb); nb &= nb - 1; buf[(size_t)(S | (1u << b2)) * m + b2] += c; }
        }
    }
    for (int j = 0; j < m; j++) e[vs[j]] = buf[(size_t)full * m + j];
}

/* all-pairs Hamiltonian path counts P[s*n+t] of T (order n), one DP per start. buf: (1<<n)*n i64. */
static void a_all_pairs(int n, const uint32_t *out, i64 *P, i64 *buf) {
    uint32_t full = (1u << n) - 1;
    for (int s = 0; s < n; s++) {
        memset(buf, 0, sizeof(i64) * (size_t)(full + 1) * n);
        buf[(size_t)(1u << s) * n + s] = 1;
        for (uint32_t S = 1; S <= full; S++) {
            if (!((S >> s) & 1u)) continue;
            for (uint32_t r = S; r; r &= r - 1) {
                int b = __builtin_ctz(r); i64 c = buf[(size_t)S * n + b]; if (!c) continue;
                uint32_t nb = out[b] & ~S;
                while (nb) { int b2 = __builtin_ctz(nb); nb &= nb - 1; buf[(size_t)(S | (1u << b2)) * n + b2] += c; }
            }
        }
        for (int t = 0; t < n; t++) P[s * n + t] = buf[(size_t)full * n + t];
    }
}

/* ------------------------------------------------------------------ canonical form (individualization-refinement) */
typedef struct {
    int n; const uint32_t *out;
    u128 best; int have;
    long leaves;
} ACan;

/* refine the ordered colouring col[] (colours 0..ncell-1, ordered) to an equitable one; returns the new ncell */
static int a_refine(int n, const uint32_t *out, int *col, int ncell) {
    for (;;) {
        int sig[MAXV][MAXV + 1];
        uint32_t cellmask[MAXV]; for (int c = 0; c < ncell; c++) cellmask[c] = 0;
        for (int v = 0; v < n; v++) cellmask[col[v]] |= 1u << v;
        for (int v = 0; v < n; v++) {
            sig[v][0] = col[v];
            for (int c = 0; c < ncell; c++) sig[v][c + 1] = __builtin_popcount(out[v] & cellmask[c]);
        }
        int ord[MAXV]; for (int i = 0; i < n; i++) ord[i] = i;
        for (int i = 1; i < n; i++) {
            int x = ord[i], j = i - 1;
            while (j >= 0 && memcmp(sig[ord[j]], sig[x], sizeof(int) * (ncell + 1)) > 0) { ord[j + 1] = ord[j]; j--; }
            ord[j + 1] = x;
        }
        int nc = 0, newcol[MAXV];
        for (int i = 0; i < n; i++) {
            if (i > 0 && memcmp(sig[ord[i]], sig[ord[i - 1]], sizeof(int) * (ncell + 1)) != 0) nc++;
            newcol[ord[i]] = nc;
        }
        nc++;
        memcpy(col, newcol, sizeof(int) * n);
        if (nc == ncell) return nc;
        ncell = nc;
    }
}

static void a_ir(ACan *C, int *col, int ncell) {
    int n = C->n;
    ncell = a_refine(n, C->out, col, ncell);
    if (ncell == n) {
        int perm[MAXV]; for (int v = 0; v < n; v++) perm[col[v]] = v;
        u128 code = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) code = (code << 1) | ((C->out[perm[i]] >> perm[j]) & 1u);
        C->leaves++;
        if (!C->have || code < C->best) { C->best = code; C->have = 1; }
        return;
    }
    /* target: first non-singleton cell */
    int cnt[MAXV] = {0};
    for (int v = 0; v < n; v++) cnt[col[v]]++;
    int tc = 0; while (cnt[tc] == 1) tc++;
    for (int v = 0; v < n; v++) if (col[v] == tc) {
        int c2[MAXV];
        for (int u = 0; u < n; u++) c2[u] = col[u] < tc ? col[u] : (u == v ? tc : col[u] + 1);
        a_ir(C, c2, ncell + 1);
    }
}

/* canonical code of T with an initial ordered colouring init[] (NULL = all equal) */
static u128 a_canon_col(int n, const uint32_t *out, const int *init) {
    ACan C; C.n = n; C.out = out; C.have = 0; C.best = 0; C.leaves = 0;
    int col[MAXV], ncell = 1;
    if (init) {
        /* compress initial colours to 0..k-1 preserving order */
        int mx = 0; for (int v = 0; v < n; v++) if (init[v] > mx) mx = init[v];
        int map[64]; for (int c = 0; c <= mx; c++) map[c] = -1;
        for (int v = 0; v < n; v++) map[init[v]] = 1;
        int k = 0; for (int c = 0; c <= mx; c++) if (map[c] == 1) map[c] = k++;
        for (int v = 0; v < n; v++) col[v] = map[init[v]];
        ncell = k;
    } else for (int v = 0; v < n; v++) col[v] = 0;
    a_ir(&C, col, ncell);
    return C.best;
}
static u128 a_canon(int n, const uint32_t *out) { return a_canon_col(n, out, NULL); }

/* decode: inverse of the leaf code (bits for pairs i<j in row-major order, MSB first) */
static void a_decode(int n, u128 code, uint32_t *out) {
    int tot = n * (n - 1) / 2, b = tot - 1;
    for (int i = 0; i < n; i++) out[i] = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        int bit = (int)((code >> b) & 1u); b--;
        if (bit) out[i] |= 1u << j; else out[j] |= 1u << i;
    }
}

/* ------------------------------------------------------------------ hash set of u128 */
typedef struct { u128 *key; char *used; size_t cap, cnt; } AHash;
static void ah_init(AHash *h, size_t cap) {
    h->cap = cap; h->cnt = 0; h->key = (u128 *)malloc(sizeof(u128) * cap); h->used = (char *)calloc(cap, 1);
    if (!h->key || !h->used) { fprintf(stderr, "hash alloc\n"); exit(1); }
}
static void ah_free(AHash *h) { free(h->key); free(h->used); }
static int ah_insert(AHash *h, u128 k) {
    if (h->cnt * 10 > h->cap * 7) { fprintf(stderr, "hash too full\n"); exit(1); }
    u64 x = (u64)k * 0xD6E8FEB86659FD93ULL ^ (u64)(k >> 64) * 0x9E3779B97F4A7C15ULL;
    x ^= x >> 29;
    size_t i = (size_t)(x % h->cap);
    while (h->used[i]) { if (h->key[i] == k) return 0; if (++i == h->cap) i = 0; }
    h->used[i] = 1; h->key[i] = k; h->cnt++; return 1;
}

static u64 a_fib(int k) { if (k < 0) return 0; u64 a = 0, b = 1; for (int i = 0; i < k; i++) { u64 t = a + b; a = b; b = t; } return a; }

#endif
