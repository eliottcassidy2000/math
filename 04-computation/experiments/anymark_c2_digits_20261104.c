/*
 * anymark_c2_digits_20261104.c  (opus-2026-11-04-S77)
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
 * Any-mark C2 census (opus-2026-11-04-S77).  All tournament classes of orders 1..NMAX (strong or not), every vertex a:
 *   h = H(B), sigma_a = [P A_a](B) - h (the C2 child of (B, a) has H = 3h + 2 sigma_a, for every B and a).
 * Reports:
 *   - the any-mark readiness spectra R_s = {h : some (B, a) with sigma_a = s}, s = 0, 1, 2, against all odd h <= LIM;
 *   - the one-step any-mark sweep V = {3h + 2 sigma_a} and the odd numbers <= LIM it misses;
 *   - for orders <= NCHILD: the digits carried by every C2 child (all marks of the child), and whether digits 1 and 2 occur.
 */
static void tables(int n, const uint32_t *out, long long *dp, long long *Hs) {
    int full = (1 << n) - 1;
    for (size_t i = 0; i < (size_t)(full + 1) * n; i++) dp[i] = 0;
    for (int x = 0; x < n; x++) dp[(size_t)(1 << x) * n + x] = 1;
    for (int S = 1; S <= full; S++) for (int x = 0; x < n; x++) { long long c = dp[(size_t)S * n + x]; if (!c) continue;
        uint32_t nb = out[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << y)) * n + y] += c; } }
    Hs[0] = 1; for (int S = 1; S <= full; S++) { long long t = 0; for (int x = 0; x < n; x++) t += dp[(size_t)S * n + x]; Hs[S] = t; }
}
/* sigma for every mark: [P A_a] = sum_{S not containing a} h(S) * s_a(V - S); s_a via DP from a */
static void all_sigmas(int n, const uint32_t *out, long long *dp, long long *Hs, long long *sig, long long *hout) {
    tables(n, out, dp, Hs);
    int full = (1 << n) - 1; *hout = Hs[full];
    for (int a = 0; a < n; a++) {
        for (size_t i = 0; i < (size_t)(full + 1) * n; i++) dp[i] = 0;
        dp[(size_t)(1 << a) * n + a] = 1;
        long long pa = 0;
        for (int S = 1; S <= full; S++) { if (!((S >> a) & 1)) continue;
            for (int x = 0; x < n; x++) { long long c = dp[(size_t)S * n + x]; if (!c) continue;
                uint32_t nb = out[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << y)) * n + y] += c; } } }
        for (int S = 1; S <= full; S++) if ((S >> a) & 1) { long long st = 0; for (int x = 0; x < n; x++) st += dp[(size_t)S * n + x]; pa += Hs[full & ~S] * st; }
        sig[a] = pa - Hs[full];
    }
}

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 8, NCHILD = argc > 2 ? atoi(argv[2]) : 6;
    long long LIM = argc > 3 ? atoll(argv[3]) : 2000;
    printf("# anymark_c2_digits_20261104 NMAX=%d NCHILD=%d LIM=%lld\n", NMAX, NCHILD, LIM);
    char *R0 = calloc(LIM + 1, 1), *R1 = calloc(LIM + 1, 1), *R2 = calloc(LIM + 1, 1), *Hset = calloc(LIM + 1, 1);
    long long VL = 3 * LIM + 64; char *V = calloc(VL + 1, 1);
    u128 *cur = malloc(sizeof(u128) * 2); size_t ncur;
    { uint32_t one[1] = {0}; cur[0] = canon(1, one); ncur = 1; }
    long long *dp = malloc(sizeof(long long) * ((size_t)1 << (NMAX + 3)) * (NMAX + 3)), *Hs = malloc(sizeof(long long) << (NMAX + 3));
    long long sig[32], sig2[32];
    long child_cases = 0, child_has1 = 0, child_has2 = 0, child_has12 = 0;
    for (int n = 1; n <= NMAX; n++) {
        for (size_t i = 0; i < ncur; i++) {
            uint32_t out[MAXN]; decode(n, cur[i], out);
            long long h; all_sigmas(n, out, dp, Hs, sig, &h);
            if (h <= LIM) Hset[h] = 1;
            for (int a = 0; a < n; a++) {
                if (h <= LIM) { if (sig[a] == 0) R0[h] = 1; if (sig[a] == 1) R1[h] = 1; if (sig[a] == 2) R2[h] = 1; }
                long long v = 3 * h + 2 * sig[a]; if (v <= VL) V[v] = 1;
                if (n <= NCHILD) {
                    /* the C2 child of (B, a): x = n (x -> a only), y = n+1 (all old beat y, y -> x) */
                    uint32_t o2[MAXN]; for (int u = 0; u < n; u++) o2[u] = out[u];
                    int x = n, y = n + 1; o2[x] = 1u << a; o2[y] = 1u << x;
                    for (int u = 0; u < n; u++) { if (u != a) o2[u] |= 1u << x; o2[u] |= 1u << y; }
                    long long h2; all_sigmas(n + 2, o2, dp, Hs, sig2, &h2);
                    int has1 = 0, has2 = 0; for (int u = 0; u < n + 2; u++) { if (sig2[u] == 1) has1 = 1; if (sig2[u] == 2) has2 = 1; }
                    child_cases++; child_has1 += has1; child_has2 += has2; child_has12 += (has1 && has2);
                }
            }
        }
        long odd = 0, r1 = 0, r2 = 0, r12 = 0, hs = 0; long long lim_n = LIM;
        for (long long v = 1; v <= lim_n; v += 2) { if (!Hset[v]) continue; hs++; if (R1[v]) r1++; if (R2[v]) r2++; if (R1[v] && R2[v]) r12++; }
        printf("n<=%d: classes this order %zu | realised odd h<=%lld: %ld; with a digit-1 mark: %ld; digit-2: %ld; both: %ld\n",
               n, ncur, LIM, hs, r1, r2, r12);
        fflush(stdout);
        if (n == NMAX) break;
        HSet hsx; hs_init(&hsx, ncur * (1u << n) / 2 + 1024);
        size_t capn = ncur * (1u << n) / 4 + 1024, nn = 0; u128 *nxt = malloc(sizeof(u128) * capn);
        for (size_t i = 0; i < ncur; i++) {
            uint32_t out[MAXN]; decode(n, cur[i], out);
            for (uint32_t O = 0; O < (1u << n); O++) {
                uint32_t o2[MAXN]; for (int v = 0; v < n; v++) o2[v] = out[v] | (((O >> v) & 1) ? 0u : (1u << n)); o2[n] = O;
                u128 c = canon(n + 1, o2);
                if (hs_insert(&hsx, c, 0)) { if (nn == capn) { capn *= 2; nxt = realloc(nxt, sizeof(u128) * capn); } nxt[nn++] = c; }
            }
        }
        free(hsx.key); free(hsx.hc); free(hsx.used); free(cur); cur = nxt; ncur = nn;
    }
    printf("C2 children of all (B, a) with |B| <= %d: %ld cases; child has a digit-1 mark: %ld; digit-2 mark: %ld; both: %ld\n",
           NCHILD, child_cases, child_has1, child_has2, child_has12);
    printf("odd v <= %lld missed by the one-step any-mark sweep 3h + 2 sigma_a (orders <= %d):", LIM, NMAX);
    long miss = 0; for (long long v = 1; v <= LIM; v += 2) if (!V[v]) { miss++; if (miss <= 80) printf(" %lld", v); }
    printf("  [total %ld]\n", miss);
    printf("odd h <= 400 realised at orders <= %d but with no digit-1 mark:", NMAX);
    for (long long v = 1; v <= 400; v += 2) if (Hset[v] && !R1[v]) printf(" %lld", v);
    printf("\nodd h <= 400 realised at orders <= %d but with no digit-2 mark:", NMAX);
    for (long long v = 1; v <= 400; v += 2) if (Hset[v] && !R2[v]) printf(" %lld", v);
    printf("\n");
    return 0;
}
