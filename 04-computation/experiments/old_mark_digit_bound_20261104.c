/*
 * old_mark_digit_bound_20261104.c  (opus-2026-11-04-S78)
 *
 * The old-mark digit lower bound for C2 children, checked on every tournament class B of order <= NMAX and every mark a.
 * T1 = C2(B, a): x -> a only (all other old vertices beat x); y beaten by all of B; y -> x.
 * Digits sigma_v = B_v - h, B_v = sum_{S subset V - v} h(T[S]) s_v(T[V - S]); d = out-degree.
 * Checks (failure counts are printed; all should be 0):
 *   (D)    sigma_v >= d(d - 1) at every vertex of B and of T1;
 *   native sigma_x(T1) = 0, sigma_y(T1) = h - e_a, H(T1) = h + 2 B_a;
 *   sink   for the sink z != a of B: sigma_z(T1) = 3 h(B - z) + 2 B_a(B - z) - e_a(B - z);
 *   old    every old w != a with d_B(w) >= 1 has sigma_w(T1) >= (d_B(w) + 2)(d_B(w) + 1);
 *   min    every old w != a has sigma_w(T1) >= 4 (|B| >= 2), and sigma_w(T1) = 4 only when B is transitive, w its sink
 *          and a its penultimate vertex;
 *   mark   sigma_a(T1) = 0 if d_B(a) = 0; = e_w(B - a) + B_w(B - a) if a -> w is its only out-arc;
 *          >= (d + 1) d if d = d_B(a) >= 2.
 * Usage: omb NMAX
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

/* digits of every vertex: sig[v] = B_v - h, and e[v] */
static void all_stats(int n, const uint32_t *out, long long *dp, long long *Hs, long long *st, long long *sig, long long *ev, long long *hout) {
    int full = (1 << n) - 1;
    for (size_t i = 0; i < (size_t)(full + 1) * n; i++) dp[i] = 0;
    for (int x = 0; x < n; x++) dp[(size_t)(1 << x) * n + x] = 1;
    for (int S = 1; S <= full; S++) for (int x = 0; x < n; x++) { long long c = dp[(size_t)S * n + x]; if (!c) continue;
        uint32_t nb = out[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << y)) * n + y] += c; } }
    Hs[0] = 1; for (int S = 1; S <= full; S++) { long long t = 0; for (int x = 0; x < n; x++) t += dp[(size_t)S * n + x]; Hs[S] = t; }
    for (int v = 0; v < n; v++) ev[v] = dp[(size_t)full * n + v];
    uint32_t inn[MAXN]; for (int v = 0; v < n; v++) { inn[v] = 0; for (int u = 0; u < n; u++) if ((out[u] >> v) & 1) inn[v] |= 1u << u; }
    for (size_t i = 0; i < (size_t)(full + 1) * n; i++) st[i] = 0;
    for (int x = 0; x < n; x++) st[(size_t)(1 << x) * n + x] = 1;
    for (int S = 1; S <= full; S++) for (int x = 0; x < n; x++) { long long c = st[(size_t)S * n + x]; if (!c) continue;
        uint32_t nb = inn[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; st[(size_t)(S | (1 << y)) * n + y] += c; } }
    *hout = Hs[full];
    for (int v = 0; v < n; v++) { long long Bv = 0; for (int S = 0; S <= full; S++) if (!((S >> v) & 1)) Bv += Hs[S] * st[(size_t)(full & ~S) * n + v]; sig[v] = Bv - Hs[full]; }
}
/* stats of T - r: h, and for vertex q: e_q, B_q */
static void sub_stats(int n, const uint32_t *out, int r, int q, long long *dp, long long *Hs, long long *st, long long *h, long long *eq, long long *Bq) {
    uint32_t o2[MAXN]; int idx[MAXN], m = 0, qi = -1;
    for (int v = 0; v < n; v++) if (v != r) { if (v == q) qi = m; idx[m++] = v; }
    for (int i = 0; i < m; i++) { o2[i] = 0; for (int j = 0; j < m; j++) if ((out[idx[i]] >> idx[j]) & 1) o2[i] |= 1u << j; }
    long long sg[MAXN], ev[MAXN], hh;
    if (m == 0) { *h = 1; *eq = 0; *Bq = 0; return; }
    all_stats(m, o2, dp, Hs, st, sg, ev, &hh);
    *h = hh; *eq = ev[qi]; *Bq = sg[qi] + hh;
}
int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 7;
    printf("# old_mark_digit_bound_20261104 NMAX=%d (orders of B)\n", NMAX);
    u128 *cur = malloc(sizeof(u128) * 2); size_t ncur;
    { uint32_t one[1] = {0}; cur[0] = canon(1, one); ncur = 1; }
    for (int n = 1; n <= NMAX; n++) {
        long pairs = 0, fD = 0, fnat = 0, fsink = 0, fdeg = 0, fmin = 0, feq = 0, fa0 = 0, fa1 = 0, fa2 = 0, eq4 = 0, sinkcases = 0, a1cases = 0;
        long long minold = -1, mina1 = -1, minnonsink = -1;
        #pragma omp parallel reduction(+:pairs,fD,fnat,fsink,fdeg,fmin,feq,fa0,fa1,fa2,eq4,sinkcases,a1cases)
        {
            int N2 = n + 2; size_t sz = ((size_t)1 << N2) * N2;
            long long *dp = malloc(sizeof(long long) * sz), *st = malloc(sizeof(long long) * sz), *Hs = malloc(sizeof(long long) << N2);
            long long lmin = -1, lmina1 = -1, lminns = -1;
            #pragma omp for schedule(dynamic, 8)
            for (size_t i = 0; i < ncur; i++) {
                uint32_t out[MAXN]; decode(n, cur[i], out);
                long long sB[MAXN], eB[MAXN], h;
                all_stats(n, out, dp, Hs, st, sB, eB, &h);
                int deg[MAXN], trans = 1, seen[MAXN];
                for (int v = 0; v < n; v++) seen[v] = 0;
                for (int v = 0; v < n; v++) { deg[v] = __builtin_popcount(out[v]); if (seen[deg[v]]) trans = 0; seen[deg[v]] = 1; }
                for (int v = 0; v < n; v++) if (sB[v] < (long long)deg[v] * (deg[v] - 1)) fD++;
                for (int a = 0; a < n; a++) {
                    pairs++;
                    uint32_t o2[MAXN]; for (int u = 0; u < n; u++) o2[u] = out[u];
                    int x = n, y = n + 1; o2[x] = 1u << a; o2[y] = 1u << x;
                    for (int u = 0; u < n; u++) { if (u != a) o2[u] |= 1u << x; o2[u] |= 1u << y; }
                    long long s1[MAXN], e1[MAXN], H1;
                    all_stats(N2, o2, dp, Hs, st, s1, e1, &H1);
                    if (s1[x] != 0 || s1[y] != h - eB[a] || H1 != h + 2 * (sB[a] + h)) fnat++;
                    for (int v = 0; v < N2; v++) { int dd = __builtin_popcount(o2[v]); if (s1[v] < (long long)dd * (dd - 1)) fD++; }
                    long long mo = -1;
                    for (int w = 0; w < n; w++) if (w != a) {
                        if (mo < 0 || s1[w] < mo) mo = s1[w];
                        if (deg[w] >= 1) { if (s1[w] < (long long)(deg[w] + 2) * (deg[w] + 1)) fdeg++; if (lminns < 0 || s1[w] < lminns) lminns = s1[w]; }
                        else {
                            sinkcases++;
                            long long hp, ep, Bp; sub_stats(n, out, w, a, dp, Hs, st, &hp, &ep, &Bp);
                            if (s1[w] != 3 * hp + 2 * Bp - ep) fsink++;
                        }
                        if (s1[w] == 4) { eq4++; int ok = trans && deg[w] == 0 && deg[a] == 1; if (!ok) feq++; }
                    }
                    if (n >= 2) { if (mo < 4) fmin++; if (lmin < 0 || mo < lmin) lmin = mo; }
                    if (deg[a] == 0) { if (s1[a] != 0) fa0++; }
                    else if (deg[a] == 1) {
                        a1cases++; int w = __builtin_ctz(out[a]);
                        long long hu, eu, Bu; sub_stats(n, out, a, w, dp, Hs, st, &hu, &eu, &Bu);
                        if (s1[a] != eu + Bu) fa1++;
                        if (lmina1 < 0 || s1[a] < lmina1) lmina1 = s1[a];
                    } else { if (s1[a] < (long long)(deg[a] + 1) * deg[a]) fa2++; }
                }
            }
            #pragma omp critical
            { if (lmin >= 0 && (minold < 0 || lmin < minold)) minold = lmin; if (lmina1 >= 0 && (mina1 < 0 || lmina1 < mina1)) mina1 = lmina1;
              if (lminns >= 0 && (minnonsink < 0 || lminns < minnonsink)) minnonsink = lminns; }
            free(dp); free(st); free(Hs);
        }
        printf("|B|=%d: classes %zu, (B,a) pairs %ld | failures: (D) %ld, native %ld, sink formula %ld (of %ld), non-sink old (d+2)(d+1) %ld, old min>=4 %ld, digit-4 equality case %ld (digit-4 old marks %ld), mark d=0 %ld, mark d=1 formula %ld (of %ld), mark d>=2 %ld | min old digit %lld, min non-sink old digit %lld, min mark digit at d=1 %lld\n",
               n, ncur, pairs, fD, fnat, fsink, sinkcases, fdeg, fmin, feq, eq4, fa0, fa1, a1cases, fa2, minold, minnonsink, mina1);
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
    return 0;
}
