/*
 * anymark_c2_digits_20261104.c  (opus-2026-11-04-S77)
 *
 * Any-mark C2 census. Generates every tournament class of order <= NMAX and computes, for every mark a, the digit
 * sigma_a = [P A_a] - h. The C2 child of (B, a) has H = 3h + 2 sigma_a at every mark (RESULT law (4); THM-4865).
 * Reports
 *   - the any-mark readiness spectra R_s = {h : some (B, a) with sigma_a = s}, s = 0, 1, 2, among odd h <= LIM;
 *   - the one-step sweep {3h + 2 sigma_a} and the odd numbers <= LIM it misses;
 *   - for orders <= NCHILD, the digits carried by the marks of every C2 child.
 *
 * Usage: any NMAX NCHILD LIM   (any 8 6 1500 reproduces anymark_c2_digits_20261104.out)
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
