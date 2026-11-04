/*
 * scf_crevasses_20261104.c  (opus-2026-11-04-S73)
 *
 * The single-cycle floor and the crevasses of the Hamiltonian-path spectrum.
 *
 * Mode "census NMAX" (orders 3..NMAX <= 10, every isomorphism class, OpenMP):
 *   per class: H (Hamiltonian paths), hc (directed Hamiltonian cycles), self-converse flag; for even n, the Pfaffian of
 *   the +-1 skew-adjacency matrix (S_ij = +1 if i -> j, -1 if j -> i).
 *   Reports: the strong floor f(n) = min H over strong classes (strong <=> hc >= 1, Camion); the Hamiltonian cycle
 *   counts of the floor classes (single-cycle floor); the cycle ceiling H >= f(n) + 2(hc - 1) over every strong class
 *   (minimum slack, equality cases); the k-cycle floors f_k(n) = min{H : hc = k}; the self-converse ("real points")
 *   spectrum and its gaps; Pfaffian statistics against H and hc.
 * Mode "chains NMAX": every interval chain IC(p) (Moon's S_n; THM-4744) of order <= NMAX <= 20: H against
 *   prod (2^(p_i - 1) + 1) and the number of Hamiltonian cycles (claimed: exactly one).
 * Mode "joins KA KB VMAX": the one-arc-reversed join J = A => B with the arc a -> b reversed (b -> a).  With the rooted
 *   triple x(A, a) = (h(A), R_a(A), e_a(A)), R_a(A) = sum over X subset V(A)-a of p(X) * (paths of A - X starting at a)
 *   (klein's root response B_a, THM-4865), and the end-rooted triple y(B, b) = (h(B), R'_b(B), s_b(B)) = x(B^op, b):
 *       H(J) = x_h y_h + x_R y_R - x_e y_e      (a Lorentzian inner product of signature (2, 1)).
 *   Enumerates all distinct triples of orders <= max(KA, KB), counts representations r(N) of every odd N <= VMAX by
 *   pairs (x of order <= KA, y of order <= KB), verifies the formula on random joins by direct path counting, and
 *   reports coverage gaps, dyadic-window minima of r and the mean of r by residue classes.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define MAXN 10
typedef unsigned long long u64;

/* ---------------- reentrant canonical form ---------------- */
typedef struct {
    int n; const uint16_t *out;
    int cellof[MAXN], cellverts[MAXN][MAXN], cellsize[MAXN];
    int used[MAXN], pos[MAXN];
    u64 best; int have, cnt;
} CCtx;

static void canon_rec(CCtx *C, int q, u64 code) {
    int n = C->n;
    if (q == n) {
        if (!C->have || code < C->best) { C->best = code; C->have = 1; C->cnt = 1; }
        else if (code == C->best) C->cnt++;
        return;
    }
    int c = C->cellof[q];
    int nb = (q + 1) * q / 2, tot = n * (n - 1) / 2;
    for (int k = 0; k < C->cellsize[c]; k++) {
        int v = C->cellverts[c][k];
        if (C->used[v]) continue;
        u64 nc = code;
        for (int p = 0; p < q; p++) nc = (nc << 1) | ((C->out[C->pos[p]] >> v) & 1u);
        if (C->have && nc > (C->best >> (tot - nb))) continue;
        C->used[v] = 1; C->pos[q] = v;
        canon_rec(C, q + 1, nc);
        C->used[v] = 0;
    }
}

static u64 canon(int n, const uint16_t *out, int *aut) {
    CCtx C; memset(&C, 0, sizeof C);
    int col[MAXN], ncol, key[MAXN];
    for (int v = 0; v < n; v++) {
        int s = __builtin_popcount(out[v]);
        unsigned inv = (~out[v]) & ((1u << n) - 1) & ~(1u << v);
        int c3 = 0;
        for (int u = 0; u < n; u++) if ((out[v] >> u) & 1) c3 += __builtin_popcount(out[u] & inv);
        key[v] = s * 1000 + c3;
    }
    {
        int tmp[MAXN]; memcpy(tmp, key, sizeof(int) * n);
        for (int i = 1; i < n; i++) { int x = tmp[i], j = i - 1; while (j >= 0 && tmp[j] > x) { tmp[j + 1] = tmp[j]; j--; } tmp[j + 1] = x; }
        int d[MAXN], nd = 0;
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
    if (aut) *aut = C.cnt;
    return C.best;
}

static void decode(int n, u64 code, uint16_t *out) {
    int tot = n * (n - 1) / 2, b = tot - 1;
    for (int i = 0; i < n; i++) out[i] = 0;
    for (int q = 1; q < n; q++) for (int p = 0; p < q; p++) {
        int bit = (code >> b) & 1; b--;
        if (bit) out[p] |= 1u << q; else out[q] |= 1u << p;
    }
}
static void converse(int n, const uint16_t *out, uint16_t *op) {
    for (int i = 0; i < n; i++) op[i] = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if ((out[i] >> j) & 1) op[j] |= 1u << i;
}
static void delete_vertex(int n, const uint16_t *out, int x, uint16_t *res) {
    int k = 0;
    for (int i = 0; i < n; i++) {
        if (i == x) continue;
        uint16_t r = 0; int kk = 0;
        for (int j = 0; j < n; j++) { if (j == x) continue; if ((out[i] >> j) & 1) r |= 1u << kk; kk++; }
        res[k++] = r;
    }
}
/* dp buffer of (1<<n)*n entries supplied by caller */
static void ham_counts(int n, const uint16_t *out, u64 *hs, u64 *he, u64 *dp) {
    int full = (1 << n) - 1;
    uint16_t op[MAXN]; converse(n, out, op);
    for (int pass = 0; pass < 2; pass++) {
        const uint16_t *g = pass ? op : out;
        memset(dp, 0, sizeof(u64) * (size_t)(full + 1) * n);
        for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
        for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) {
            u64 c = dp[(size_t)S * n + v]; if (!c) continue;
            unsigned nb = g[v] & ~S;
            while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << w)) * n + w] += c; }
        }
        for (int v = 0; v < n; v++) (pass ? hs : he)[v] = dp[(size_t)full * n + v];
    }
}

/* ---------------- hash map (read-only after build) ---------------- */
typedef struct { u64 *key; int *val; size_t cap; } HMap;
static void hm_init(HMap *h, size_t cap) { h->cap = cap; h->key = malloc(sizeof(u64) * cap); h->val = malloc(sizeof(int) * cap); for (size_t i = 0; i < cap; i++) h->key[i] = ~0ULL; }
static size_t hm_slot(const HMap *h, u64 k) { u64 x = k * 0x9E3779B97F4A7C15ULL; size_t i = (size_t)(x >> 17) % h->cap; while (h->key[i] != ~0ULL && h->key[i] != k) { i++; if (i == h->cap) i = 0; } return i; }
static int hm_get(const HMap *h, u64 k) { size_t i = hm_slot(h, k); return h->key[i] == k ? h->val[i] : -1; }
static void hm_put(HMap *h, u64 k, int v) { size_t i = hm_slot(h, k); h->key[i] = k; h->val[i] = v; }

typedef struct { int nclass, nmerged, nsc; u64 *code; int *mid; HMap cls; } Order;
static Order ord_[MAXN + 1];

static int cmp_u64(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }

/* generate the sorted distinct canonical codes of order n from order n-1 (parallel) */
static u64 *generate(int n, int *count) {
    Order *P = &ord_[n - 1];
    size_t per = (size_t)1 << (n - 1), tot = (size_t)P->nclass * per;
    u64 *all = malloc(sizeof(u64) * (tot ? tot : 1));
    if (n == 1) { all[0] = 0; *count = 1; return all; }
#pragma omp parallel for schedule(dynamic, 16)
    for (int c = 0; c < P->nclass; c++) {
        uint16_t base[MAXN], t[MAXN];
        decode(n - 1, P->code[c], base);
        for (unsigned mask = 0; mask < per; mask++) {
            for (int i = 0; i < n - 1; i++) t[i] = base[i] | (((mask >> i) & 1) ? 0 : (1u << (n - 1)));
            t[n - 1] = (uint16_t)mask;
            all[(size_t)c * per + mask] = canon(n, t, NULL);
        }
    }
    qsort(all, tot, sizeof(u64), cmp_u64);
    size_t m = 0;
    for (size_t i = 0; i < tot; i++) if (i == 0 || all[i] != all[i - 1]) all[m++] = all[i];
    *count = (int)m;
    return realloc(all, sizeof(u64) * m);
}

static void build_small(int n) {   /* exact order with merged ids, for n <= 9 */
    Order *O = &ord_[n];
    int cnt; O->code = generate(n, &cnt); O->nclass = cnt;
    size_t cap = 4; while (cap < (size_t)cnt * 4) cap <<= 1;
    hm_init(&O->cls, cap);
    for (int i = 0; i < cnt; i++) hm_put(&O->cls, O->code[i], i);
    O->mid = malloc(sizeof(int) * cnt);
    int *conv = malloc(sizeof(int) * cnt);
#pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < cnt; i++) {
        uint16_t t[MAXN], op[MAXN];
        decode(n, O->code[i], t); converse(n, t, op);
        conv[i] = hm_get(&O->cls, canon(n, op, NULL));
    }
    int m = 0, sc = 0;
    for (int i = 0; i < cnt; i++) {
        if (conv[i] < 0) { fprintf(stderr, "converse missing\n"); exit(1); }
        if (conv[i] == i) sc++;
        if (conv[i] < i) O->mid[i] = O->mid[conv[i]]; else O->mid[i] = m++;
    }
    O->nmerged = m; O->nsc = sc;
    free(conv);
}


/* ---------------- S73: paths, cycles, Pfaffians ---------------- */

/* H and the number of directed Hamiltonian cycles (each counted once), n <= 16, adjacency as bitmasks */
static u64 big_paths(int n, const uint32_t *out, u64 *hc_out) {
    int full = (1 << n) - 1;
    u64 *dp = calloc((size_t)(full + 1) * n, sizeof(u64));
    for (int v = 0; v < n; v++) dp[(size_t)(1u << v) * n + v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) {
        u64 c = dp[(size_t)S * n + v]; if (!c) continue;
        uint32_t nb = out[v] & ~(uint32_t)S;
        while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << w)) * n + w] += c; }
    }
    u64 h = 0; for (int v = 0; v < n; v++) h += dp[(size_t)full * n + v];
    if (hc_out) {
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
        u64 hc = 0;
        if (n >= 3) for (int v = 1; v < n; v++) if (out[v] & 1u) hc += dp[(size_t)full * n + v];
        if (n == 1) hc = 0;
        *hc_out = hc;
    }
    free(dp);
    return h;
}

/* Pfaffian of the +-1 skew matrix of a tournament of even order n <= 10 (subset DP over even subsets) */
static long long pfaffian(int n, const uint16_t *out) {
    static __thread long long pf[1 << MAXN];
    int full = (1 << n) - 1;
    pf[0] = 1;
    for (int S = 1; S <= full; S++) {
        if (__builtin_popcount(S) & 1) continue;
        int i = __builtin_ctz(S);
        long long s = 0; int pos = 0;
        for (int j = i + 1; j < n; j++) {
            if (!((S >> j) & 1)) continue;
            pos++;   /* position of j inside S, counting i as position 0 */
            int a = ((out[i] >> j) & 1) ? 1 : -1;
            long long term = a * pf[S & ~(1 << i) & ~(1 << j)];
            s += (pos & 1) ? term : -term;   /* sign (-1)^(pos - 1) */
        }
        pf[S] = s;
    }
    return pf[full];
}

static void to32(int n, const uint16_t *t, uint32_t *o) { for (int i = 0; i < n; i++) o[i] = t[i]; }

/* ---------------- census mode ---------------- */
static void census(int NMAX) {
    for (int n = 1; n <= NMAX; n++) {
        double t0 = omp_get_wtime();
        u64 *codes; int cnt;
        if (n < NMAX || n <= 9) { build_small(n); codes = ord_[n].code; cnt = ord_[n].nclass; }
        else codes = generate(n, &cnt);
        if (n < 3) continue;
        u64 *H = malloc(sizeof(u64) * cnt), *HC = malloc(sizeof(u64) * cnt);
        char *sc = malloc(cnt); long long *PF = calloc(cnt, sizeof(long long));
#pragma omp parallel for schedule(dynamic, 256)
        for (int i = 0; i < cnt; i++) {
            uint16_t t[MAXN], op[MAXN]; uint32_t o[MAXN];
            decode(n, codes[i], t); to32(n, t, o);
            H[i] = big_paths(n, o, &HC[i]);
            converse(n, t, op);
            sc[i] = canon(n, op, NULL) == codes[i];
            if (!(n & 1)) PF[i] = pfaffian(n, t);
        }
        /* strong floor and single-cycle floor */
        u64 f = ~0ULL; for (int i = 0; i < cnt; i++) if (HC[i] >= 1 && H[i] < f) f = H[i];
        long long nfloor = 0, nfloor_hc1 = 0;
        for (int i = 0; i < cnt; i++) if (HC[i] >= 1 && H[i] == f) { nfloor++; if (HC[i] == 1) nfloor_hc1++; }
        /* cycle ceiling */
        long long minslack = (long long)1 << 62, neq = 0, nviol = 0, nstrong = 0;
        for (int i = 0; i < cnt; i++) if (HC[i] >= 1) {
            nstrong++;
            long long s = (long long)H[i] - (long long)f - 2 * ((long long)HC[i] - 1);
            if (s < minslack) minslack = s;
            if (s == 0) neq++;
            if (s < 0) { if (nviol < 5) printf("    CEILING VIOLATION n=%d code=%llu H=%llu hc=%llu\n", n, codes[i], H[i], HC[i]); nviol++; }
        }
        printf("n=%d classes=%d strong=%lld | strong floor f=%llu: %lld floor classes, %lld with exactly one Hamiltonian cycle | cycle ceiling H >= f + 2(hc-1): min slack %lld, equality classes %lld, violations %lld\n",
               n, cnt, nstrong, f, nfloor, nfloor_hc1, minslack, neq, nviol);
        if (neq) {
            printf("   equality cases (H, hc):");
            int shown = 0;
            for (int i = 0; i < cnt && shown < 12; i++) if (HC[i] >= 1 && (long long)H[i] - (long long)f - 2 * ((long long)HC[i] - 1) == 0) { printf(" (%llu,%llu)%s", H[i], HC[i], sc[i] ? "sc" : ""); shown++; }
            printf("\n");
        }
        /* the uniquely Hamiltonian stratum (Douglas's class; Garey: F(2n-6) classes) and the single-cycle band */
        {
            long long nu = 0; u64 umax = 0, multimin = ~0ULL;
            for (int i = 0; i < cnt; i++) {
                if (HC[i] == 1) { nu++; if (H[i] > umax) umax = H[i]; }
                if (HC[i] >= 2 && H[i] < multimin) multimin = H[i];
            }
            char *u = calloc(umax + 2, 1), *s = calloc(umax + 2, 1);
            for (int i = 0; i < cnt; i++) { if (HC[i] == 1) u[H[i]] = 1; if (HC[i] >= 1 && H[i] <= umax) s[H[i]] = 1; }
            long long nd = 0; for (u64 v = 1; v <= umax; v += 2) nd += u[v];
            printf("   uniquely Hamiltonian stratum: %lld classes, %lld distinct H in [%llu, %llu]; single-cycle band [f, %llu): every strong class there has one Hamiltonian cycle\n",
                   nu, nd, f, umax, multimin);
            printf("   odd values in [f, %llu] attained by strong classes but by no uniquely Hamiltonian one:", multimin < umax ? multimin : umax);
            int k = 0; for (u64 v = f; v <= umax && v <= multimin; v += 2) if (s[v] && !u[v]) { if (k < 12) printf(" %llu", v); k++; }
            printf(" (%d); odd values in [f, %llu] attained by no strong class:", k, multimin < umax ? multimin : umax);
            k = 0; for (u64 v = f; v <= umax && v <= multimin; v += 2) if (!s[v]) { if (k < 12) printf(" %llu", v); k++; }
            printf(" (%d)\n", k);
            /* is every uniquely Hamiltonian class a one-arc join (a vertex cut A | V-A crossed backward by exactly one arc)? */
            long long njoin = 0;
            for (int i = 0; i < cnt; i++) {
                if (HC[i] != 1) continue;
                uint16_t t[MAXN]; decode(n, codes[i], t);
                int ok = 0;
                for (int A = 1; A < (1 << n) - 1 && !ok; A++) {
                    int back = 0;
                    for (int y = 0; y < n && back < 2; y++) if (!((A >> y) & 1)) back += __builtin_popcount(t[y] & A);
                    if (back == 1) ok = 1;
                }
                njoin += ok;
            }
            printf("   one-arc join test: %lld of %lld uniquely Hamiltonian classes have a cut crossed backward by exactly one arc\n", njoin, nu);
            if (n <= 8) {
                printf("   H-values of the uniquely Hamiltonian classes:");
                for (u64 v = 1; v <= umax; v += 2) if (u[v]) printf(" %llu", v);
                printf("\n");
            }
            free(u); free(s);
        }
        /* k-cycle floors */
        printf("   k-cycle floors f_k (k: min H [classes at min]):");
        for (int k = 1; k <= 16; k++) {
            u64 m = ~0ULL; long long c = 0;
            for (int i = 0; i < cnt; i++) if (HC[i] == (u64)k) { if (H[i] < m) { m = H[i]; c = 1; } else if (H[i] == m) c++; }
            if (m != ~0ULL) printf(" %d:%llu[%lld]", k, m, c);
        }
        printf("\n");
        /* spectrum gaps: all classes and self-converse classes */
        u64 hmax = 0; for (int i = 0; i < cnt; i++) if (H[i] > hmax) hmax = H[i];
        char *seen = calloc(hmax + 2, 1), *seensc = calloc(hmax + 2, 1);
        u64 hmaxsc = 0;
        for (int i = 0; i < cnt; i++) { seen[H[i]] = 1; if (sc[i]) { seensc[H[i]] = 1; if (H[i] > hmaxsc) hmaxsc = H[i]; } }
        printf("   odd gaps of the order-%d spectrum below its max %llu:", n, hmax);
        int g = 0; for (u64 v = 1; v <= hmax; v += 2) if (!seen[v]) { if (g < 16) printf(" %llu", v); g++; }
        printf(" (%d)\n", g);
        printf("   self-converse (real-point) spectrum: max %llu; odd gaps below 200:", hmaxsc);
        g = 0; for (u64 v = 1; v <= 200 && v <= hmaxsc; v += 2) if (!seensc[v]) { printf(" %llu", v); g++; }
        printf(" (%d below 200)\n", g);
        free(seen); free(seensc);
        /* Pfaffian statistics, even n */
        if (!(n & 1)) {
            long long tab4[4][4] = {{0}}, agree4 = 0, agree8 = 0, aggr_hc = 0;
            for (int i = 0; i < cnt; i++) {
                long long p = PF[i] < 0 ? -PF[i] : PF[i];
                tab4[p % 4][H[i] % 4]++;
                if ((u64)(p % 4) == H[i] % 4) agree4++;
                if ((u64)(p % 8) == H[i] % 8) agree8++;
                if ((p % 4 == 1) == (HC[i] % 2 == 0)) aggr_hc++;
            }
            printf("   |Pf| vs H (n even): |Pf| = H mod 4 in %lld/%d classes, mod 8 in %lld/%d; table [|Pf| mod 4][H mod 4]: 1->(%lld,%lld) 3->(%lld,%lld)\n",
                   agree4, cnt, agree8, cnt, tab4[1][1], tab4[1][3], tab4[3][1], tab4[3][3]);
            printf("   |Pf| mod 4 == 1 iff hc even: %lld/%d\n", aggr_hc, cnt);
            if (n <= 6) {
                printf("   (H, hc, |Pf|) by class:");
                for (int i = 0; i < cnt; i++) printf(" (%llu,%llu,%lld)", H[i], HC[i], PF[i] < 0 ? -PF[i] : PF[i]);
                printf("\n");
            }
        }
        printf("   [%.1f s]\n", omp_get_wtime() - t0);
        fflush(stdout);
        free(H); free(HC); free(sc); free(PF);
        if (n == NMAX && !(n < NMAX || n <= 9)) free(codes);
    }
}

/* ---------------- chains mode ---------------- */
static int g_parts[32], g_np, g_bad, g_count;
static void chain_check(int N) {
    int n = N + 1; uint32_t o[32];
    for (int i = 0; i < n; i++) { o[i] = 0; for (int j = i + 1; j < n; j++) o[i] |= 1u << j; }
    int l = 0; u64 prod = 1;
    for (int k = 0; k < g_np; k++) {
        int l2 = l + g_parts[k];
        o[l] &= ~(1u << l2); o[l2] |= 1u << l;     /* back arc x_(l2) -> x_l */
        prod *= ((u64)1 << (g_parts[k] - 1)) + 1;
        l = l2;
    }
    u64 hc; u64 h = big_paths(n, o, &hc);
    g_count++;
    if (h != prod || hc != 1) { g_bad++; if (g_bad <= 5) { printf("   FAIL parts"); for (int k = 0; k < g_np; k++) printf(" %d", g_parts[k]); printf(": H=%llu prod=%llu hc=%llu\n", h, prod, hc); } }
}
static void compositions(int rem, int N) {
    if (rem == 0) { chain_check(N); return; }
    for (int p = 2; p <= rem; p++) { g_parts[g_np++] = p; compositions(rem - p, N); g_np--; }
}
static void chains(int NMAX) {
    for (int n = 3; n <= NMAX; n++) {
        g_bad = 0; g_count = 0;
        compositions(n - 1, n - 1);
        printf("order %d: %d interval chains; H = prod(2^(p-1)+1) and exactly one Hamiltonian cycle in all but %d\n", n, g_count, g_bad);
        fflush(stdout);
    }
}

/* ---------------- joins mode ---------------- */
typedef struct { u64 h, R, e; int ord; } Trip;
static int cmp_trip(const void *a, const void *b) {
    const Trip *x = a, *y = b;
    if (x->h != y->h) return x->h < y->h ? -1 : 1;
    if (x->R != y->R) return x->R < y->R ? -1 : 1;
    if (x->e != y->e) return x->e < y->e ? -1 : 1;
    return x->ord - y->ord;
}
/* rooted triples (h, R_v, e_v) of a tournament of order n */
static void triples(int n, const uint16_t *t, u64 *h, u64 *R, u64 *e) {
    int full = (1 << n) - 1;
    static __thread u64 dE[(1 << MAXN) * MAXN], dS[(1 << MAXN) * MAXN], p[1 << MAXN];
    uint16_t op[MAXN]; converse(n, t, op);
    for (int pass = 0; pass < 2; pass++) {
        const uint16_t *g = pass ? op : t; u64 *dp = pass ? dS : dE;
        memset(dp, 0, sizeof(u64) * (size_t)(full + 1) * n);
        for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
        for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) {
            u64 c = dp[(size_t)S * n + v]; if (!c) continue;
            unsigned nb = g[v] & ~S;
            while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << w)) * n + w] += c; }
        }
    }
    p[0] = 1;
    for (int S = 1; S <= full; S++) { u64 s = 0; for (int v = 0; v < n; v++) s += dE[(size_t)S * n + v]; p[S] = s; }
    *h = p[full];
    for (int v = 0; v < n; v++) {
        e[v] = dE[(size_t)full * n + v];
        u64 r = 0;
        for (int S = 1; S <= full; S++) if ((S >> v) & 1) r += p[full & ~S] * dS[(size_t)S * n + v];   /* paths of T[S] starting at v */
        R[v] = r;
    }
}

static void joins(int KA, int KB, long VMAX) {
    int K = KA > KB ? KA : KB;
    Trip *tr = NULL; size_t ntr = 0, cap = 0;
    for (int n = 1; n <= K; n++) {
        build_small(n);
        Order *O = &ord_[n];
        size_t base = ntr; ntr += (size_t)O->nclass * n;
        if (ntr > cap) { cap = ntr * 2; tr = realloc(tr, sizeof(Trip) * cap); }
#pragma omp parallel for schedule(dynamic, 64)
        for (int c = 0; c < O->nclass; c++) {
            uint16_t t[MAXN]; u64 h, R[MAXN], e[MAXN];
            decode(n, O->code[c], t); triples(n, t, &h, R, e);
            for (int v = 0; v < n; v++) { Trip *x = &tr[base + (size_t)c * n + v]; x->h = h; x->R = R[v]; x->e = e[v]; x->ord = n; }
        }
    }
    qsort(tr, ntr, sizeof(Trip), cmp_trip);
    size_t m = 0;
    for (size_t i = 0; i < ntr; i++) if (i == 0 || tr[i].h != tr[m - 1].h || tr[i].R != tr[m - 1].R || tr[i].e != tr[m - 1].e) tr[m++] = tr[i];
    printf("joins KA=%d KB=%d VMAX=%ld: %zu distinct rooted triples (h, R, e) of orders <= %d (each kept at its least order)\n", KA, KB, VMAX, m, K);
    /* verification of the join formula by direct counting on random rooted pairs */
    {
        unsigned long long seed = 20261104ULL; int bad = 0, tested = 0;
        for (int it = 0; it < 3000; it++) {
            seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
            int na = 1 + (int)((seed >> 33) % KA), nb = 1 + (int)((seed >> 45) % KB);
            Order *OA = &ord_[na], *OB = &ord_[nb];
            seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
            int ca = (int)((seed >> 20) % OA->nclass), cb = (int)((seed >> 40) % OB->nclass);
            int a = (int)((seed >> 8) % na), b = (int)((seed >> 14) % nb);
            uint16_t A[MAXN], B[MAXN], Bop[MAXN]; decode(na, OA->code[ca], A); decode(nb, OB->code[cb], B); converse(nb, B, Bop);
            u64 hA, RA[MAXN], eA[MAXN], hB, RB[MAXN], eB[MAXN];
            triples(na, A, &hA, RA, eA); triples(nb, Bop, &hB, RB, eB);    /* y(B, b) = x(B^op, b) */
            u64 pred = hA * hB + RA[a] * RB[b] - eA[a] * eB[b];
            int n = na + nb; uint32_t o[32];
            for (int i = 0; i < na; i++) o[i] = A[i] | (((1u << nb) - 1) << na);
            for (int j = 0; j < nb; j++) o[na + j] = (uint32_t)B[j] << na;
            o[a] &= ~(1u << (na + b)); o[na + b] |= 1u << a;               /* reverse a -> b into b -> a */
            u64 direct = big_paths(n, o, NULL);
            tested++; if (direct != pred) { bad++; if (bad <= 3) printf("   formula mismatch: orders %d,%d pred %llu direct %llu\n", na, nb, pred, direct); }
        }
        printf("   join formula H(J) = x_h y_h + x_R y_R - x_e y_e verified by direct counting on %d random joins: %d mismatches\n", tested, bad);
    }
    long *r = calloc(VMAX + 1, sizeof(long));
    size_t na_ = 0, nb_ = 0;
    for (size_t i = 0; i < m; i++) { if (tr[i].ord <= KA) na_++; if (tr[i].ord <= KB) nb_++; }
#pragma omp parallel
    {
        long *rl = calloc(VMAX + 1, sizeof(long));
#pragma omp for schedule(dynamic, 64)
        for (long i = 0; i < (long)m; i++) {
            if (tr[i].ord > KA) continue;
            for (size_t j = 0; j < m; j++) {
                if (tr[j].ord > KB) continue;
                long long N = (long long)(tr[i].h * tr[j].h + tr[i].R * tr[j].R) - (long long)(tr[i].e * tr[j].e);
                if (N >= 1 && N <= VMAX) rl[N]++;
            }
        }
#pragma omp critical
        for (long v = 0; v <= VMAX; v++) r[v] += rl[v];
        free(rl);
    }
    printf("   pairs: %zu A-triples (order <= %d) x %zu B-triples (order <= %d)\n", na_, KA, nb_, KB);
    long ng = 0; printf("   odd N <= %ld with no representation:", VMAX);
    for (long v = 1; v <= VMAX; v += 2) if (!r[v]) { if (ng < 40) printf(" %ld", v); ng++; }
    printf(" (%ld total)\n", ng);
    long even = 0; for (long v = 2; v <= VMAX; v += 2) even += r[v];
    printf("   representations landing on even N: %ld (must be 0: every H is odd)\n", even);
    for (long lo = 1; lo * 2 <= VMAX + 1; lo *= 2) {
        long hi = lo * 2, mn = -1, cntv = 0; double sum = 0;
        for (long v = lo | 1; v < hi && v <= VMAX; v += 2) { if (mn < 0 || r[v] < mn) mn = r[v]; sum += r[v]; cntv++; }
        printf("   window [%ld, %ld): min r = %ld, mean r = %.1f\n", lo, hi, mn, cntv ? sum / cntv : 0.0);
    }
    /* residue bias in a saturated window */
    long lo = VMAX / 8, hi = VMAX / 4;
    int mods[4] = {3, 5, 7, 9};
    for (int q = 0; q < 4; q++) {
        int md = mods[q]; double s[16] = {0}; long c[16] = {0}; double tot = 0; long tc = 0;
        for (long v = lo | 1; v < hi; v += 2) { s[v % md] += r[v]; c[v % md]++; tot += r[v]; tc++; }
        printf("   mean r(N)/overall on [%ld, %ld) by N mod %d:", lo, hi, md);
        for (int k = 0; k < md; k++) printf(" %d:%.3f", k, c[k] ? (s[k] / c[k]) / (tot / tc) : 0.0);
        printf("\n");
    }
    free(r); free(tr);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: census NMAX | chains NMAX | joins KA KB VMAX\n"); return 1; }
    printf("scf_crevasses_20261104 mode=%s threads=%d\n", argv[1], omp_get_max_threads());
    if (!strcmp(argv[1], "census")) census(argc > 2 ? atoi(argv[2]) : 9);
    else if (!strcmp(argv[1], "chains")) chains(argc > 2 ? atoi(argv[2]) : 16);
    else if (!strcmp(argv[1], "joins")) joins(atoi(argv[2]), atoi(argv[3]), atol(argv[4]));
    return 0;
}
