/*
 * audit_census.c  -- independent audit code for THM-4866 / HYP-9457 (auditor, 2026-11-04).
 * Written from scratch; shares no code with scf_crevasses_20261104.c.
 *
 * Canonical form: individualisation-refinement (ordered partitions refined by out-neighbour colour counts,
 * first non-singleton cell individualised, minimum code over all leaves).
 * Paths: subset DP. Cycles: subset DP from a fixed start vertex. |Pf|: sqrt(det) of the +-1 skew matrix via exact
 * fraction-free (Bareiss) elimination.
 *
 * Mode "census NMAX" (NMAX <= 9): full class census, orders 1..NMAX.
 * Mode "ears": order-10 statements via one-vertex ears of the order-9 census (needs census through 9 in memory).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define MAXN 12
typedef unsigned long long u64;
typedef __int128 i128;

/* ---------- canonical form by individualisation-refinement ---------- */
static void refine(int n, const uint16_t *adj, int *col, int *k) {
    for (;;) {
        int sig[MAXN][MAXN + 1];
        for (int v = 0; v < n; v++) {
            for (int j = 0; j <= *k; j++) sig[v][j] = 0;
            sig[v][0] = col[v];
            for (int u = 0; u < n; u++) if ((adj[v] >> u) & 1) sig[v][1 + col[u]]++;
        }
        int ord[MAXN];
        for (int i = 0; i < n; i++) ord[i] = i;
        /* insertion sort by signature */
        for (int i = 1; i < n; i++) {
            int x = ord[i], j = i - 1;
            while (j >= 0 && memcmp(sig[ord[j]], sig[x], sizeof(int) * (*k + 1)) > 0) { ord[j + 1] = ord[j]; j--; }
            ord[j + 1] = x;
        }
        int nc = 0, newc[MAXN];
        for (int i = 0; i < n; i++) {
            if (i > 0 && memcmp(sig[ord[i]], sig[ord[i - 1]], sizeof(int) * (*k + 1)) != 0) nc++;
            newc[ord[i]] = nc;
        }
        nc++;
        memcpy(col, newc, sizeof(int) * n);
        if (nc == *k) return;
        *k = nc;
    }
}

static u64 leaf_code(int n, const uint16_t *adj, const int *col) {
    int pos[MAXN];
    for (int v = 0; v < n; v++) pos[col[v]] = v;
    u64 c = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) c = (c << 1) | ((adj[pos[i]] >> pos[j]) & 1);
    return c;
}

static void ir_search(int n, const uint16_t *adj, int *col, int k, u64 *best, int *have) {
    refine(n, adj, col, &k);
    if (k == n) {
        u64 c = leaf_code(n, adj, col);
        if (!*have || c < *best) { *best = c; *have = 1; }
        return;
    }
    /* first non-singleton cell */
    int cnt[MAXN] = {0};
    for (int v = 0; v < n; v++) cnt[col[v]]++;
    int c = 0; while (cnt[c] == 1) c++;
    for (int v = 0; v < n; v++) {
        if (col[v] != c) continue;
        int col2[MAXN];
        for (int u = 0; u < n; u++) col2[u] = col[u] > c ? col[u] + 1 : (col[u] == c && u != v ? c + 1 : col[u]);
        ir_search(n, adj, col2, k + 1, best, have);
    }
}

static u64 canon(int n, const uint16_t *adj) {
    int col[MAXN];
    for (int v = 0; v < n; v++) col[v] = 0;
    u64 best = 0; int have = 0;
    ir_search(n, adj, col, 1, &best, &have);
    return best;
}

static void decode(int n, u64 code, uint16_t *adj) {
    int b = n * (n - 1) / 2 - 1;
    for (int i = 0; i < n; i++) adj[i] = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        if ((code >> b) & 1) adj[i] |= 1u << j; else adj[j] |= 1u << i;
        b--;
    }
}

static void conv(int n, const uint16_t *adj, uint16_t *op) {
    for (int i = 0; i < n; i++) op[i] = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if ((adj[i] >> j) & 1) op[j] |= 1u << i;
}

/* ---------- paths and cycles ---------- */
static u64 count_paths(int n, const uint16_t *adj) {
    static __thread u64 dp[(1 << MAXN) * MAXN];
    int full = (1 << n) - 1;
    memset(dp, 0, sizeof(u64) * (size_t)(full + 1) * n);
    for (int v = 0; v < n; v++) dp[(1 << v) * n + v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) {
        u64 c = dp[S * n + v]; if (!c) continue;
        int nb = adj[v] & ~S;
        for (int w = 0; w < n; w++) if ((nb >> w) & 1) dp[(S | (1 << w)) * n + w] += c;
    }
    u64 h = 0; for (int v = 0; v < n; v++) h += dp[full * n + v];
    return h;
}
static u64 count_cycles(int n, const uint16_t *adj) {
    if (n < 3) return 0;
    static __thread u64 dp[(1 << MAXN) * MAXN];
    int full = (1 << n) - 1;
    memset(dp, 0, sizeof(u64) * (size_t)(full + 1) * n);
    dp[1 * n + 0] = 1;
    for (int S = 1; S <= full; S += 2) for (int v = 0; v < n; v++) {
        u64 c = dp[S * n + v]; if (!c) continue;
        int nb = adj[v] & ~S;
        for (int w = 0; w < n; w++) if ((nb >> w) & 1) dp[(S | (1 << w)) * n + w] += c;
    }
    u64 hc = 0; for (int v = 1; v < n; v++) if (adj[v] & 1) hc += dp[full * n + v];
    return hc;
}

/* |Pf| via exact Bareiss determinant of the +-1 skew matrix (det = Pf^2) */
static long long abs_pf(int n, const uint16_t *adj) {
    i128 M[MAXN][MAXN];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) M[i][j] = (i == j) ? 0 : (((adj[i] >> j) & 1) ? 1 : -1);
    i128 prev = 1; int sign = 1;
    for (int k = 0; k < n - 1; k++) {
        if (M[k][k] == 0) {
            int r = -1; for (int i = k + 1; i < n; i++) if (M[i][k] != 0) { r = i; break; }
            if (r < 0) return 0;
            for (int j = 0; j < n; j++) { i128 t = M[k][j]; M[k][j] = M[r][j]; M[r][j] = t; }
            sign = -sign;
        }
        for (int i = k + 1; i < n; i++) for (int j = k + 1; j < n; j++) M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) / prev;
        prev = M[k][k];
    }
    i128 det = M[n - 1][n - 1] * sign;
    long long d = (long long)det;
    if (d < 0) { fprintf(stderr, "negative det?\n"); exit(1); }
    long long r = 0; while ((r + 1) * (r + 1) <= d) r++;
    if (r * r != d) { fprintf(stderr, "det not a square\n"); exit(1); }
    return r;
}

/* strong components: returns number of components, comp[v] = index in condensation order (0 = top) */
static int strong_components(int n, const uint16_t *adj, int *comp) {
    /* reachability closure */
    uint16_t R[MAXN];
    for (int v = 0; v < n; v++) R[v] = adj[v] | (1u << v);
    for (int it = 0; it < n; it++) for (int v = 0; v < n; v++) { uint16_t r = R[v]; for (int u = 0; u < n; u++) if ((R[v] >> u) & 1) r |= R[u]; R[v] = r; }
    int nc = 0, id[MAXN]; for (int v = 0; v < n; v++) id[v] = -1;
    for (int v = 0; v < n; v++) if (id[v] < 0) { for (int u = 0; u < n; u++) if (((R[v] >> u) & 1) && ((R[u] >> v) & 1)) id[u] = nc; nc++; }
    /* order components: component reaching more vertices comes first */
    int reach[MAXN] = {0};
    for (int c = 0; c < nc; c++) { for (int v = 0; v < n; v++) if (id[v] == c) { reach[c] = __builtin_popcount(R[v]); break; } }
    for (int v = 0; v < n; v++) { int r = 0; for (int c = 0; c < nc; c++) if (reach[c] > reach[id[v]]) r++; comp[v] = r; }
    return nc;
}

static int one_arc_cut(int n, const uint16_t *adj) {
    int full = (1 << n) - 1;
    for (int A = 1; A < full; A++) {
        int back = 0;
        for (int y = 0; y < n; y++) if (!((A >> y) & 1)) back += __builtin_popcount(adj[y] & A);
        if (back == 1) return 1;
    }
    return 0;
}

/* ---------- census storage ---------- */
static u64 *cls[MAXN + 1]; static int ncls[MAXN + 1];
static int cmpu(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }

static void build(int n) {
    if (n == 1) { cls[1] = malloc(8); cls[1][0] = 0; ncls[1] = 1; return; }
    int m = ncls[n - 1]; size_t per = (size_t)1 << (n - 1);
    u64 *all = malloc(sizeof(u64) * m * per);
#pragma omp parallel for schedule(dynamic, 8)
    for (int c = 0; c < m; c++) {
        uint16_t a[MAXN]; decode(n - 1, cls[n - 1][c], a);
        for (size_t mask = 0; mask < per; mask++) {
            uint16_t t[MAXN];
            /* new vertex n-1: bit i of mask set => n-1 -> i, else i -> n-1 */
            for (int i = 0; i < n - 1; i++) t[i] = a[i] | (((mask >> i) & 1) ? 0 : (1u << (n - 1)));
            t[n - 1] = (uint16_t)mask;
            all[(size_t)c * per + mask] = canon(n, t);
        }
    }
    qsort(all, (size_t)m * per, sizeof(u64), cmpu);
    size_t k = 0;
    for (size_t i = 0; i < (size_t)m * per; i++) if (i == 0 || all[i] != all[i - 1]) all[k++] = all[i];
    cls[n] = realloc(all, sizeof(u64) * k); ncls[n] = (int)k;
}

static u64 *Hs[MAXN + 1], *HCs[MAXN + 1];

static void census_order(int n, int fnext) {
    int m = ncls[n];
    u64 *H = malloc(8 * (size_t)m), *HC = malloc(8 * (size_t)m);
    char *sc = calloc(m, 1), *cut = calloc(m, 1), *mino = calloc(m, 1), *mini = calloc(m, 1);
    long long *pf = calloc(m, sizeof(long long));
    char *dec_ok = calloc(m, 1);
#pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < m; i++) {
        uint16_t t[MAXN], op[MAXN]; decode(n, cls[n][i], t);
        H[i] = count_paths(n, t); HC[i] = count_cycles(n, t);
        conv(n, t, op); sc[i] = canon(n, op) == cls[n][i];
        int mo = n, mi = n;
        for (int v = 0; v < n; v++) { int s = __builtin_popcount(t[v]); if (s < mo) mo = s; if (n - 1 - s < mi) mi = n - 1 - s; }
        mino[i] = (char)mo; mini[i] = (char)mi;
        if (HC[i] == 1) cut[i] = (char)one_arc_cut(n, t);
        if (!(n & 1)) pf[i] = abs_pf(n, t);
        if (sc[i]) {
            /* real-point decomposition: components pair C_i ~ C_(k+1-i)^op; H = H(M) * s^2 */
            int comp[MAXN]; int k = strong_components(n, t, comp);
            u64 prod = 1, hpart[MAXN]; int sz[MAXN] = {0};
            for (int c = 0; c < k; c++) {
                uint16_t sub[MAXN]; int idx[MAXN], q = 0;
                for (int v = 0; v < n; v++) if (comp[v] == c) idx[q++] = v;
                for (int x = 0; x < q; x++) { sub[x] = 0; for (int y = 0; y < q; y++) if ((t[idx[x]] >> idx[y]) & 1) sub[x] |= 1u << y; }
                hpart[c] = count_paths(q, sub); sz[c] = q; prod *= hpart[c];
            }
            int ok = (prod == H[i]);
            for (int c = 0; c < k; c++) if (hpart[c] != hpart[k - 1 - c] || sz[c] != sz[k - 1 - c]) ok = 0;
            int msz = (k & 1) ? sz[k / 2] : 0;
            if ((msz & 1) != (n & 1)) ok = 0;
            dec_ok[i] = (char)ok;
        }
    }
    Hs[n] = H; HCs[n] = HC;
    if (n < 3) return;
    u64 f = ~0ULL; for (int i = 0; i < m; i++) if (HC[i] && H[i] < f) f = H[i];
    int nf = 0, nf1 = 0; for (int i = 0; i < m; i++) if (HC[i] && H[i] == f) { nf++; if (HC[i] == 1) nf1++; }
    long long minsl = 1LL << 60; int neq = 0, nviol = 0, nstrong = 0, eqnonfloor = 0;
    for (int i = 0; i < m; i++) if (HC[i]) {
        nstrong++;
        long long s = (long long)H[i] - (long long)f - 2 * ((long long)HC[i] - 1);
        if (s < minsl) minsl = s; if (s == 0) { neq++; if (H[i] != f) eqnonfloor++; } if (s < 0) nviol++;
    }
    printf("n=%d classes=%d strong=%d f=%llu floor classes=%d (hc=1: %d) ceiling: min slack %lld, equality %d (non-floor %d), violations %d\n",
           n, m, nstrong, f, nf, nf1, minsl, neq, eqnonfloor, nviol);
    u64 band = ~0ULL; for (int i = 0; i < m; i++) if (HC[i] >= 2 && H[i] < band) band = H[i];
    int nu = 0, ncut = 0, nsole = 0, nsoleo = 0, nsolei = 0; u64 umax = 0;
    for (int i = 0; i < m; i++) if (HC[i] == 1) {
        nu++; ncut += cut[i]; if (H[i] > umax) umax = H[i];
        if (mino[i] == 1 || mini[i] == 1) nsole++;
        if (mino[i] == 1) nsoleo++;
        if (mini[i] == 1) nsolei++;
    }
    char *uh = calloc(umax + 2, 1); int nd = 0;
    for (int i = 0; i < m; i++) if (HC[i] == 1) uh[H[i]] = 1;
    for (u64 v = 0; v <= umax; v++) nd += uh[v];
    printf("   f_2-band end min{H: hc>=2} = %llu; U(n): %d classes, %d distinct H, range [%llu, %llu]\n", band, nu, nd, f, umax);
    printf("   U(n) classes with a one-arc cut: %d / %d; with a vertex of out-degree 1 or in-degree 1: %d (out-deg 1: %d, in-deg 1: %d)\n",
           ncut, nu, nsole, nsoleo, nsolei);
    if (n <= 8) { printf("   U(n) values:"); for (u64 v = 0; v <= umax; v++) if (uh[v]) printf(" %llu", v); printf("\n"); }
    /* band check: every strong class with H < band has hc = 1 (by definition), list band holes */
    u64 hmax = 0; for (int i = 0; i < m; i++) if (H[i] > hmax) hmax = H[i];
    char *st = calloc(hmax + 2, 1), *all = calloc(hmax + 2, 1), *scs = calloc(hmax + 2, 1);
    for (int i = 0; i < m; i++) { all[H[i]] = 1; if (HC[i]) st[H[i]] = 1; if (sc[i]) scs[H[i]] = 1; }
    printf("   odd values in [f, band) missed by strong classes:"); int bh = 0;
    for (u64 v = f; v < band && v <= hmax; v += 2) if (!st[v]) { printf(" %llu", v); bh++; } printf(" (%d)\n", bh);
    if (fnext > 0) {
        printf("   odd values in [f(n), f(n+1)) = [%llu, %d) missed by strong classes:", f, fnext); int hh = 0;
        for (u64 v = f; v < (u64)fnext; v += 2) if (v > hmax || !st[v]) { printf(" %llu", v); hh++; } printf(" (%d)\n", hh);
    }
    printf("   k-cycle floors:");
    for (int k = 1; k <= 6; k++) { u64 mm = ~0ULL; int c = 0; for (int i = 0; i < m; i++) if (HC[i] == (u64)k) { if (H[i] < mm) { mm = H[i]; c = 1; } else if (H[i] == mm) c++; } if (mm != ~0ULL) printf(" f_%d=%llu[%d]", k, mm, c); }
    printf("\n");
    /* self-converse */
    int nsc = 0, ndec = 0, bad5 = 0, bad3 = 0;
    for (int i = 0; i < m; i++) if (sc[i]) { nsc++; ndec += dec_ok[i]; if (H[i] == 5 && (n & 1)) bad5++; if (H[i] == 3 && !(n & 1)) bad3++; }
    printf("   self-converse classes %d; decomposition H = H(M) s^2 with palindromic components verified in %d; H=5 at odd n: %d, H=3 at even n: %d\n", nsc, ndec, bad5, bad3);
    printf("   self-converse odd gaps below 200:"); int g = 0;
    for (u64 v = 1; v < 200; v += 2) if (v > hmax || !scs[v]) { g++; } printf(" %d\n", g);
    if (!(n & 1)) {
        long long ag4 = 0, aghc = 0;
        for (int i = 0; i < m; i++) { if ((u64)(pf[i] % 4) == H[i] % 4) ag4++; if ((pf[i] % 4 == 1) == (HC[i] % 2 == 0)) aghc++; }
        printf("   |Pf| = H mod 4 in %lld/%d; (|Pf| = 1 mod 4) iff hc even in %lld/%d\n", ag4, m, aghc, m);
    }
    free(st); free(all); free(scs); free(uh);
    free(sc); free(cut); free(mino); free(mini); free(pf); free(dec_ok);
    fflush(stdout);
}

/* ---------- order 10 via ears ---------- */
static int cmp_pair(const void *a, const void *b) { const u64 *x = a, *y = b; return x[0] < y[0] ? -1 : x[0] > y[0] ? 1 : (x[1] < y[1] ? -1 : x[1] > y[1]); }

static void ears10(void) {
    const int n = 10; const u64 LIM = 205;  /* all strong 10-classes with H <= 205; plus all uniquely Hamiltonian */
    int m = ncls[9];
    size_t cap = 1 << 20, cnt = 0; u64 *rec = malloc(sizeof(u64) * 3 * cap);   /* (canon, H, hc) */
#pragma omp parallel
    {
        size_t lcap = 1 << 14, lc = 0; u64 *loc = malloc(sizeof(u64) * 3 * lcap);
#pragma omp for schedule(dynamic, 16)
        for (int c = 0; c < m; c++) {
            if (!HCs[9][c]) continue;
            int uniq = HCs[9][c] == 1;
            if (Hs[9][c] > LIM && !uniq) continue;
            uint16_t a[MAXN]; decode(9, cls[9][c], a);
            for (int mask = 0; mask < 512; mask++) {
                uint16_t t[MAXN];
                for (int i = 0; i < 9; i++) t[i] = a[i] | (((mask >> i) & 1) ? 0 : (1u << 9));
                t[9] = (uint16_t)mask;
                u64 h = count_paths(n, t);
                if (h > LIM && !uniq) continue;
                u64 hc = count_cycles(n, t);
                if (!hc) continue;
                if (!(h <= LIM || hc == 1)) continue;
                if (lc == lcap) { lcap *= 2; loc = realloc(loc, sizeof(u64) * 3 * lcap); }
                loc[3 * lc] = canon(n, t); loc[3 * lc + 1] = h; loc[3 * lc + 2] = hc; lc++;
            }
        }
#pragma omp critical
        {
            if (cnt + lc > cap) { while (cnt + lc > cap) cap *= 2; rec = realloc(rec, sizeof(u64) * 3 * cap); }
            memcpy(rec + 3 * cnt, loc, sizeof(u64) * 3 * lc); cnt += lc;
        }
        free(loc);
    }
    qsort(rec, cnt, sizeof(u64) * 3, cmp_pair);
    size_t k = 0;
    for (size_t i = 0; i < cnt; i++) if (i == 0 || rec[3 * i] != rec[3 * (k - 1)]) { memmove(rec + 3 * k, rec + 3 * i, 24); k++; }
    /* report */
    u64 f = ~0ULL; for (size_t i = 0; i < k; i++) if (rec[3 * i + 1] < f) f = rec[3 * i + 1];
    int nf = 0, nf1 = 0; for (size_t i = 0; i < k; i++) if (rec[3 * i + 1] == f) { nf++; if (rec[3 * i + 2] == 1) nf1++; }
    long long minsl = 1LL << 60; int neq = 0, nviol = 0, nsmall = 0;
    for (size_t i = 0; i < k; i++) {
        u64 h = rec[3 * i + 1], hc = rec[3 * i + 2];
        if (h > 153) continue;    /* violators/equality need hc <= 15 and H <= 153 (H >= n hc) */
        nsmall++;
        long long s = (long long)h - (long long)f - 2 * ((long long)hc - 1);
        if (s < minsl) minsl = s; if (s == 0) neq++; if (s < 0) nviol++;
    }
    printf("n=10 (ears of strong 9-classes with H <= 205, and of all uniquely Hamiltonian 9-classes): %zu classes kept\n", k);
    printf("   f(10) = %llu, floor classes %d (hc = 1: %d); ceiling on the %d strong classes with H <= 153: min slack %lld, equality %d, violations %d\n",
           f, nf, nf1, nsmall, minsl, neq, nviol);
    u64 band = ~0ULL; int c2 = 0;
    for (size_t i = 0; i < k; i++) if (rec[3 * i + 2] >= 2 && rec[3 * i + 1] < band) band = rec[3 * i + 1];
    for (size_t i = 0; i < k; i++) if (rec[3 * i + 2] == 2 && rec[3 * i + 1] == band) c2++;
    u64 f2 = ~0ULL; for (size_t i = 0; i < k; i++) if (rec[3 * i + 2] == 2 && rec[3 * i + 1] < f2) f2 = rec[3 * i + 1];
    printf("   min{H : hc >= 2} = %llu; f_2(10) = %llu (%d classes with hc = 2 at the minimum)\n", band, f2, c2);
    int nu = 0, ncut = 0, nsole = 0; u64 umax = 0; char uh[2000] = {0}; int nd = 0;
    char st[2000] = {0};
    for (size_t i = 0; i < k; i++) {
        u64 h = rec[3 * i + 1], hc = rec[3 * i + 2];
        if (h < 2000) st[h] = 1;
        if (hc != 1) continue;
        nu++; if (h > umax) umax = h; uh[h] = 1;
        uint16_t t[MAXN]; decode(10, rec[3 * i], t);
        ncut += one_arc_cut(10, t);
        int mo = 10, mi = 10; for (int v = 0; v < 10; v++) { int s = __builtin_popcount(t[v]); if (s < mo) mo = s; if (9 - s < mi) mi = 9 - s; }
        if (mo == 1 || mi == 1) nsole++;
    }
    for (int v = 0; v < 2000; v++) nd += uh[v];
    printf("   U(10): %d classes, %d distinct H, range [%llu, %llu]; one-arc cut in %d; vertex of out- or in-degree 1 in %d\n", nu, nd, f, umax, ncut, nsole);
    printf("   odd values in [125, %llu) missed by strong 10-classes:", band); int bh = 0;
    for (u64 v = 125; v < band; v += 2) if (!st[v]) { printf(" %llu", v); bh++; } printf(" (%d)\n", bh);
    printf("   odd values in [125, 205] missed by strong 10-classes:"); bh = 0;
    for (u64 v = 125; v <= 205; v += 2) if (!st[v]) { bh++; } printf(" %d\n", bh);
    free(rec);
}

#ifndef NO_MAIN
int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 8;
    int do_ears = argc > 2 && !strcmp(argv[2], "ears");
    static const int fv[] = {0, 1, 0, 3, 5, 9, 15, 25, 45, 75, 125, 225};
    for (int n = 1; n <= NMAX; n++) {
        double t0 = omp_get_wtime();
        build(n);
        census_order(n, n + 1 <= 11 ? fv[n + 1] : 0);
        if (n >= 3) printf("   [%.1f s]\n", omp_get_wtime() - t0);
    }
    if (do_ears && NMAX == 9) ears10();
    return 0;
}
#endif
