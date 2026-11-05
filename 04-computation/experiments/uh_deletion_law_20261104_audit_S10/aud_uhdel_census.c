/* aud_uhdel_census.c -- independent audit code (uh_deletion_law_20261104_audit_S10). Shares no code with the session's
 * programs or with the HYP-9464 auditor's aud_fewcycle.c (different canonical form: individualisation-refinement
 * instead of minimum over Hamiltonian cycles; hc of the finished classes recomputed by an independent cycle DP).
 *
 * Mode "ext NMAX": every strong tournament class with 1 <= hc <= 5 at orders 3..NMAX by vertex extension.
 *   Completeness: a strong T of order n >= 4 has a vertex v with T - v strong (Moon 1966), and hc(T - v) <= hc(T)
 *   (insert v into a cycle of T - v at a switch x -> v -> y; deleting v recovers the cycle, so the map is injective).
 *   hc of an extension by z: cycles through z <-> parent Hamiltonian paths from N+(z) to N-(z).
 * Mode "lab NMAX" (NMAX <= 7): the same class counts from ALL labelled tournaments (hc by direct DP), classes counted by
 *   canonical form; used to validate the extension census and the canonical form at small orders.
 * For every sole-exit-free class (minimum out-degree >= 2) with hc in {3,4,5} it prints a DUMP line
 *   DUMP mode n hc aut row_0 ... row_{n-1} | hc(T - 0) ... hc(T - (n-1))
 * (rows = out-neighbour bitmasks in canonical labelling) and the joint counts of (k1, k2, |Aut|),
 * k1 = #{v : hc(T - v) = 1}, k2 = #{v : hc(T - v) = 2}. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define MAXN 12
#define KMAX 5
typedef unsigned long long u64;
typedef struct { uint16_t r[MAXN]; } Key;   /* canonical rows; unused rows = 0 */

static int popc(unsigned x) { return __builtin_popcount(x); }

/* ---------- canonical form by individualisation-refinement ---------- */
typedef struct { int n; const uint16_t *out; Key best; int have; int cnt; } CanonState;

/* refine colours col[] (values 0..ncol-1, label-invariant) to the coarsest equitable refinement; returns ncol */
static int refine(int n, const uint16_t *out, int *col) {
    int ncol = 0; for (int v = 0; v < n; v++) if (col[v] + 1 > ncol) ncol = col[v] + 1;
    for (;;) {
        /* signature of v: (col[v], number of out-neighbours in each colour) */
        int sig[MAXN][MAXN + 1];
        for (int v = 0; v < n; v++) {
            sig[v][0] = col[v];
            for (int c = 0; c < ncol; c++) sig[v][c + 1] = 0;
            for (int w = 0; w < n; w++) if ((out[v] >> w) & 1) sig[v][col[w] + 1]++;
        }
        /* rank distinct signatures in lexicographic order */
        int idx[MAXN]; for (int v = 0; v < n; v++) idx[v] = v;
        for (int a = 1; a < n; a++) { int x = idx[a], b = a - 1;
            while (b >= 0 && memcmp(sig[idx[b]], sig[x], sizeof(int) * (ncol + 1)) > 0) { idx[b + 1] = idx[b]; b--; }
            idx[b + 1] = x; }
        int nc[MAXN], r = 0;
        for (int a = 0; a < n; a++) {
            if (a > 0 && memcmp(sig[idx[a]], sig[idx[a - 1]], sizeof(int) * (ncol + 1))) r++;
            nc[idx[a]] = r;
        }
        int newn = r + 1;
        memcpy(col, nc, sizeof(int) * n);
        if (newn == ncol) return ncol;
        ncol = newn;
    }
}
static void search(CanonState *S, int *col) {
    int n = S->n;
    int ncol = refine(n, S->out, col);
    if (ncol == n) {
        Key k; memset(&k, 0, sizeof k);
        for (int v = 0; v < n; v++) { uint16_t row = 0; for (int w = 0; w < n; w++) if ((S->out[v] >> w) & 1) row |= (uint16_t)(1u << col[w]); k.r[col[v]] = row; }
        int c = S->have ? memcmp(&k, &S->best, sizeof k) : -1;
        if (c < 0) { S->best = k; S->have = 1; S->cnt = 1; } else if (c == 0) S->cnt++;
        return;
    }
    /* target cell: smallest colour with >= 2 vertices */
    int size[MAXN] = {0}; for (int v = 0; v < n; v++) size[col[v]]++;
    int tc = 0; while (size[tc] < 2) tc++;
    for (int x = 0; x < n; x++) {
        if (col[x] != tc) continue;
        int c2[MAXN]; for (int v = 0; v < n; v++) c2[v] = 2 * col[v] + ((v == x) ? 0 : 1);
        /* compress */
        int used[2 * MAXN] = {0}; for (int v = 0; v < n; v++) used[c2[v]] = 1;
        int map[2 * MAXN], r = 0; for (int c = 0; c < 2 * n; c++) if (used[c]) map[c] = r++;
        for (int v = 0; v < n; v++) c2[v] = map[c2[v]];
        search(S, c2);
    }
}
static Key canon(int n, const uint16_t *out, int *aut) {
    CanonState S; S.n = n; S.out = out; S.have = 0; S.cnt = 0;
    int col[MAXN]; for (int v = 0; v < n; v++) col[v] = popc(out[v]);
    /* compress scores to ranks */
    int used[MAXN + 1] = {0}; for (int v = 0; v < n; v++) used[col[v]] = 1;
    int map[MAXN + 1], r = 0; for (int c = 0; c <= n; c++) if (used[c]) map[c] = r++;
    for (int v = 0; v < n; v++) col[v] = map[col[v]];
    search(&S, col);
    if (aut) *aut = S.cnt;
    return S.best;
}

/* ---------- Hamiltonian counts ---------- */
/* P[s][t] = Hamiltonian paths from s to t (s != t), n >= 2 */
static void path_matrix(int n, const uint16_t *out, u64 P[MAXN][MAXN], u64 *dp) {
    unsigned full = (1u << n) - 1;
    for (int s = 0; s < n; s++) {
        memset(dp, 0, sizeof(u64) * ((size_t)full + 1) * n);
        dp[(size_t)(1u << s) * n + s] = 1;
        for (unsigned S = 1; S <= full; S++) {
            if (!((S >> s) & 1)) continue;
            for (int v = 0; v < n; v++) {
                u64 c = dp[(size_t)S * n + v]; if (!c) continue;
                unsigned nb = out[v] & ~S & full;
                while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1u << w)) * n + w] += c; }
            }
        }
        for (int t = 0; t < n; t++) P[s][t] = (t == s) ? 0 : dp[(size_t)full * n + t];
    }
}
/* number of directed Hamiltonian cycles: DP over paths starting at the LAST vertex n-1 (different anchor than usual) */
static u64 hcycles(int n, const uint16_t *out, u64 *dp) {
    if (n < 3) return 0;
    unsigned full = (1u << n) - 1; int a = n - 1;
    memset(dp, 0, sizeof(u64) * ((size_t)full + 1) * n);
    dp[(size_t)(1u << a) * n + a] = 1;
    for (unsigned S = 1; S <= full; S++) {
        if (!((S >> a) & 1)) continue;
        for (int v = 0; v < n; v++) {
            u64 c = dp[(size_t)S * n + v]; if (!c) continue;
            unsigned nb = out[v] & ~S & full;
            while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1u << w)) * n + w] += c; }
        }
    }
    u64 h = 0; for (int v = 0; v < n; v++) if (v != a && ((out[v] >> a) & 1)) h += dp[(size_t)full * n + v];
    return h;
}
static void delete_vertex(int n, const uint16_t *out, int x, uint16_t *res) {
    int k = 0;
    for (int i = 0; i < n; i++) {
        if (i == x) continue;
        uint16_t r = 0; int kk = 0;
        for (int j = 0; j < n; j++) { if (j == x) continue; if ((out[i] >> j) & 1) r |= (uint16_t)(1u << kk); kk++; }
        res[k++] = r;
    }
}
static int cmpkey(const void *a, const void *b) { return memcmp(a, b, sizeof(Key)); }

/* ---------- statistics over a list of classes ---------- */
static void stats(int n, Key *cls, long ncls, const char *tag) {
    u64 *dp = malloc(sizeof(u64) * ((size_t)1 << n) * n);
    long t[KMAX + 1] = {0}, u[KMAX + 1] = {0};
    static long joint[KMAX + 1][MAXN + 1][MAXN + 1][8];
    memset(joint, 0, sizeof joint);
    long autother = 0;
    for (long c = 0; c < ncls; c++) {
        const uint16_t *o = cls[c].r;
        u64 h = hcycles(n, o, dp);
        if (h < 1 || h > KMAX) { printf("ERROR: class with hc %llu\n", h); exit(1); }
        t[h]++;
        int mo = 99; for (int v = 0; v < n; v++) if (popc(o[v]) < mo) mo = popc(o[v]);
        if (mo < 2) continue;
        u[h]++;
        if (h < 3) continue;
        int aut; Key kk = canon(n, o, &aut);
        if (memcmp(&kk, &cls[c], sizeof(Key))) { printf("ERROR: stored form not canonical\n"); exit(1); }
        int k1 = 0, k2 = 0; u64 hv[MAXN];
        for (int v = 0; v < n; v++) {
            uint16_t r[MAXN]; delete_vertex(n, o, v, r);
            hv[v] = hcycles(n - 1, r, dp);
            if (hv[v] == 1) k1++;
            if (hv[v] == 2) k2++;
        }
        int ai = aut < 7 ? aut : 7; if (aut >= 7) autother++;
        joint[h][k1][k2][ai]++;
        printf("DUMP %s %d %llu %d", tag, n, h, aut);
        for (int v = 0; v < n; v++) printf(" %u", (unsigned)o[v]);
        printf(" |");
        for (int v = 0; v < n; v++) printf(" %llu", hv[v]);
        printf("\n");
    }
    printf("%s order %d: %ld strong classes with 1 <= hc <= %d; t_k:", tag, n, ncls, KMAX);
    for (int k = 1; k <= KMAX; k++) printf(" %ld", t[k]);
    printf("; u_k (min out-degree >= 2):");
    for (int k = 1; k <= KMAX; k++) printf(" %ld", u[k]);
    printf("\n");
    for (int h = 3; h <= KMAX; h++) {
        printf("%s JOINT n=%d hc=%d (k1,k2,|Aut|):", tag, n, h);
        for (int a = 0; a <= MAXN; a++) for (int b = 0; b <= MAXN; b++) for (int c = 0; c < 8; c++) if (joint[h][a][b][c])
            printf(" (%d,%d,%d%s):%ld", a, b, c, c == 7 ? "+" : "", joint[h][a][b][c]);
        printf("\n");
    }
    fflush(stdout);
    free(dp);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s ext|lab NMAX [threads]\n", argv[0]); return 1; }
    int NMAX = atoi(argv[2]);
    if (argc > 3) omp_set_num_threads(atoi(argv[3]));
    if (!strcmp(argv[1], "lab")) {
        if (NMAX > 7) NMAX = 7;
        for (int n = 3; n <= NMAX; n++) {
            int np = n * (n - 1) / 2, pi[32], pj[32], c = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { pi[c] = i; pj[c] = j; c++; }
            long cap = 1 << 16, cnt = 0; Key *keys = malloc(sizeof(Key) * cap);
            u64 *dp = malloc(sizeof(u64) * ((size_t)1 << n) * n);
            for (u64 bits = 0; bits < (1ull << np); bits++) {
                uint16_t out[MAXN]; memset(out, 0, sizeof out);
                for (int b = 0; b < np; b++) { if ((bits >> b) & 1) out[pi[b]] |= (uint16_t)(1u << pj[b]); else out[pj[b]] |= (uint16_t)(1u << pi[b]); }
                u64 h = hcycles(n, out, dp);
                if (h < 1 || h > KMAX) continue;
                if (cnt == cap) { cap *= 2; keys = realloc(keys, sizeof(Key) * cap); }
                keys[cnt++] = canon(n, out, NULL);
            }
            qsort(keys, cnt, sizeof(Key), cmpkey);
            long u = 0; for (long i = 0; i < cnt; i++) if (i == 0 || memcmp(&keys[i], &keys[i - 1], sizeof(Key))) keys[u++] = keys[i];
            printf("lab order %d: %ld labelled tournaments with 1 <= hc <= %d\n", n, cnt, KMAX);
            stats(n, keys, u, "lab");
            free(keys); free(dp);
        }
        return 0;
    }
    /* extension census */
    Key *cur = calloc(1, sizeof(Key)); cur[0].r[0] = 2; cur[0].r[1] = 4; cur[0].r[2] = 1; long ncur = 1;
    { int a; Key k = canon(3, cur[0].r, &a); cur[0] = k; }
    stats(3, cur, ncur, "ext");
    for (int n = 4; n <= NMAX; n++) {
        int m = n - 1; double t0 = omp_get_wtime();
        long cap = 1 << 20, cnt = 0; Key *cand = malloc(sizeof(Key) * cap);
        #pragma omp parallel
        {
            u64 *dp = malloc(sizeof(u64) * ((size_t)1 << m) * m);
            long lcap = 1 << 14, lcnt = 0; Key *loc = malloc(sizeof(Key) * lcap);
            #pragma omp for schedule(dynamic, 4)
            for (long c = 0; c < ncur; c++) {
                const uint16_t *pr = cur[c].r; u64 P[MAXN][MAXN];
                path_matrix(m, pr, P, dp);
                for (unsigned mask = 1; mask + 1 < (1u << m); mask++) {   /* mask = N+(z); 0 and full give no cycle */
                    u64 h = 0;
                    for (int a = 0; a < m && h <= KMAX; a++) if ((mask >> a) & 1) for (int b = 0; b < m; b++) if (!((mask >> b) & 1)) h += P[a][b];
                    if (h < 1 || h > KMAX) continue;
                    uint16_t t[MAXN]; memset(t, 0, sizeof t);
                    for (int i = 0; i < m; i++) t[i] = pr[i] | (((mask >> i) & 1) ? 0 : (uint16_t)(1u << m));
                    t[m] = (uint16_t)mask;
                    if (lcnt == lcap) { lcap *= 2; loc = realloc(loc, sizeof(Key) * lcap); }
                    loc[lcnt++] = canon(n, t, NULL);
                }
            }
            #pragma omp critical
            {
                if (cnt + lcnt > cap) { while (cnt + lcnt > cap) cap *= 2; cand = realloc(cand, sizeof(Key) * cap); }
                memcpy(cand + cnt, loc, sizeof(Key) * lcnt); cnt += lcnt;
            }
            free(loc); free(dp);
        }
        qsort(cand, cnt, sizeof(Key), cmpkey);
        long u = 0; for (long i = 0; i < cnt; i++) if (i == 0 || memcmp(&cand[i], &cand[i - 1], sizeof(Key))) cand[u++] = cand[i];
        free(cur); cur = cand; ncur = u;
        printf("ext order %d: %ld candidates, %ld classes [%.1f s]\n", n, cnt, ncur, omp_get_wtime() - t0);
        stats(n, cur, ncur, "ext");
    }
    return 0;
}
