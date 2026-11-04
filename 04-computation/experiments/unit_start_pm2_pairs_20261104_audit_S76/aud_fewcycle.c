/*
 * aud_fewcycle.c  (independent audit of THM-4872 (B), opus-2026-11-04-S76 audit)
 *
 * Written from scratch.  Generates every strong tournament with at most K = 5 Hamiltonian cycles, order by order:
 *   - Moon: every strong tournament of order n >= 4 has a vertex v with T - v strong;
 *   - insertion: if T - v has a Hamiltonian cycle C and v has both in- and out-neighbours, then some consecutive x, y on
 *     C have x -> v -> y, and C' = C with v inserted is a Hamiltonian cycle of T from which C is recovered; so
 *     hc(T) >= hc(T - v).
 * Hence every strong T with hc <= 5 extends a strong parent with hc <= 5 by one vertex, and the level-by-level extension
 * of all classes is complete.  For a parent P the cycles of the extension by v with in-set X are exactly
 *     hc = sum_{y not in X, x in X} M[y][x],   M[y][x] = # Hamiltonian paths of P from y to x,
 * evaluated for all X by a Gray code.  M is obtained by enumerating every Hamiltonian path of P (DFS with two sound
 * prunings: every unvisited vertex keeps an in-neighbour among the unvisited vertices and the current one; at most one
 * unvisited vertex lacks an unvisited out-neighbour).
 *
 * Canonical form: own individualisation-refinement (ordered equitable partition by out-neighbour counts per cell,
 * branching on every vertex of the first non-singleton cell, no pruning); the canonical code is the lexicographically
 * least relabelled adjacency over all leaves, and |Aut| = number of leaves attaining it.  |Aut| is cross-checked by a
 * separate backtracking automorphism counter on every sole-exit-free class with hc in {3,4,5}, and hc of every such class
 * is recomputed by a subset DP (paths from vertex 0 returning to 0).
 *
 * A sole exit is a vertex of out-degree 1; "sole-exit-free" = no such vertex (THM-4869 notation).
 * Build: gcc -O3 -march=native -fopenmp -o aud_fewcycle aud_fewcycle.c ; run: ./aud_fewcycle NMAX
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define MAXN 16
#define K 5
typedef struct { uint16_t row[MAXN]; uint8_t hc; uint8_t aut; } Cls;

/* ---------- canonical form by individualisation-refinement ---------- */
typedef struct { int n; const uint16_t *out; uint16_t best[MAXN]; int have; int nbest; } CanonCtx;

static void refine(int n, const uint16_t *out, int *cell, int *ncell) {
    for (;;) {
        int k = *ncell;
        uint16_t cm[MAXN]; memset(cm, 0, sizeof cm);
        for (int v = 0; v < n; v++) cm[cell[v]] |= (uint16_t)(1u << v);
        /* key[v] = (cell, counts per cell) */
        int key[MAXN][MAXN + 1];
        for (int v = 0; v < n; v++) { key[v][0] = cell[v]; for (int c = 0; c < k; c++) key[v][c + 1] = __builtin_popcount(out[v] & cm[c]); }
        int ord[MAXN]; for (int v = 0; v < n; v++) ord[v] = v;
        /* insertion sort by key */
        for (int i = 1; i < n; i++) { int t = ord[i], j = i - 1;
            while (j >= 0 && memcmp(key[ord[j]], key[t], sizeof(int) * (k + 1)) > 0) { ord[j + 1] = ord[j]; j--; }
            ord[j + 1] = t; }
        int nc = 0, newcell[MAXN];
        for (int i = 0; i < n; i++) { if (i > 0 && memcmp(key[ord[i]], key[ord[i - 1]], sizeof(int) * (k + 1)) != 0) nc++; newcell[ord[i]] = nc; }
        nc++;
        memcpy(cell, newcell, sizeof(int) * n);
        if (nc == k) { *ncell = nc; return; }
        *ncell = nc;
    }
}

static void canon_rec(CanonCtx *C, int *cell, int ncell) {
    int n = C->n;
    refine(n, C->out, cell, &ncell);
    if (ncell == n) {
        uint16_t rl[MAXN]; memset(rl, 0, sizeof rl);
        for (int v = 0; v < n; v++) { uint16_t r = 0; uint16_t o = C->out[v]; while (o) { int w = __builtin_ctz(o); o &= o - 1; r |= (uint16_t)(1u << cell[w]); } rl[cell[v]] = r; }
        int cmp = C->have ? memcmp(rl, C->best, sizeof(uint16_t) * n) : -1;
        /* memcmp on little-endian uint16 is a fixed total order; any fixed total order works */
        if (cmp < 0) { memcpy(C->best, rl, sizeof(uint16_t) * n); C->have = 1; C->nbest = 1; }
        else if (cmp == 0) C->nbest++;
        return;
    }
    /* first non-singleton cell */
    int size[MAXN] = {0}; for (int v = 0; v < n; v++) size[cell[v]]++;
    int tc = 0; while (size[tc] == 1) tc++;
    for (int v = 0; v < n; v++) if (cell[v] == tc) {
        int c2[MAXN];
        for (int w = 0; w < n; w++) c2[w] = cell[w] > tc ? cell[w] + 1 : cell[w];
        for (int w = 0; w < n; w++) if (cell[w] == tc && w != v) c2[w] = tc + 1;
        c2[v] = tc;
        canon_rec(C, c2, ncell + 1);
    }
}

static int canon(int n, const uint16_t *out, uint16_t *res) {
    CanonCtx C; C.n = n; C.out = out; C.have = 0; C.nbest = 0;
    int cell[MAXN]; for (int v = 0; v < n; v++) cell[v] = 0;
    canon_rec(&C, cell, 1);
    memset(res, 0, sizeof(uint16_t) * MAXN);
    memcpy(res, C.best, sizeof(uint16_t) * n);
    return C.nbest;
}

/* ---------- independent automorphism counter (plain backtracking) ---------- */
static int aut_rec(int n, const uint16_t *out, int *img, uint16_t used, int v) {
    if (v == n) return 1;
    int cnt = 0, dv = __builtin_popcount(out[v]);
    for (int w = 0; w < n; w++) {
        if ((used >> w) & 1u) continue;
        if (__builtin_popcount(out[w]) != dv) continue;
        int ok = 1;
        for (int u = 0; u < v && ok; u++) {
            int a = (out[u] >> v) & 1, b = (out[img[u]] >> w) & 1;
            if (a != b) ok = 0;
        }
        if (!ok) continue;
        img[v] = w; cnt += aut_rec(n, out, img, (uint16_t)(used | (1u << w)), v + 1);
    }
    return cnt;
}
static int aut_count(int n, const uint16_t *out) { int img[MAXN]; return aut_rec(n, out, img, 0, 0); }

/* ---------- Hamiltonian cycles by DP (independent recount) ---------- */
static long long hc_dp(int n, const uint16_t *out, long long *dp) {
    unsigned full = (1u << n) - 1;
    memset(dp, 0, sizeof(long long) * ((size_t)1 << n) * n);
    dp[1 * n + 0] = 1;
    for (unsigned S = 1; S <= full; S += 2) for (int x = 0; x < n; x++) {
        long long c = dp[(size_t)S * n + x]; if (!c) continue;
        unsigned nb = out[x] & ~S;
        while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1u << y)) * n + y] += c; }
    }
    long long h = 0; for (int x = 1; x < n; x++) if ((out[x] >> 0) & 1u) h += dp[(size_t)full * n + x];
    return h;
}

/* ---------- all Hamiltonian paths of P: M[y][x] ---------- */
typedef struct { int n; const uint16_t *out; uint16_t in[MAXN]; long long M[MAXN][MAXN]; long long npaths; } PathCtx;
static void hp_rec(PathCtx *P, int start, int cur, uint16_t unv) {
    if (!unv) { P->M[start][cur]++; P->npaths++; return; }
    /* prune: every unvisited u needs an in-neighbour in unv | {cur}; at most one unvisited u lacks an out-neighbour in unv */
    uint16_t avail = (uint16_t)(unv | (1u << cur));
    int noout = 0; uint16_t t = unv;
    while (t) { int u = __builtin_ctz(t); t &= t - 1;
        if (!(P->in[u] & avail)) return;
        if (!(P->out[u] & unv)) { if (++noout > 1) return; } }
    uint16_t nb = P->out[cur] & unv;
    while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; hp_rec(P, start, y, (uint16_t)(unv & ~(1u << y))); }
}

static int cmpcls(const void *a, const void *b) { return memcmp(((const Cls *)a)->row, ((const Cls *)b)->row, sizeof(uint16_t) * MAXN); }

static long long fib(int k) { long long a = 0, b = 1; for (int i = 0; i < k; i++) { long long t = a + b; a = b; b = t; } return a; }

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 11;
    int VERIFY_ALL_HC = argc > 2 ? atoi(argv[2]) : 11;   /* recount hc of every class by DP up to this order */
    if (NMAX > 15) NMAX = 15;
    printf("# aud_fewcycle NMAX=%d (strong tournaments with hc <= %d by vertex extension; own canonical form)\n", NMAX, K);
    Cls *L = malloc(sizeof(Cls)); size_t nL = 1;
    memset(L, 0, sizeof(Cls));
    L[0].row[0] = 1u << 1; L[0].row[1] = 1u << 2; L[0].row[2] = 1u << 0; L[0].hc = 1; L[0].aut = 3;
    { uint16_t c[MAXN]; L[0].aut = (uint8_t)canon(3, L[0].row, c); memcpy(L[0].row, c, sizeof c); }
    for (int n = 3; n <= NMAX; n++) {
        if (n > 3) {
            int m = n - 1;    /* parent order */
            int nth = omp_get_max_threads();
            Cls **buf = calloc(nth, sizeof(Cls *)); size_t *bn = calloc(nth, sizeof(size_t)), *bc = calloc(nth, sizeof(size_t));
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                #pragma omp for schedule(dynamic, 8)
                for (size_t pi = 0; pi < nL; pi++) {
                    const uint16_t *po = L[pi].row;
                    PathCtx P; memset(&P, 0, sizeof P); P.n = m; P.out = po;
                    for (int x = 0; x < m; x++) { P.in[x] = 0; for (int y = 0; y < m; y++) if (y != x && ((po[y] >> x) & 1u)) P.in[x] |= (uint16_t)(1u << y); }
                    uint16_t fullm = (uint16_t)((1u << m) - 1);
                    for (int s = 0; s < m; s++) hp_rec(&P, s, s, (uint16_t)(fullm & ~(1u << s)));
                    /* Gray code over X (in-set of the new vertex) */
                    uint16_t X = 0; long long hcv = 0;
                    for (uint32_t g = 1; g < (1u << m); g++) {
                        int t = __builtin_ctz(g);
                        if (!((X >> t) & 1u)) {           /* add t to X */
                            long long add = 0, rem = 0;
                            for (int y = 0; y < m; y++) if (y != t && !((X >> y) & 1u)) add += P.M[y][t];
                            for (int x = 0; x < m; x++) if ((X >> x) & 1u) rem += P.M[t][x];
                            hcv += add - rem; X |= (uint16_t)(1u << t);
                        } else {                           /* remove t from X */
                            X &= (uint16_t)~(1u << t);
                            long long add = 0, rem = 0;
                            for (int y = 0; y < m; y++) if (y != t && !((X >> y) & 1u)) rem += P.M[y][t];
                            for (int x = 0; x < m; x++) if ((X >> x) & 1u) add += P.M[t][x];
                            hcv += add - rem;
                        }
                        if (X == fullm) continue;          /* v would have no out-neighbour */
                        if (hcv < 1 || hcv > K) continue;
                        uint16_t T[MAXN]; memset(T, 0, sizeof T);
                        for (int x = 0; x < m; x++) T[x] = po[x];
                        int v = m;
                        for (int x = 0; x < m; x++) { if ((X >> x) & 1u) T[x] |= (uint16_t)(1u << v); else T[v] |= (uint16_t)(1u << x); }
                        Cls c; memset(&c, 0, sizeof c);
                        c.aut = (uint8_t)canon(n, T, c.row); c.hc = (uint8_t)hcv;
                        if (bn[tid] == bc[tid]) { bc[tid] = bc[tid] ? 2 * bc[tid] : 4096; buf[tid] = realloc(buf[tid], sizeof(Cls) * bc[tid]); }
                        buf[tid][bn[tid]++] = c;
                    }
                }
            }
            size_t tot = 0; for (int t = 0; t < nth; t++) tot += bn[t];
            Cls *all = malloc(sizeof(Cls) * (tot ? tot : 1)); size_t k = 0;
            for (int t = 0; t < nth; t++) { memcpy(all + k, buf[t], sizeof(Cls) * bn[t]); k += bn[t]; free(buf[t]); }
            free(buf); free(bn); free(bc);
            qsort(all, tot, sizeof(Cls), cmpcls);
            size_t u = 0; long incons = 0;
            for (size_t i = 0; i < tot; i++) {
                if (u > 0 && memcmp(all[i].row, all[u - 1].row, sizeof(uint16_t) * MAXN) == 0) {
                    if (all[i].hc != all[u - 1].hc || all[i].aut != all[u - 1].aut) incons++;
                    continue;
                }
                all[u++] = all[i];
            }
            free(L); L = all; nL = u;
            printf("n=%2d extensions with hc<=5: %zu -> classes %zu (inconsistent duplicates %ld)\n", n, tot, nL, incons);
        } else printf("n=%2d classes 1 (C3)\n", n);
        /* statistics */
        long t[K + 1] = {0}, uu[K + 1] = {0}, autdist[K + 1][8], autall[K + 1][8];
        memset(autdist, 0, sizeof autdist); memset(autall, 0, sizeof autall);
        long hcbad = 0, autbad = 0, hcchecked = 0;
        long long *dp = malloc(sizeof(long long) * ((size_t)1 << n) * n);
        for (size_t i = 0; i < nL; i++) {
            Cls *c = &L[i]; int k2 = c->hc;
            t[k2]++;
            int a = c->aut; int ai = a == 1 ? 0 : a == 3 ? 1 : a == 5 ? 2 : a == 7 ? 3 : a == 9 ? 4 : a == 21 ? 5 : 6;
            autall[k2][ai]++;
            int sef = 1; for (int x = 0; x < n; x++) if (__builtin_popcount(c->row[x]) == 1) { sef = 0; break; }
            if (n <= VERIFY_ALL_HC || (sef && k2 >= 3)) { hcchecked++; if (hc_dp(n, c->row, dp) != k2) hcbad++; }
            if (sef) {
                uu[k2]++; autdist[k2][ai]++;
                if (k2 >= 3 && aut_count(n, c->row) != a) autbad++;
            }
        }
        free(dp);
        printf("   classes %zu; t_k (k=1..5): %ld %ld %ld %ld %ld; hc recounted by DP on %ld classes: %ld mismatches\n", nL, t[1], t[2], t[3], t[4], t[5], hcchecked, hcbad);
        printf("   sole-exit-free u_k (k=1..5): %ld %ld %ld %ld %ld   [F(2n-9) = %lld, 2F(2n-11) = %lld, 4F(2n-12)-3F(2n-14) = %lld]\n",
               uu[1], uu[2], uu[3], uu[4], uu[5], n >= 5 ? fib(2 * n - 9) : 0, n >= 6 ? 2 * fib(2 * n - 11) : 0, n >= 7 ? 4 * fib(2 * n - 12) - 3 * fib(2 * n - 14) : 0);
        for (int k2 = 3; k2 <= 5; k2++)
            printf("   u_%d by |Aut| (1,3,5,7,9,21,other): %ld %ld %ld %ld %ld %ld %ld%s\n", k2, autdist[k2][0], autdist[k2][1], autdist[k2][2], autdist[k2][3],
                   autdist[k2][4], autdist[k2][5], autdist[k2][6],
                   k2 == 3 && n >= 6 ? "" : "");
        if (n >= 6) printf("   predicted u_3 split: Z_3 F(2n-10) = %lld, asymmetric F(2n-11) = %lld\n", fib(2 * n - 10), fib(2 * n - 11));
        printf("   all classes by hc and |Aut| (1,3,5,7,9,21,other):");
        for (int k2 = 1; k2 <= 5; k2++) printf("  k=%d: %ld %ld %ld %ld %ld %ld %ld", k2, autall[k2][0], autall[k2][1], autall[k2][2], autall[k2][3], autall[k2][4], autall[k2][5], autall[k2][6]);
        printf("\n   independent backtracking |Aut| on the sole-exit-free hc = 3,4,5 classes: %ld mismatches\n", autbad);
        fflush(stdout);
    }
    return 0;
}
