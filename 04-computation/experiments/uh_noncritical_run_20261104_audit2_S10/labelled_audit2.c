/* labelled_audit2.c -- independent audit (audit2, 2026-11-04) of THM-4868 (G) and (H).
 *
 * Canonical-form-free check on LABELLED tournaments.  Nothing here uses the author's programs,
 * the author's generator, or the theorem being audited: uniqueness of the Hamiltonian cycle is
 * decided by a direct DFS count, non-criticality by direct strong-connectivity tests, e_p by a
 * direct subset DP.
 *
 * Usage:
 *   labelled_audit2 A n   all 2^C(n,2) labelled tournaments on {0..n-1}       (n <= 7)
 *   labelled_audit2 B n   all tournaments containing the fixed Hamiltonian cycle
 *                         0->1->...->n-1->0 (the other C(n,2)-n arcs free)    (n <= 9)
 *
 * For every tournament T with exactly one Hamiltonian cycle (UH) it checks:
 *   (C-input) NC(T) is one cyclic run of the cycle, its first vertex y is a sole entry and the
 *             unique non-critical sole entry (these are only inputs for (G),(H); audited elsewhere);
 *   (G)       adm(T) = {p : e_p(T) = 1} equals {y, y-} plus y-- iff y- is a sole entry;
 *   (L3)      hc(T^p) = e_p(T) for EVERY vertex p (exact DP count on T^p, n+1 vertices);
 *   (H)       for p in adm(T): T^p is UH and r(T^p) - r(T) = 0, 1, 2 for p = y, y-, y--; the new
 *             vertex is the run start of T^p with predecessor p; plus the four facts used in the proof
 *             (p critical in T^p; x != p; x+ critical in T^p; y- non-critical in T^p when p = y--);
 *   labelled bijection counts:  mode A: L(n+1) = (n+1) * sum_U |adm(U)|;
 *                                mode B: n * L0(n+1) = (n+1) * sum_U |adm(U)|  (compare across runs).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define MAXV 12
typedef uint64_t u64;
typedef struct { int n; uint32_t out[MAXV], in[MAXV]; } Tour;

static void finish(Tour *T) {
    uint32_t full = (1u << T->n) - 1;
    for (int v = 0; v < T->n; v++) T->in[v] = full & ~T->out[v] & ~(1u << v);
}

/* DFS count of directed Hamiltonian cycles (each counted once, as a sequence starting at 0), capped */
typedef struct { const Tour *T; int cnt, cap; int path[MAXV]; int cyc[MAXV]; } Dfs;
static int dfs_rec(Dfs *c, int cur, uint32_t vis, int depth) {
    const Tour *T = c->T;
    if (depth == T->n) {
        if (T->out[cur] & 1u) { if (c->cnt == 0) memcpy(c->cyc, c->path, sizeof(int) * T->n); c->cnt++; }
        return c->cnt >= c->cap;
    }
    uint32_t cand = T->out[cur] & ~vis;
    while (cand) {
        int w = __builtin_ctz(cand); cand &= cand - 1;
        c->path[depth] = w;
        if (dfs_rec(c, w, vis | (1u << w), depth + 1)) return 1;
    }
    return 0;
}
static int hc_capped(const Tour *T, int cap, int *cyc) {
    Dfs c; c.T = T; c.cnt = 0; c.cap = cap; c.path[0] = 0;
    dfs_rec(&c, 0, 1u, 1);
    if (c.cnt > 0 && cyc) memcpy(cyc, c.cyc, sizeof(int) * T->n);
    return c.cnt;
}

/* exact Hamiltonian cycle count by DP from vertex 0 */
static u64 hc_exact(const Tour *T) {
    int n = T->n; size_t N = (size_t)1 << (n - 1);
    u64 *f = calloc(N * n, sizeof(u64));
    for (int v = 1; v < n; v++) if ((T->out[0] >> v) & 1) f[((size_t)1 << (v - 1)) * n + v] = 1;
    for (size_t S = 1; S < N; S++)
        for (int v = 1; v < n; v++) {
            if (!((S >> (v - 1)) & 1)) continue;
            u64 val = f[S * n + v]; if (!val) continue;
            uint32_t cand = (T->out[v] >> 1) & ~(uint32_t)S;
            while (cand) { int b = __builtin_ctz(cand); cand &= cand - 1; f[(S | ((size_t)1 << b)) * n + b + 1] += val; }
        }
    u64 tot = 0;
    for (int v = 1; v < n; v++) if (T->out[v] & 1u) tot += f[(N - 1) * n + v];
    free(f);
    return tot;
}

/* e[p] = number of Hamiltonian paths ending at p (any start) */
static void endcounts(const Tour *T, u64 *e) {
    int n = T->n; size_t N = (size_t)1 << n;
    u64 *f = calloc(N * n, sizeof(u64));
    for (int v = 0; v < n; v++) f[((size_t)1 << v) * n + v] = 1;
    for (size_t S = 1; S < N; S++)
        for (int v = 0; v < n; v++) {
            if (!((S >> v) & 1)) continue;
            u64 val = f[S * n + v]; if (!val) continue;
            uint32_t cand = T->out[v] & ~(uint32_t)S;
            while (cand) { int w = __builtin_ctz(cand); cand &= cand - 1; f[(S | ((size_t)1 << w)) * n + w] += val; }
        }
    for (int p = 0; p < n; p++) e[p] = f[(N - 1) * n + p];
    free(f);
}

static int strong_sub(const Tour *T, uint32_t M) {
    if (!M) return 1;
    int s = __builtin_ctz(M);
    uint32_t seen = 1u << s, fr = seen;
    while (fr) { int u = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = T->out[u] & M & ~seen; seen |= nb; fr |= nb; }
    if (seen != M) return 0;
    seen = 1u << s; fr = seen;
    while (fr) { int u = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = T->in[u] & M & ~seen; seen |= nb; fr |= nb; }
    return seen == M;
}

static uint32_t nc_mask(const Tour *T) {
    uint32_t full = (1u << T->n) - 1, m = 0;
    for (int v = 0; v < T->n; v++) if (strong_sub(T, full & ~(1u << v))) m |= 1u << v;
    return m;
}

/* child T^p: new vertex n, p -> n, n -> every other vertex */
static void child(const Tour *T, int p, Tour *C) {
    int n = T->n; C->n = n + 1;
    for (int v = 0; v < n; v++) C->out[v] = T->out[v];
    C->out[p] |= 1u << n;
    C->out[n] = ((1u << n) - 1) & ~(1u << p);
    finish(C);
}

typedef struct {
    long long uh, bad_cinput, bad_G, bad_L3, bad_H, bad_facts, bad_childstart;
    long long sum_adm, adm2, adm3, adm_other;
    long long Nr[MAXV + 2], N2r[MAXV + 2];
    long long delta_cnt[3];
} Stats;

static void stats_add(Stats *a, const Stats *b) {
    a->uh += b->uh; a->bad_cinput += b->bad_cinput; a->bad_G += b->bad_G; a->bad_L3 += b->bad_L3;
    a->bad_H += b->bad_H; a->bad_facts += b->bad_facts; a->bad_childstart += b->bad_childstart;
    a->sum_adm += b->sum_adm; a->adm2 += b->adm2; a->adm3 += b->adm3; a->adm_other += b->adm_other;
    for (int i = 0; i < MAXV + 2; i++) { a->Nr[i] += b->Nr[i]; a->N2r[i] += b->N2r[i]; }
    for (int i = 0; i < 3; i++) a->delta_cnt[i] += b->delta_cnt[i];
}

/* analysis of a UH tournament whose unique cycle is cyc[] */
typedef struct { int ok; int y, x, r; uint32_t NC; int succ[MAXV + 1], pred[MAXV + 1]; } Run;

static int analyse_run(const Tour *T, const int *cyc, Run *R) {
    int n = T->n;
    for (int i = 0; i < n; i++) { R->succ[cyc[i]] = cyc[(i + 1) % n]; R->pred[cyc[i]] = cyc[(i + n - 1) % n]; }
    R->NC = nc_mask(T);
    R->r = __builtin_popcount(R->NC);
    uint32_t full = (1u << n) - 1;
    if (R->NC == 0 || R->NC == full) return 0;
    int starts = 0, y = -1;
    for (int v = 0; v < n; v++) if (((R->NC >> v) & 1) && !((R->NC >> R->pred[v]) & 1)) { starts++; y = v; }
    if (starts != 1) return 0;
    int x = y; while ((R->NC >> R->succ[x]) & 1) x = R->succ[x];
    R->y = y; R->x = x;
    if (__builtin_popcount(T->in[y]) != 1) return 0;
    if (__builtin_popcount(T->out[x]) != 1) return 0;
    for (int v = 0; v < n; v++)  /* unique non-critical sole entry */
        if (v != y && ((R->NC >> v) & 1) && __builtin_popcount(T->in[v]) == 1) return 0;
    return 1;
}

static void process_uh(const Tour *T, const int *cyc, Stats *st) {
    int n = T->n;
    Run R;
    st->uh++;
    if (!analyse_run(T, cyc, &R)) { st->bad_cinput++; return; }
    int y = R.y, ym = R.pred[y], ymm = R.pred[ym];
    int ym_sole = (__builtin_popcount(T->in[ym]) == 1);
    int nsole = 0; for (int v = 0; v < n; v++) if (__builtin_popcount(T->in[v]) == 1) nsole++;
    st->Nr[R.r]++;
    if (nsole == 2) st->N2r[R.r]++;
    u64 e[MAXV];
    endcounts(T, e);
    uint32_t adm = 0;
    for (int p = 0; p < n; p++) if (e[p] == 1) adm |= 1u << p;
    uint32_t pred_adm = (1u << y) | (1u << ym) | (ym_sole ? (1u << ymm) : 0u);
    if (adm != pred_adm) st->bad_G++;
    int a = __builtin_popcount(adm);
    st->sum_adm += a;
    if (a == 2) st->adm2++; else if (a == 3) st->adm3++; else st->adm_other++;
    for (int p = 0; p < n; p++) {
        Tour C; child(T, p, &C);
        u64 h = hc_exact(&C);
        if (h != e[p]) st->bad_L3++;
        if (!((adm >> p) & 1)) continue;
        /* (H) */
        int ccyc[MAXV + 1];
        int hc = hc_capped(&C, 2, ccyc);
        if (hc != 1) { st->bad_H++; continue; }
        Run RC;
        if (!analyse_run(&C, ccyc, &RC)) { st->bad_H++; continue; }
        int expect = (p == y) ? 0 : (p == ym) ? 1 : 2;
        if (RC.r - R.r != expect) st->bad_H++; else st->delta_cnt[expect]++;
        if (RC.y != n || RC.pred[n] != p) st->bad_childstart++;
        /* facts used in the proof */
        int f_ok = 1;
        if ((RC.NC >> p) & 1) f_ok = 0;                       /* p critical in child */
        if (R.x == p) f_ok = 0;                                /* x_U != p */
        if ((RC.NC >> R.succ[R.x]) & 1) f_ok = 0;              /* x_U+ critical in child */
        if (__builtin_popcount(C.out[R.x]) != 1) f_ok = 0;     /* x_U still a sole exit */
        if (expect == 2 && !((RC.NC >> ym) & 1)) f_ok = 0;     /* y_U- non-critical in child */
        for (int v = 0; v < n; v++)                            /* NC(U) - p stays NC */
            if (v != p && ((R.NC >> v) & 1) && !((RC.NC >> v) & 1)) f_ok = 0;
        if (!f_ok) st->bad_facts++;
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s A|B n\n", argv[0]); return 1; }
    char mode = argv[1][0];
    int n = atoi(argv[2]);
    int pi[64], pj[64], M = 0;
    if (mode == 'A') {
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { pi[M] = i; pj[M] = j; M++; }
    } else {
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            if (j == i + 1) continue;
            if (i == 0 && j == n - 1) continue;
            pi[M] = i; pj[M] = j; M++;
        }
    }
    long long total = 1LL << M;
    Stats glob; memset(&glob, 0, sizeof glob);
    double t0 = omp_get_wtime();
#pragma omp parallel
    {
        Stats loc; memset(&loc, 0, sizeof loc);
#pragma omp for schedule(dynamic, 4096)
        for (long long code = 0; code < total; code++) {
            Tour T; T.n = n; memset(T.out, 0, sizeof T.out);
            if (mode == 'B') { for (int i = 0; i + 1 < n; i++) T.out[i] |= 1u << (i + 1); T.out[n - 1] |= 1u; }
            for (int b = 0; b < M; b++) {
                if ((code >> b) & 1) T.out[pi[b]] |= 1u << pj[b]; else T.out[pj[b]] |= 1u << pi[b];
            }
            finish(&T);
            int cyc[MAXV];
            int hc = hc_capped(&T, 2, cyc);
            if (hc != 1) continue;
            process_uh(&T, cyc, &loc);
        }
#pragma omp critical
        stats_add(&glob, &loc);
    }
    double t1 = omp_get_wtime();
    printf("mode %c n=%d tournaments=%lld  UH(labelled%s)=%lld\n", mode, n, total, mode == 'B' ? ", fixed cycle" : "", glob.uh);
    printf("  bad: C-input=%lld G=%lld L3=%lld H=%lld facts=%lld childstart=%lld\n",
           glob.bad_cinput, glob.bad_G, glob.bad_L3, glob.bad_H, glob.bad_facts, glob.bad_childstart);
    printf("  |adm|=2: %lld  |adm|=3: %lld  other: %lld  sum|adm| = %lld\n", glob.adm2, glob.adm3, glob.adm_other, glob.sum_adm);
    printf("  child run-length deltas 0/1/2: %lld %lld %lld\n", glob.delta_cnt[0], glob.delta_cnt[1], glob.delta_cnt[2]);
    printf("  labelled N(n,r) r=2..n-1:");
    for (int r = 2; r <= n - 1; r++) printf(" %lld", glob.Nr[r]);
    printf("\n  labelled N2(n,r) (two sole entries) r=2..n-1:");
    for (int r = 2; r <= n - 1; r++) printf(" %lld", glob.N2r[r]);
    printf("\n  time %.2fs\n", t1 - t0);
    return 0;
}
