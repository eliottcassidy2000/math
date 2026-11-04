/*
 * aud_census.c  (independent audit of THM-4869 / HYP-9461; opus-2026-11-04-S75 audit)
 *
 * Two generation modes, both with the audit's own canonical form (aud_common.h):
 *   exh NEX          : EXHAUSTIVE enumeration of all tournament classes of orders 1..NEX by adding a vertex in all 2^m
 *                      ways to every class of order m (trivially complete; no Moon theorem, no hc monotonicity).
 *                      Writes the strong classes with hc <= 40 of order NEX to aud_strong_hcle40_<NEX>.bin.
 *   ext FILE K NMAX  : Moon extension from the strong classes with hc <= K in FILE (order m0) to orders m0+1..NMAX.
 *                      Completeness: Moon 1966 gives v with T - v strong, and hc(T - v) <= hc(T).
 *                      Child hc = sum over w in O, u notin O of P(w -> u), computed by a lowest-bit recursion over O,
 *                      and re-checked by a direct DP on every kept child and on a sample of rejected children.
 * Per order, on the strong classes with hc <= KS (KS = 40 in exh mode, K in ext mode):
 *   - t_k, u_k (u: no vertex of out-degree 1), for k <= min(KS, 40);
 *   - 2-strong classes and 2-arc-strong classes with a cut vertex, by hc;
 *   - Theorem A: for every 2-arc-strong class with a cut vertex v (every such v), hc = S+ * prod H(D_i) * S-, S+, S- >= 3;
 *   - forced arcs (arcs on every Hamiltonian cycle) of 2-arc-strong classes (exh mode: ALL strong classes of the order);
 *   - for sole-exit-free classes with hc <= 5: lambda, all one-arc-leaving sets, minimal ones, uniqueness, product law,
 *     s_l = 1, prefix Q = T/D uniquely Hamiltonian with z a sole exit, a sole entry in T; core and prefix codes;
 *   - rooted pairs (Q, z), Q uniquely Hamiltonian, z a sole exit: R(n);
 *   - direct core enumeration (every strong D with hc <= 5 and every root u), c_k(n), and blocks (W, A) with hc(W) <= 2.
 */
#include "aud_common.h"
#include <omp.h>

#define KT 40
typedef struct { u128 code; int hc; } Rec;

/* -------------------------------------------------------------- per-arc cycle counts */
/* returns hc; cnt[a*n+b] = number of Hamiltonian cycles using arc a->b. F, G: (1<<n)*n i64 each */
static i64 arc_cycle_counts(int n, const uint32_t *out, i64 *cnt, i64 *F, i64 *G) {
    uint32_t full = (1u << n) - 1;
    /* F[S][v]: paths from 0 with vertex set S (0 in S) ending at v */
    memset(F, 0, sizeof(i64) * (size_t)(full + 1) * n);
    F[(size_t)1 * n + 0] = 1;
    for (uint32_t S = 1; S <= full; S += 2) {
        for (uint32_t r = S; r; r &= r - 1) {
            int v = __builtin_ctz(r); i64 c = F[(size_t)S * n + v]; if (!c) continue;
            uint32_t nb = out[v] & ~S;
            while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; F[(size_t)(S | (1u << w)) * n + w] += c; }
        }
    }
    /* G[S][v]: paths with vertex set S (0 not in S) starting at v and ending at some w with w -> 0 */
    memset(G, 0, sizeof(i64) * (size_t)(full + 1) * n);
    for (int v = 1; v < n; v++) if (out[v] & 1u) G[(size_t)(1u << v) * n + v] = 1;
    for (uint32_t S = 2; S <= full; S += 2) {
        if (__builtin_popcount(S) < 2) continue;
        for (uint32_t r = S; r; r &= r - 1) {
            int v = __builtin_ctz(r); i64 s = 0;
            uint32_t rest = S & ~(1u << v);
            for (uint32_t q = out[v] & rest; q; q &= q - 1) { int w = __builtin_ctz(q); s += G[(size_t)rest * n + w]; }
            G[(size_t)S * n + v] = s;
        }
    }
    i64 hc = 0;
    for (int v = 1; v < n; v++) if (out[v] & 1u) hc += F[(size_t)full * n + v];
    for (int i = 0; i < n * n; i++) cnt[i] = 0;
    for (int a = 0; a < n; a++) for (uint32_t r = out[a]; r; r &= r - 1) {
        int b = __builtin_ctz(r);
        if (a == 0) cnt[a * n + b] = G[(size_t)(full & ~1u) * n + b];
        else if (b == 0) cnt[a * n + b] = F[(size_t)full * n + a];
        else {
            i64 s = 0;
            /* S contains 0 and a, not b */
            uint32_t freebits = full & ~(1u | (1u << a) | (1u << b));
            for (uint32_t T = freebits;; T = (T - 1) & freebits) {
                uint32_t S = T | 1u | (1u << a);
                i64 f = F[(size_t)S * n + a];
                if (f) s += f * G[(size_t)(full & ~S) * n + b];
                if (!T) break;
            }
            cnt[a * n + b] = s;
        }
    }
    return hc;
}

/* strong components of T[V], in order (all arcs forward). Returns count; comp[i] = vertex mask */
static int components(int n, const uint32_t *out, uint32_t V, uint32_t *comp) {
    uint32_t reach[MAXV];
    for (uint32_t r = V; r; r &= r - 1) {
        int x = __builtin_ctz(r);
        uint32_t seen = 1u << x, fr = seen;
        while (fr) { int y = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = out[y] & V & ~seen; seen |= nb; fr |= nb; }
        reach[x] = seen;
    }
    int nc = 0; uint32_t left = V;
    while (left) {
        /* first component: vertices whose reach is all of left (they reach everything remaining) */
        uint32_t c = 0;
        for (uint32_t r = left; r; r &= r - 1) { int x = __builtin_ctz(r); if ((reach[x] & left) == left) c |= 1u << x; }
        if (!c) { fprintf(stderr, "components: no source component\n"); exit(1); }
        comp[nc++] = c; left &= ~c;
    }
    (void)n;
    return nc;
}

/* -------------------------------------------------------------- statistics */
typedef struct {
    u64 t[KT + 1], u[KT + 1], s2[KT + 1], cv[KT + 1], big;
    u64 thmA_checked, thmA_fail_formula, thmA_fail_bound, thmA_min_hc, thmA_minS;
    u64 forced_classes, forced_arcs, forced_no_cutvertex, forced_type[5];
    u64 dec_classes, dec_lambda2, dec_nmin_ne1, dec_small, dec_notstrong, dec_prodfail, dec_sl_ne1, dec_q_not_uh, dec_z_not_se, dec_no_entry;
    u64 dec_by_kd[6][MAXV + 1];
    u64 uh, uh_rooted, dec_overflow;
} Stats;

typedef struct { int k, d; u128 core; } CoreRec;
typedef struct { int m; u128 pre; } PreRec;

static void stats_add(Stats *a, const Stats *b) {
    for (int k = 0; k <= KT; k++) { a->t[k] += b->t[k]; a->u[k] += b->u[k]; a->s2[k] += b->s2[k]; a->cv[k] += b->cv[k]; }
    a->big += b->big;
    a->thmA_checked += b->thmA_checked; a->thmA_fail_formula += b->thmA_fail_formula; a->thmA_fail_bound += b->thmA_fail_bound;
    if (b->thmA_min_hc && (!a->thmA_min_hc || b->thmA_min_hc < a->thmA_min_hc)) a->thmA_min_hc = b->thmA_min_hc;
    if (b->thmA_minS && (!a->thmA_minS || b->thmA_minS < a->thmA_minS)) a->thmA_minS = b->thmA_minS;
    a->forced_classes += b->forced_classes; a->forced_arcs += b->forced_arcs; a->forced_no_cutvertex += b->forced_no_cutvertex;
    for (int i = 0; i < 5; i++) a->forced_type[i] += b->forced_type[i];
    a->dec_classes += b->dec_classes; a->dec_lambda2 += b->dec_lambda2; a->dec_nmin_ne1 += b->dec_nmin_ne1; a->dec_small += b->dec_small;
    a->dec_notstrong += b->dec_notstrong; a->dec_prodfail += b->dec_prodfail; a->dec_sl_ne1 += b->dec_sl_ne1;
    a->dec_q_not_uh += b->dec_q_not_uh; a->dec_z_not_se += b->dec_z_not_se; a->dec_no_entry += b->dec_no_entry;
    for (int k = 0; k < 6; k++) for (int d = 0; d <= MAXV; d++) a->dec_by_kd[k][d] += b->dec_by_kd[k][d];
    a->uh += b->uh; a->uh_rooted += b->uh_rooted; a->dec_overflow += b->dec_overflow;
}

/* sub-tournament on mask V relabelled 0..m-1 */
static int subtour(const uint32_t *out, uint32_t V, uint32_t *o2, int *vs) {
    int m = 0; for (uint32_t r = V; r; r &= r - 1) vs[m++] = __builtin_ctz(r);
    for (int i = 0; i < m; i++) { o2[i] = 0; for (int j = 0; j < m; j++) if ((out[vs[i]] >> vs[j]) & 1u) o2[i] |= 1u << j; }
    return m;
}

/* print one class compactly */
static void print_rows(int n, const uint32_t *out) {
    for (int v = 0; v < n; v++) { printf(" "); for (int w = 0; w < n; w++) printf("%d", (out[v] >> w) & 1); }
}

/* per-class analysis; hc known; KS = bound for "few-cycle" stats. forced: compute forced arcs for 2-arc-strong. */
static void analyze_class(int n, const uint32_t *out, i64 hc, int forced, Stats *S, i64 *buf, i64 *buf2,
                          CoreRec **cores, size_t *ncores, size_t *capcores, PreRec **pres, size_t *npres, size_t *cappres,
                          int printforced) {
    uint32_t full = (1u << n) - 1;
    int se = 0, sen = 0;
    for (int v = 0; v < n; v++) { if (__builtin_popcount(out[v]) == 1) se++; if (n - 1 - __builtin_popcount(out[v]) == 1) sen++; }
    int s2 = a_two_strong(n, out);
    int lam2 = s2 ? 1 : a_lambda2(n, out);
    int cutv = !s2;
    if (hc <= KT) {
        S->t[hc]++; if (!se) S->u[hc]++;
        if (s2) S->s2[hc]++; else if (lam2) S->cv[hc]++;
    } else S->big++;
    /* Theorem A on 2-arc-strong classes with a cut vertex (any hc) */
    if (lam2 && cutv) {
        for (int v = 0; v < n; v++) {
            uint32_t V = full & ~(1u << v);
            if (a_strong(n, out, V)) continue;
            uint32_t comp[MAXV]; int m = components(n, out, V, comp);
            i64 e[MAXV], st[MAXV];
            /* S+ over D_1 */
            a_end_counts(out, comp[0], st, buf, 1);    /* converse end counts = start counts */
            i64 Sp = 0; for (uint32_t r = comp[0] & out[v]; r; r &= r - 1) Sp += st[__builtin_ctz(r)];
            a_end_counts(out, comp[m - 1], e, buf, 0);
            i64 Sm = 0; for (uint32_t r = comp[m - 1]; r; r &= r - 1) { int b = __builtin_ctz(r); if ((out[b] >> v) & 1u) Sm += e[b]; }
            i64 mid = 1;
            for (int i = 1; i < m - 1; i++) { a_end_counts(out, comp[i], e, buf, 0); i64 H = 0; for (uint32_t r = comp[i]; r; r &= r - 1) H += e[__builtin_ctz(r)]; mid *= H; }
            S->thmA_checked++;
            if (Sp * mid * Sm != hc) S->thmA_fail_formula++;
            if (Sp < 3 || Sm < 3) S->thmA_fail_bound++;
            u64 mn = (u64)(Sp < Sm ? Sp : Sm);
            if (!S->thmA_minS || mn < S->thmA_minS) S->thmA_minS = mn;
        }
        if (!S->thmA_min_hc || (u64)hc < S->thmA_min_hc) S->thmA_min_hc = hc;
    }
    /* forced arcs of 2-arc-strong classes */
    if (forced && lam2) {
        i64 cnt[MAXV * MAXV];
        i64 h2 = arc_cycle_counts(n, out, cnt, buf, buf2);
        if (h2 != hc) { fprintf(stderr, "hc mismatch in arc counts\n"); exit(1); }
        int nf = 0; int fa[64], fb[64];
        for (int a = 0; a < n; a++) for (uint32_t r = out[a]; r; r &= r - 1) { int b = __builtin_ctz(r); if (cnt[a * n + b] == hc && nf < 64) { fa[nf] = a; fb[nf] = b; nf++; } }
        if (nf) {
            S->forced_classes++; S->forced_arcs += nf;
            if (!cutv) S->forced_no_cutvertex++;
            /* classify relative to the first cut vertex */
            int v0 = -1; for (int v = 0; v < n && v0 < 0; v++) if (!a_strong(n, out, full & ~(1u << v))) v0 = v;
            if (v0 >= 0) {
                uint32_t comp[MAXV]; int m = components(n, out, full & ~(1u << v0), comp);
                for (int i = 0; i < nf; i++) {
                    int a = fa[i], b = fb[i], ia = -1, ib = -1;
                    for (int c = 0; c < m; c++) { if ((comp[c] >> a) & 1u) ia = c; if ((comp[c] >> b) & 1u) ib = c; }
                    int type;
                    if (a == v0 || b == v0) type = 4;                 /* incident to v: claimed impossible */
                    else if (ia == ib && ia == 0) type = 0;           /* inside D_1 */
                    else if (ia == ib && ia == m - 1) type = 1;       /* inside D_m */
                    else if (ib == ia + 1) type = 2;                  /* junction */
                    else type = 3;                                    /* inside a middle component, or a skip arc */
                    S->forced_type[type]++;
                }
            }
            if (printforced) {
                #pragma omp critical
                {
                    printf("    forced-arc class n=%d hc=%lld cutvertex=%d forced arcs:", n, (long long)hc, cutv);
                    for (int i = 0; i < nf; i++) printf(" %d->%d", fa[i], fb[i]);
                    printf(" | scores"); for (int v = 0; v < n; v++) printf(" %d", __builtin_popcount(out[v]));
                    printf(" | rows"); print_rows(n, out); printf("\n");
                }
            }
        }
    }
    /* rooted uniquely Hamiltonian pairs (Q, z) with z a sole exit */
    if (hc == 1) {
        S->uh++;
        for (int z = 0; z < n; z++) if (__builtin_popcount(out[z]) == 1) {
            int col[MAXV]; for (int v = 0; v < n; v++) col[v] = (v == z) ? 0 : 1;
            u128 c = a_canon_col(n, out, col);
            if (*npres == *cappres) { *cappres *= 2; *pres = realloc(*pres, sizeof(PreRec) * *cappres); }
            (*pres)[*npres].m = n; (*pres)[*npres].pre = c; (*npres)++;
        }
    }
    /* decomposition of sole-exit-free classes with hc <= 5 */
    if (!se && hc <= 5) {
        S->dec_classes++;
        if (lam2) { S->dec_lambda2++; return; }
        /* all one-arc-leaving sets via d+(X) recursion */
        static __thread unsigned char dplus[1 << 16];
        uint32_t list[4096]; int nl = 0;
        dplus[0] = 0;
        for (uint32_t X = 1; X < full; X++) {
            int x = __builtin_ctz(X); uint32_t Xp = X & (X - 1);
            int add = __builtin_popcount(out[x] & ~X & full);
            int sub = 0; for (uint32_t r = Xp; r; r &= r - 1) if ((out[__builtin_ctz(r)] >> x) & 1u) sub++;
            int d = dplus[Xp] + add - sub;
            dplus[X] = (unsigned char)d;
            if (d == 1) { if (nl < 4096) list[nl++] = X; else S->dec_overflow++; }
        }
        int nmin = 0; uint32_t D = 0;
        for (int i = 0; i < nl; i++) {
            int minimal = 1;
            for (int j = 0; j < nl && minimal; j++) if (j != i && (list[j] & list[i]) == list[j]) minimal = 0;
            if (minimal) { nmin++; D = list[i]; }
        }
        if (nmin != 1) { S->dec_nmin_ne1++; return; }
        int d = __builtin_popcount(D);
        if (d < 3) S->dec_small++;
        if (!a_strong(n, out, D)) S->dec_notstrong++;
        int uu = -1, ll = -1;
        for (uint32_t r = D; r; r &= r - 1) { int x = __builtin_ctz(r); uint32_t o = out[x] & ~D & full; if (o) { uu = x; ll = __builtin_ctz(o); } }
        i64 e[MAXV], st[MAXV];
        a_end_counts(out, D, e, buf, 0);
        i64 eu = e[uu];
        uint32_t R = full & ~D; i64 sl;
        if (__builtin_popcount(R) == 1) sl = 1;
        else { a_end_counts(out, R, st, buf, 1); sl = st[ll]; }
        if (sl * eu != hc) S->dec_prodfail++;
        if (sl != 1) S->dec_sl_ne1++;
        if (!sen) S->dec_no_entry++;
        if (eu <= 5 && d <= MAXV) S->dec_by_kd[eu][d]++;
        /* core code: T[D] rooted at u */
        {
            uint32_t o2[MAXV]; int vs[MAXV]; int m = subtour(out, D, o2, vs);
            int col[MAXV]; for (int i = 0; i < m; i++) col[i] = (vs[i] == uu) ? 0 : 1;
            u128 c = a_canon_col(m, o2, col);
            if (*ncores == *capcores) { *capcores *= 2; *cores = realloc(*cores, sizeof(CoreRec) * *capcores); }
            (*cores)[*ncores].k = (int)eu; (*cores)[*ncores].d = d; (*cores)[*ncores].core = c; (*ncores)++;
        }
        /* prefix Q = T/D with z -> l only */
        if (__builtin_popcount(R) >= 2) {
            uint32_t q[MAXV]; int vs[MAXV]; int m = subtour(out, R, q, vs);
            int z = m; q[z] = 0;
            for (int i = 0; i < m; i++) { if (vs[i] == ll) q[z] |= 1u << i; else q[i] |= 1u << z; }
            int mq = m + 1;
            if (a_hc_sub(q, (1u << mq) - 1, buf) != 1) S->dec_q_not_uh++;
            if (__builtin_popcount(q[z]) != 1) S->dec_z_not_se++;
        }
    }
}

/* direct core enumeration over strong D with hc <= 5 (each root u), and blocks over strong W with hc <= 2 */
typedef struct { int k, d; u128 code; int asz; } KeyRec;

static int cmp_key(const void *a, const void *b) {
    const KeyRec *x = a, *y = b;
    if (x->k != y->k) return x->k - y->k;
    if (x->d != y->d) return x->d - y->d;
    if (x->asz != y->asz) return x->asz - y->asz;
    if (x->code < y->code) return -1;
    if (x->code > y->code) return 1;
    return 0;
}

static void cores_and_blocks(int n, Rec *list, size_t nl, u64 *corecnt, u64 *blockcnt, int printblocks) {
    KeyRec *keys = malloc(sizeof(KeyRec) * 1024); size_t nk = 0, ck = 1024;
    KeyRec *bk = malloc(sizeof(KeyRec) * 1024); size_t nb = 0, cb = 1024;
    #pragma omp parallel
    {
        i64 *buf = malloc(sizeof(i64) * ((size_t)1 << n) * n);
        KeyRec *lk = malloc(sizeof(KeyRec) * 256); size_t lnk = 0, lck = 256;
        KeyRec *lb = malloc(sizeof(KeyRec) * 256); size_t lnb = 0, lcb = 256;
        #pragma omp for schedule(dynamic, 64)
        for (size_t i = 0; i < nl; i++) {
            if (list[i].hc > 5) continue;
            uint32_t out[MAXV]; a_decode(n, list[i].code, out);
            uint32_t full = (1u << n) - 1;
            i64 e[MAXV]; a_end_counts(out, full, e, buf, 0);
            /* cores */
            for (int u = 0; u < n; u++) {
                if (e[u] < 1 || e[u] > 5) continue;
                int ok = 1;
                for (int x = 0; x < n && ok; x++) if (x != u && __builtin_popcount(out[x]) == 1) ok = 0;
                uint32_t Du = full & ~(1u << u);
                for (uint32_t Z = Du; Z && ok; Z = (Z - 1) & Du) {
                    int arcs = 0;
                    for (uint32_t r = Z; r && arcs < 2; r &= r - 1) arcs += __builtin_popcount(out[__builtin_ctz(r)] & ~Z & full);
                    if (arcs == 1) ok = 0;
                }
                if (!ok) continue;
                int col[MAXV]; for (int v = 0; v < n; v++) col[v] = (v == u) ? 0 : 1;
                u128 c = a_canon_col(n, out, col);
                if (lnk == lck) { lck *= 2; lk = realloc(lk, sizeof(KeyRec) * lck); }
                lk[lnk].k = (int)e[u]; lk[lnk].d = n; lk[lnk].code = c; lk[lnk].asz = 1; lnk++;
            }
            /* blocks (W, A): hc(W) <= 2 */
            if (list[i].hc <= 2) {
                for (uint32_t A = 1; A <= full; A++) {
                    if (__builtin_popcount(A) < 2) continue;
                    i64 s = 0; for (uint32_t r = A; r; r &= r - 1) s += e[__builtin_ctz(r)];
                    if (s > 5) continue;
                    int ok = 1;
                    for (int x = 0; x < n && ok; x++) if (!((A >> x) & 1u) && __builtin_popcount(out[x]) < 2) ok = 0;
                    for (uint32_t Z = full; Z && ok; Z = (Z - 1) & full) {
                        int arcs = __builtin_popcount(Z & A);
                        for (uint32_t r = Z; r && arcs < 2; r &= r - 1) arcs += __builtin_popcount(out[__builtin_ctz(r)] & ~Z & full);
                        if (arcs == 1) ok = 0;
                    }
                    if (!ok) continue;
                    int col[MAXV]; for (int v = 0; v < n; v++) col[v] = ((A >> v) & 1u) ? 0 : 1;
                    u128 c = a_canon_col(n, out, col);
                    if (lnb == lcb) { lcb *= 2; lb = realloc(lb, sizeof(KeyRec) * lcb); }
                    lb[lnb].k = (int)s; lb[lnb].d = n; lb[lnb].code = c; lb[lnb].asz = __builtin_popcount(A); lnb++;
                    if (printblocks) {
                        #pragma omp critical
                        {
                            printf("    block n=%d k=%lld hc(W)=%d |A|=%d end counts", n, (long long)s, list[i].hc, __builtin_popcount(A));
                            for (int v = 0; v < n; v++) printf(" %lld", (long long)e[v]);
                            printf(" A:"); for (int v = 0; v < n; v++) if ((A >> v) & 1u) printf(" %d", v);
                            printf(" rows"); print_rows(n, out); printf("\n");
                        }
                    }
                }
            }
        }
        #pragma omp critical
        {
            for (size_t j = 0; j < lnk; j++) { if (nk == ck) { ck *= 2; keys = realloc(keys, sizeof(KeyRec) * ck); } keys[nk++] = lk[j]; }
            for (size_t j = 0; j < lnb; j++) { if (nb == cb) { cb *= 2; bk = realloc(bk, sizeof(KeyRec) * cb); } bk[nb++] = lb[j]; }
        }
        free(lk); free(lb); free(buf);
    }
    qsort(keys, nk, sizeof(KeyRec), cmp_key);
    for (int k = 0; k <= 5; k++) corecnt[k] = 0;
    for (size_t j = 0; j < nk; j++) if (j == 0 || cmp_key(&keys[j], &keys[j - 1]) != 0) corecnt[keys[j].k]++;
    qsort(bk, nb, sizeof(KeyRec), cmp_key);
    for (int k = 0; k <= 5; k++) blockcnt[k] = 0;
    for (size_t j = 0; j < nb; j++) if (j == 0 || cmp_key(&bk[j], &bk[j - 1]) != 0) blockcnt[bk[j].k]++;
    free(keys); free(bk);
}

static u64 core_formula(int k, int d) {
    if (k == 3) return d == 5 ? 2 : (d >= 6 ? 3ull << (d - 6) : 0);
    if (k == 4) return d == 5 ? 2 : (d >= 6 ? 1ull << (d - 5) : 0);
    if (k == 5) return d == 6 ? 4 : (d >= 7 ? 5ull << (d - 7) : 0);
    return 0;
}

static int cmp_u128(const void *a, const void *b) { u128 x = *(const u128 *)a, y = *(const u128 *)b; return x < y ? -1 : x > y; }

/* global records across orders */
static u64 Robs[MAXV + 2];          /* rooted (UH, sole exit) pairs by order, observed */
static u64 core_seen[6][MAXV + 2];  /* distinct cores (k, d) seen as terminal cores of u-classes */

static void report_order(int n, Rec *list, size_t nl, int KS, int forced, int printforced, int do_cores) {
    Stats tot; memset(&tot, 0, sizeof tot);
    CoreRec *allc = malloc(sizeof(CoreRec) * 1024); size_t nall = 0, call = 1024;
    PreRec *allp = malloc(sizeof(PreRec) * 1024); size_t nallp = 0, callp = 1024;
    #pragma omp parallel
    {
        Stats S; memset(&S, 0, sizeof S);
        i64 *buf = malloc(sizeof(i64) * ((size_t)1 << n) * n), *buf2 = malloc(sizeof(i64) * ((size_t)1 << n) * n);
        CoreRec *cores = malloc(sizeof(CoreRec) * 256); size_t nc = 0, cc = 256;
        PreRec *pres = malloc(sizeof(PreRec) * 256); size_t np = 0, cp = 256;
        #pragma omp for schedule(dynamic, 256)
        for (size_t i = 0; i < nl; i++) {
            if (list[i].hc < 1) continue;          /* not strong */
            uint32_t out[MAXV]; a_decode(n, list[i].code, out);
            analyze_class(n, out, list[i].hc, forced, &S, buf, buf2, &cores, &nc, &cc, &pres, &np, &cp, printforced);
        }
        #pragma omp critical
        {
            stats_add(&tot, &S);
            for (size_t j = 0; j < nc; j++) { if (nall == call) { call *= 2; allc = realloc(allc, sizeof(CoreRec) * call); } allc[nall++] = cores[j]; }
            for (size_t j = 0; j < np; j++) { if (nallp == callp) { callp *= 2; allp = realloc(allp, sizeof(PreRec) * callp); } allp[nallp++] = pres[j]; }
        }
        free(buf); free(buf2); free(cores); free(pres);
    }
    int kmax = KS < KT ? KS : KT;
    {
        size_t ns = 0, nk = 0;
        for (size_t i = 0; i < nl; i++) if (list[i].hc >= 1) { ns++; if (list[i].hc <= KS) nk++; }
        printf("n=%2d classes in list %zu: strong %zu, strong with hc<=%d: %zu\n", n, nl, ns, KS, nk);
    }
    printf("   t_k (k=1..%d):", kmax); for (int k = 1; k <= kmax; k++) printf(" %llu", tot.t[k]); printf("\n");
    printf("   u_k (k=1..%d):", kmax); for (int k = 1; k <= kmax; k++) printf(" %llu", tot.u[k]); printf("\n");
    printf("   2-strong by k:"); int first = 0; for (int k = 1; k <= kmax; k++) if (tot.s2[k]) { printf(" %d:%llu", k, tot.s2[k]); if (!first) first = k; }
    printf("   (min hc of 2-strong <= %d: %d)\n", kmax, first);
    printf("   2-arc-strong with cut vertex by k:"); first = 0; for (int k = 1; k <= kmax; k++) if (tot.cv[k]) { printf(" %d:%llu", k, tot.cv[k]); if (!first) first = k; }
    printf("   (min: %d)\n", first);
    printf("   Theorem A: (class, cut vertex) pairs checked %llu, formula failures %llu, S+ or S- < 3: %llu, min hc %llu, min min(S+,S-) %llu\n",
           tot.thmA_checked, tot.thmA_fail_formula, tot.thmA_fail_bound, tot.thmA_min_hc, tot.thmA_minS);
    if (forced) printf("   forced arcs among 2-arc-strong classes: %llu classes, %llu arcs, without cut vertex %llu; types (D1, Dm, junction, middle/skip, at v): %llu %llu %llu %llu %llu\n",
           tot.forced_classes, tot.forced_arcs, tot.forced_no_cutvertex, tot.forced_type[0], tot.forced_type[1], tot.forced_type[2], tot.forced_type[3], tot.forced_type[4]);
    printf("   sole-exit-free hc<=5 classes %llu: lambda>=2 %llu, #minimal!=1 %llu, |D|<3 %llu, D not strong %llu, product failures %llu, s_l!=1 %llu, Q not UH %llu, z not sole exit %llu, no sole entry %llu, list overflow %llu\n",
           tot.dec_classes, tot.dec_lambda2, tot.dec_nmin_ne1, tot.dec_small, tot.dec_notstrong, tot.dec_prodfail, tot.dec_sl_ne1, tot.dec_q_not_uh, tot.dec_z_not_se, tot.dec_no_entry, tot.dec_overflow);
    /* rooted prefixes R(n) */
    {
        u128 *codes = malloc(sizeof(u128) * (nallp + 1)); size_t m = 0;
        for (size_t j = 0; j < nallp; j++) codes[m++] = allp[j].pre;
        qsort(codes, m, sizeof(u128), cmp_u128);
        u64 r = 0; for (size_t j = 0; j < m; j++) if (j == 0 || codes[j] != codes[j - 1]) r++;
        Robs[n] = r; free(codes);
        printf("   uniquely Hamiltonian classes %llu (F(2n-6) = %llu); rooted (UH, sole exit) pairs R(%d) = %llu (F(2n-5) = %llu)\n",
               tot.uh, a_fib(2 * n - 6), n, r, a_fib(2 * n - 5));
    }
    /* terminal cores seen, and the class table against c_k(d) R(n - d + 1) */
    {
        for (int k = 3; k <= 5; k++) {
            /* distinct cores of each order among this order's u-classes */
            printf("   k=%d: classes by core order d:", k);
            int okall = 1;
            for (int d = 3; d <= n; d++) {
                u64 c = tot.dec_by_kd[k][d]; if (!c) continue;
                u64 Rm = (n - d + 1 == 2) ? 1 : Robs[n - d + 1];
                u64 pred = core_formula(k, d) * Rm;
                printf(" d=%d:%llu%s", d, c, c == pred ? "" : "(!)");
                if (c != pred) okall = 0;
            }
            printf("  -> %s c_k(d) R(n-d+1)\n", okall ? "all equal" : "MISMATCH with");
        }
        /* distinct cores by (k, d) */
        u128 *codes = malloc(sizeof(u128) * (nall + 1));
        for (int k = 3; k <= 5; k++) for (int d = 3; d <= n; d++) {
            size_t m = 0;
            for (size_t j = 0; j < nall; j++) if (allc[j].k == k && allc[j].d == d) codes[m++] = allc[j].core;
            if (!m) continue;
            qsort(codes, m, sizeof(u128), cmp_u128);
            u64 r = 0; for (size_t j = 0; j < m; j++) if (j == 0 || codes[j] != codes[j - 1]) r++;
            if (r > core_seen[k][d]) core_seen[k][d] = r;
        }
        free(codes);
        printf("   distinct terminal cores seen so far, c_k(d) for d = 5..%d:\n", n - 1);
        for (int k = 3; k <= 5; k++) {
            printf("     k=%d:", k);
            for (int d = 5; d <= n - 1; d++) printf(" %llu%s", core_seen[k][d], core_seen[k][d] == core_formula(k, d) ? "" : "(!)");
            printf("\n");
        }
    }
    if (do_cores) {
        u64 cc[6], bc[6];
        cores_and_blocks(n, list, nl, cc, bc, n <= 6);
        printf("   DIRECT core enumeration at order d=%d: c_2..c_5 = %llu %llu %llu %llu (formula %llu %llu %llu %llu); blocks (W,A) of order %d, k=2..5: %llu %llu %llu %llu\n",
               n, cc[2], cc[3], cc[4], cc[5], core_formula(2, n), core_formula(3, n), core_formula(4, n), core_formula(5, n), n, bc[2], bc[3], bc[4], bc[5]);
    }
    fflush(stdout);
    free(allc); free(allp);
}

/* -------------------------------------------------------------- generation */
static Rec *extend_all(int m, Rec *cur, size_t ncur, size_t *nout) {
    /* every class of order m+1 = class of order m plus a vertex with any out-set O */
    size_t est = ncur * 64 + 1024; size_t cap = 1; while (cap < est * 2) cap <<= 1;
    if (m + 1 == 10) cap = (size_t)1 << 24;
    AHash h; ah_init(&h, cap);
    Rec *nxt = malloc(sizeof(Rec) * 1024); size_t nn = 0, cn = 1024;
    size_t CH = 2048;
    for (size_t base = 0; base < ncur; base += CH) {
        size_t hi = base + CH < ncur ? base + CH : ncur;
        u128 *buf = malloc(sizeof(u128) * (hi - base) << m);
        #pragma omp parallel for schedule(dynamic, 16)
        for (size_t i = base; i < hi; i++) {
            uint32_t out[MAXV]; a_decode(m, cur[i].code, out);
            for (uint32_t O = 0; O < (1u << m); O++) {
                uint32_t o2[MAXV];
                for (int v = 0; v < m; v++) o2[v] = out[v] | (((O >> v) & 1u) ? 0u : (1u << m));
                o2[m] = O;
                buf[((i - base) << m) + O] = a_canon(m + 1, o2);
            }
        }
        size_t tot = (hi - base) << m;
        for (size_t j = 0; j < tot; j++) if (ah_insert(&h, buf[j])) {
            if (nn == cn) { cn *= 2; nxt = realloc(nxt, sizeof(Rec) * cn); }
            nxt[nn].code = buf[j]; nxt[nn].hc = -1; nn++;
        }
        free(buf);
    }
    ah_free(&h);
    *nout = nn; return nxt;
}

static void compute_hc(int n, Rec *list, size_t nl) {
    #pragma omp parallel
    {
        i64 *buf = malloc(sizeof(i64) * ((size_t)1 << n) * n);
        #pragma omp for schedule(dynamic, 1024)
        for (size_t i = 0; i < nl; i++) {
            uint32_t out[MAXV]; a_decode(n, list[i].code, out);
            uint32_t full = (1u << n) - 1;
            if (n < 3 || !a_strong(n, out, full)) { list[i].hc = 0; continue; }
            i64 h = a_hc_sub(out, full, buf);
            list[i].hc = h > 1000000 ? 1000000 : (int)h;
        }
        free(buf);
    }
}

static Rec *extend_moon(int m, Rec *cur, size_t ncur, int K, size_t *nout, u64 *rechecked, u64 *sampled) {
    size_t cap = 1; while (cap < ncur * 16 + 4096) cap <<= 1;
    AHash h; ah_init(&h, cap);
    Rec *nxt = malloc(sizeof(Rec) * 1024); size_t nn = 0, cn = 1024;
    u64 rc = 0, sc = 0;
    #pragma omp parallel
    {
        i64 *P = malloc(sizeof(i64) * MAXV * MAXV), *buf = malloc(sizeof(i64) * ((size_t)1 << (m + 1)) * (m + 1));
        i64 *f = malloc(sizeof(i64) * ((size_t)1 << m));
        Rec *loc = malloc(sizeof(Rec) * 1024); size_t nl = 0, cl = 1024;
        u64 lrc = 0, lsc = 0;
        #pragma omp for schedule(dynamic, 8)
        for (size_t i = 0; i < ncur; i++) {
            uint32_t out[MAXV]; a_decode(m, cur[i].code, out);
            a_all_pairs(m, out, P, buf);
            i64 start[MAXV]; for (int s = 0; s < m; s++) { start[s] = 0; for (int t = 0; t < m; t++) start[s] += P[s * m + t]; }
            uint32_t full = (1u << m) - 1;
            f[0] = 0;
            for (uint32_t O = 1; O <= full; O++) {
                int x = __builtin_ctz(O); uint32_t Op = O & (O - 1);
                i64 val = f[Op] + start[x];
                for (uint32_t r = Op; r; r &= r - 1) { int w = __builtin_ctz(r); val -= P[w * m + x] + P[x * m + w]; }
                f[O] = val;
                if (O == full) continue;
                int keep = (val >= 1 && val <= K);
                int sample = !keep && ((O * 2654435761u + (uint32_t)i) % 997u == 0u);
                if (!keep && !sample) continue;
                uint32_t o2[MAXV];
                for (int v = 0; v < m; v++) o2[v] = out[v] | (((O >> v) & 1u) ? 0u : (1u << m));
                o2[m] = O;
                i64 direct = a_hc_sub(o2, (1u << (m + 1)) - 1, buf);
                if (direct != val) {
                    #pragma omp critical
                    { fprintf(stderr, "child hc formula mismatch: %lld vs %lld\n", (long long)val, (long long)direct); exit(1); }
                }
                if (sample) { lsc++; continue; }
                lrc++;
                u128 c = a_canon(m + 1, o2);
                if (nl == cl) { cl *= 2; loc = realloc(loc, sizeof(Rec) * cl); }
                loc[nl].code = c; loc[nl].hc = (int)val; nl++;
            }
        }
        #pragma omp critical
        {
            for (size_t j = 0; j < nl; j++) if (ah_insert(&h, loc[j].code)) {
                if (nn == cn) { cn *= 2; nxt = realloc(nxt, sizeof(Rec) * cn); }
                nxt[nn++] = loc[j];
            }
            rc += lrc; sc += lsc;
        }
        free(P); free(buf); free(f); free(loc);
    }
    ah_free(&h);
    *rechecked = rc; *sampled = sc;
    *nout = nn; return nxt;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: aud_census exh NEX | ext FILE K NMAX\n"); return 1; }
    printf("# aud_census %s threads=%d\n", argv[1], omp_get_max_threads());
    if (!strcmp(argv[1], "exh")) {
        int NEX = atoi(argv[2]);
        Rec *cur = malloc(sizeof(Rec)); size_t ncur = 1; cur[0].code = 0; cur[0].hc = 0;
        static const u64 A000568[] = {0, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056, 903753248ull};
        static const u64 A051337[] = {0, 1, 0, 1, 1, 6, 35, 353, 6008, 178133, 9355949};
        for (int m = 1; m < NEX; m++) {
            size_t nn; Rec *nx = extend_all(m, cur, ncur, &nn);
            free(cur); cur = nx; ncur = nn;
            int n = m + 1;
            compute_hc(n, cur, ncur);
            size_t nstrong = 0; for (size_t i = 0; i < ncur; i++) if (cur[i].hc >= 1) nstrong++;
            printf("CLASSES n=%d: %zu (A000568 %llu), strong %zu (A051337 %llu)\n", n, ncur, A000568[n], nstrong, n <= 10 ? A051337[n] : 0ull);
            fflush(stdout);
            if (n >= 3) {
                /* stats on all strong classes (forced arcs for every 2-arc-strong class) */
                report_order(n, cur, ncur, KT, 1, n == 9 || n == 10, 0);
                /* direct cores/blocks on the strong hc <= 5 sublist */
                size_t k5 = 0; Rec *l5 = malloc(sizeof(Rec) * (ncur + 1));
                for (size_t i = 0; i < ncur; i++) if (cur[i].hc >= 1 && cur[i].hc <= 5) l5[k5++] = cur[i];
                u64 cc[6], bc[6];
                cores_and_blocks(n, l5, k5, cc, bc, n <= 6);
                printf("   DIRECT core enumeration at order d=%d: c_2..c_5 = %llu %llu %llu %llu (formula %llu %llu %llu %llu); blocks (W,A) of order %d, k=2..5: %llu %llu %llu %llu\n",
                       n, cc[2], cc[3], cc[4], cc[5], core_formula(2, n), core_formula(3, n), core_formula(4, n), core_formula(5, n), n, bc[2], bc[3], bc[4], bc[5]);
                free(l5);
                fflush(stdout);
            }
            if (n == NEX) {
                char fn[256]; sprintf(fn, "aud_strong_hcle40_%d.bin", n);
                FILE *fp = fopen(fn, "wb"); size_t w = 0;
                for (size_t i = 0; i < ncur; i++) if (cur[i].hc >= 1 && cur[i].hc <= KT) { fwrite(&cur[i], sizeof(Rec), 1, fp); w++; }
                fclose(fp);
                printf("wrote %zu strong classes with hc <= %d to %s\n", w, KT, fn);
            }
        }
        free(cur);
    } else if (!strcmp(argv[1], "ext")) {
        const char *fn = argv[2]; int K = atoi(argv[3]); int NMAX = atoi(argv[4]);
        int m0 = 0; { const char *p = strrchr(fn, '_'); m0 = atoi(p + 1); }
        FILE *fp = fopen(fn, "rb"); if (!fp) { fprintf(stderr, "no file\n"); return 1; }
        Rec *cur = malloc(sizeof(Rec) * 1024); size_t ncur = 0, cap = 1024; Rec r;
        while (fread(&r, sizeof(Rec), 1, fp) == 1) if (r.hc <= K) { if (ncur == cap) { cap *= 2; cur = realloc(cur, sizeof(Rec) * cap); } cur[ncur++] = r; }
        fclose(fp);
        printf("read %zu strong classes of order %d with hc <= %d\n", ncur, m0, K);
        /* R(m) for m < m0: the exh run (aud_census_exh_10.out) found R(m) = F(2m-5) at every m <= 10; reuse those values */
        for (int m = 3; m < m0; m++) Robs[m] = a_fib(2 * m - 5);
        /* R(m) for m <= m0 must be known for the class-table check: recompute from the start list (hc = 1 classes) */
        report_order(m0, cur, ncur, K, 0, 0, 1);
        for (int m = m0; m < NMAX; m++) {
            size_t nn; u64 rc, sc;
            Rec *nx = extend_moon(m, cur, ncur, K, &nn, &rc, &sc);
            free(cur); cur = nx; ncur = nn;
            printf("EXTENDED to n=%d: %zu classes with hc <= %d (children re-checked by direct DP %llu, rejected children sampled and re-checked %llu)\n", m + 1, ncur, K, rc, sc);
            report_order(m + 1, cur, ncur, K, 0, 0, K <= 5);
        }
        free(cur);
    }
    return 0;
}
