/*
 * aud_uh.c -- independent audit code for THM-4867 / HYP-9458 (opus-2026-11-04-S74 audit).
 * Written from scratch; shares no code with uh_one_arc_cut_20261104.c.
 *
 * Three independent generators of U(n) = uniquely Hamiltonian tournaments up to isomorphism:
 *   (i)   BRUTE: every chord orientation of the fixed cycle 0->1->...->n-1->0 (no lemma used), n <= 9;
 *   (ii)  EXT:   every one-vertex extension of U(n-1) by an ARBITRARY out-set O (not only cyclic intervals),
 *                filtered by hc = 1 computed as sum_{w in O, u notin O} P(w->u); completeness needs only Moon's
 *                deletion theorem + "T - v strong => hc(T - v) <= hc(T)"; also records whether O is a cyclic interval;
 *   (iii) SHELOR: the five-letter grammar as recorded in session103 section 1 (seed S4, letters a..e).
 * Isomorphism: an isomorphism of uniquely Hamiltonian tournaments maps the unique cycle onto itself, so the minimum
 * over the n rotations of the chord code read along the cycle is a complete invariant (rotation-min code).
 * Hamiltonian cycle / path counts by subset DP saturated at 2 (bitmask "at least 1" / "at least 2").
 * Cut arcs: a cycle arc u -> s is a cut arc iff u cannot reach s in T - (u->s) (T - e keeps the path s ... u).
 *
 * Usage: aud_uh brute NMAX | aud_uh ext NMAX | aud_uh shelor NMAX_FULL NMAX_CUT | aud_uh all NB NE NS NC
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

typedef unsigned long long u64;
typedef unsigned __int128 u128;
typedef uint32_t msk;
#define MAXN 24

static int popc(msk x) { return __builtin_popcount(x); }

/* rotation-min code along a given Hamiltonian cycle cyc[0..n-1] (requires n <= 16 for 120 bits) */
static u128 rotcode(int n, const msk *out, const int *cyc) {
    u128 best = 0; int have = 0;
    for (int r = 0; r < n; r++) {
        u128 c = 0;
        for (int i = 0; i < n; i++) {
            int vi = cyc[(r + i) % n];
            for (int j = i + 1; j < n; j++) {
                int vj = cyc[(r + j) % n];
                c = (c << 1) | ((out[vi] >> vj) & 1u);
            }
        }
        if (!have || c < best) { best = c; have = 1; }
    }
    return best;
}
/* decode a rotation-min code to out[] with the cycle 0 -> 1 -> ... -> n-1 -> 0 */
static void decode(int n, u128 c, msk *out) {
    int tot = n * (n - 1) / 2, b = tot - 1;
    for (int i = 0; i < n; i++) out[i] = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        int bit = (int)((c >> b) & 1); b--;
        if (bit) out[i] |= 1u << j; else out[j] |= 1u << i;
    }
}

/* saturating Hamiltonian-cycle count (0, 1, 2 meaning >= 2) using vertex 0 as start */
static int hc_sat(int n, const msk *out, msk *d1, msk *d2) {
    msk in[MAXN];
    for (int v = 0; v < n; v++) { in[v] = 0; }
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if ((out[u] >> v) & 1) in[v] |= 1u << u;
    int full = (1 << n) - 1;
    memset(d1, 0, sizeof(msk) * (size_t)(full + 1));
    memset(d2, 0, sizeof(msk) * (size_t)(full + 1));
    d1[1] = 1;
    for (int S = 1; S <= full; S += 2) {
        if (!d1[S]) continue;
        msk e1 = d1[S], e2 = d2[S];
        msk rest = (msk)full & ~(msk)S;
        while (rest) {
            int w = __builtin_ctz(rest); rest &= rest - 1;
            msk c1 = e1 & in[w]; if (!c1) continue;
            msk c2 = e2 & in[w];
            int T = S | (1 << w);
            d1[T] |= 1u << w;
            if (c2 || (c1 & (c1 - 1))) d2[T] |= 1u << w;
        }
    }
    msk c1 = d1[full] & in[0], c2 = d2[full] & in[0];
    if (c2 || (c1 & (c1 - 1))) return 2;
    return c1 ? 1 : 0;
}

/* all-pairs saturating Hamiltonian path counts: P1[s] = {t : P(s->t) >= 1}, P2[s] = {t : >= 2} */
static void allpairs_sat(int m, const msk *out, msk *P1, msk *P2, msk *d1, msk *d2) {
    msk in[MAXN];
    for (int v = 0; v < m; v++) in[v] = 0;
    for (int u = 0; u < m; u++) for (int v = 0; v < m; v++) if ((out[u] >> v) & 1) in[v] |= 1u << u;
    int full = (1 << m) - 1;
    for (int s = 0; s < m; s++) {
        memset(d1, 0, sizeof(msk) * (size_t)(full + 1));
        memset(d2, 0, sizeof(msk) * (size_t)(full + 1));
        d1[1 << s] = 1u << s;
        for (int S = 1; S <= full; S++) {
            if (!((S >> s) & 1) || !d1[S]) continue;
            msk e1 = d1[S], e2 = d2[S];
            msk rest = (msk)full & ~(msk)S;
            while (rest) {
                int w = __builtin_ctz(rest); rest &= rest - 1;
                msk c1 = e1 & in[w]; if (!c1) continue;
                msk c2 = e2 & in[w];
                int T = S | (1 << w);
                d1[T] |= 1u << w;
                if (c2 || (c1 & (c1 - 1))) d2[T] |= 1u << w;
            }
        }
        P1[s] = d1[full]; P2[s] = d2[full];
    }
}

/* number of cut arcs on the cycle given by succ[] */
static int cut_arcs(int n, const msk *out, const int *succ) {
    int k = 0;
    for (int u = 0; u < n; u++) {
        int s = succ[u];
        msk reach = 1u << u, frontier = reach;
        while (frontier) {
            msk nx = 0, f = frontier;
            while (f) { int x = __builtin_ctz(f); f &= f - 1; msk o = out[x]; if (x == u) o &= ~(1u << s); nx |= o; }
            nx &= ~reach; reach |= nx; frontier = nx;
        }
        if (!((reach >> s) & 1)) k++;
    }
    return k;
}

static int cmp128(const void *a, const void *b) { u128 x = *(const u128 *)a, y = *(const u128 *)b; return x < y ? -1 : x > y; }
static size_t uniq128(u128 *a, size_t n) { if (!n) return 0; qsort(a, n, sizeof(u128), cmp128); size_t k = 1; for (size_t i = 1; i < n; i++) if (a[i] != a[k - 1]) a[k++] = a[i]; return k; }

static long long fibo(int k) { long long a = 0, b = 1; for (int i = 0; i < k; i++) { long long t = a + b; a = b; b = t; } return a; }

/* ------------------------------ (i) brute force over chord orientations ------------------------------ */
static u128 *brute(int n, size_t *cnt) {
    int ch_i[64], ch_j[64], nc = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        if (j == i + 1 || (i == 0 && j == n - 1)) continue;
        ch_i[nc] = i; ch_j[nc] = j; nc++;
    }
    long long total = 1LL << nc;
    int nth = omp_get_max_threads();
    u128 **res = calloc(nth, sizeof(u128 *)); size_t *rn = calloc(nth, sizeof(size_t)), *rc = calloc(nth, sizeof(size_t));
    long long survivors = 0;
#pragma omp parallel reduction(+:survivors)
    {
        int t = omp_get_thread_num();
        msk *d1 = malloc(sizeof(msk) << n), *d2 = malloc(sizeof(msk) << n);
        rc[t] = 1024; res[t] = malloc(sizeof(u128) * rc[t]);
        int cyc[MAXN]; for (int i = 0; i < n; i++) cyc[i] = i;
#pragma omp for schedule(dynamic, 4096)
        for (long long mask = 0; mask < total; mask++) {
            msk out[MAXN];
            for (int i = 0; i < n; i++) out[i] = 1u << ((i + 1) % n);
            for (int k = 0; k < nc; k++) { if ((mask >> k) & 1) out[ch_i[k]] |= 1u << ch_j[k]; else out[ch_j[k]] |= 1u << ch_i[k]; }
            /* quick necessary filter: relocating one vertex c_i into another slot gives a second cycle */
            int bad = 0;
            for (int i = 0; i < n && !bad; i++) {
                int pv = (i + n - 1) % n, nv = (i + 1) % n;
                if (!((out[pv] >> nv) & 1)) continue;          /* need c_(i-1) -> c_(i+1) to shortcut c_i */
                for (int j = 0; j < n; j++) {                   /* slot c_j -> c_(j+1), both != c_i, j != i-1 */
                    int a = j, b = (j + 1) % n;
                    if (a == i || b == i || a == pv) continue;
                    if (((out[a] >> i) & 1) && ((out[i] >> b) & 1)) { bad = 1; break; }
                }
            }
            if (bad) continue;
            survivors++;
            if (hc_sat(n, out, d1, d2) != 1) continue;
            if (rn[t] == rc[t]) { rc[t] *= 2; res[t] = realloc(res[t], sizeof(u128) * rc[t]); }
            res[t][rn[t]++] = rotcode(n, out, cyc);
        }
        free(d1); free(d2);
    }
    size_t tot = 0; for (int t = 0; t < nth; t++) tot += rn[t];
    u128 *all = malloc(sizeof(u128) * (tot + 1)); size_t p = 0;
    for (int t = 0; t < nth; t++) { memcpy(all + p, res[t], sizeof(u128) * rn[t]); p += rn[t]; free(res[t]); }
    size_t labelled = tot;
    *cnt = uniq128(all, tot);
    printf("BRUTE n=%d: chords %d, orientations %lld, filter survivors %lld, cycle-labelled uniquely Hamiltonian %zu, classes %zu (Garey F(2n-6) = %lld)\n",
           n, nc, total, survivors, labelled, *cnt, n >= 4 ? fibo(2 * n - 6) : 1);
    free(res); free(rn); free(rc);
    return all;
}

/* ------------------------------ (ii) extension by arbitrary out-sets ------------------------------ */
typedef struct { u128 *codes; size_t n; } Level;
static Level ext_level(int n, const Level *par, long long *e1, long long *e2, long long *nonint, long long *cand) {
    int m = n - 1;
    int nth = omp_get_max_threads();
    u128 **res = calloc(nth, sizeof(u128 *)); size_t *rn = calloc(nth, sizeof(size_t)), *rc = calloc(nth, sizeof(size_t));
    long long E1 = 0, E2 = 0, NI = 0, C = 0;
#pragma omp parallel reduction(+:E1, E2, NI, C)
    {
        int t = omp_get_thread_num();
        msk *d1 = malloc(sizeof(msk) << m), *d2 = malloc(sizeof(msk) << m);
        rc[t] = 1024; res[t] = malloc(sizeof(u128) * rc[t]);
#pragma omp for schedule(dynamic, 1)
        for (long r = 0; r < (long)par->n; r++) {
            msk tp[MAXN]; decode(m, par->codes[r], tp);       /* cycle 0 -> 1 -> ... -> m-1 -> 0 */
            msk P1[MAXN], P2[MAXN];
            allpairs_sat(m, tp, P1, P2, d1, d2);
            msk fullm = (1u << m) - 1;
            for (msk O = 1; O < fullm; O++) {
                C++;
                msk I = fullm & ~O;
                int hc = 0;
                msk o = O;
                while (o && hc < 2) {
                    int w = __builtin_ctz(o); o &= o - 1;
                    if (P2[w] & I) { hc = 2; break; }
                    hc += popc(P1[w] & I);
                }
                if (hc != 1) continue;
                /* switches along the parent cycle: positions i with c_i in I (c_i -> v) and c_(i+1) in O (v -> c_(i+1)) */
                int sw = 0, pos = -1;
                for (int i = 0; i < m; i++) { int a = i, b = (i + 1) % m; if (((I >> a) & 1) && ((O >> b) & 1)) { sw++; pos = i; } }
                if (sw != 1) { NI++; continue; }
                int po = popc(O);
                if (po == 1 || po == m - 1) E1++; else E2++;
                msk tt[MAXN];
                for (int i = 0; i < m; i++) tt[i] = tp[i] | (((I >> i) & 1) ? (1u << m) : 0);
                tt[m] = O;
                int cyc[MAXN], k = 0;
                for (int i = 0; i < m; i++) { cyc[k++] = i; if (i == pos) cyc[k++] = m; }
                u128 code = rotcode(n, tt, cyc);
                if (rn[t] == rc[t]) { rc[t] *= 2; res[t] = realloc(res[t], sizeof(u128) * rc[t]); }
                res[t][rn[t]++] = code;
            }
        }
        free(d1); free(d2);
    }
    size_t tot = 0; for (int t = 0; t < nth; t++) tot += rn[t];
    Level L; L.codes = malloc(sizeof(u128) * (tot + 1)); size_t p = 0;
    for (int t = 0; t < nth; t++) { memcpy(L.codes + p, res[t], sizeof(u128) * rn[t]); p += rn[t]; free(res[t]); }
    L.n = uniq128(L.codes, tot);
    free(res); free(rn); free(rc);
    *e1 = E1; *e2 = E2; *nonint = NI; *cand = C;
    return L;
}

/* ------------------------------ (iii) Shelor's five-letter grammar ------------------------------ */
typedef struct { int n; msk out[MAXN]; int succ[MAXN], pred[MAXN]; int type, a, a0, a1; int nU, nD; } GS;
static void gs_seed(GS *g) {
    memset(g, 0, sizeof *g);
    g->n = 4;
    g->out[0] = (1u << 1) | (1u << 2); g->out[1] = (1u << 2) | (1u << 3); g->out[2] = 1u << 3; g->out[3] = 1u << 0;
    for (int i = 0; i < 4; i++) { g->succ[i] = (i + 1) % 4; g->pred[i] = (i + 3) % 4; }
    g->type = 1; g->a0 = 0; g->a1 = 1;   /* type B: a0 -> a1 */
}
/* letter L in 'a'..'e'; returns 0 if illegal */
static int gs_apply(const GS *g, char L, GS *h) {
    int p, nt, na = -1, na0 = -1, na1 = -1, x = g->n, dU = 0, dD = 0;
    if (g->type == 0) {
        if (L == 'a') { p = g->pred[g->a]; nt = 0; }
        else if (L == 'b') { p = g->a; nt = 1; dU = 1; }
        else return 0;
    } else {
        if (L == 'c') { p = g->pred[g->a0]; nt = 0; dD = 1; }
        else if (L == 'd') { p = g->a0; nt = 1; }
        else if (L == 'e') { p = g->a1; nt = 1; dU = 1; }
        else return 0;
    }
    *h = *g;
    h->n = x + 1;
    h->out[p] |= 1u << x;
    h->out[x] = ((1u << x) - 1) & ~(1u << p);
    int s = g->succ[p];
    h->succ[x] = s; h->pred[s] = x; h->succ[p] = x; h->pred[x] = p;
    h->type = nt;
    if (nt == 0) { na = x; h->a = na; }
    else { na0 = p; na1 = x; h->a0 = na0; h->a1 = na1; }
    h->nU = g->nU + dU; h->nD = g->nD + dD;
    return 1;
}

typedef struct {
    long long words, hc_bad, cut_formula_bad, noin1, noout1, inserted_not_in1;
    long long hist[MAXN + 1];
    u128 *codes; size_t nc, cap;
} GAcc;

static void gs_visit(const GS *g, int target, int dofull, GAcc *A, msk *d1, msk *d2) {
    if (g->n == target) {
        int n = g->n;
        A->words++;
        int k = cut_arcs(n, g->out, g->succ);
        if (k != 3 + g->nU - g->nD) A->cut_formula_bad++;
        A->hist[k]++;
        int in1 = 0, out1 = 0;
        for (int v = 0; v < n; v++) { int od = popc(g->out[v]); if (od == 1) out1 = 1; if (n - 1 - od == 1) in1 = 1; }
        if (!in1) A->noin1++;
        if (!out1) A->noout1++;
        if (n > 4 && popc(g->out[n - 1]) != n - 2) A->inserted_not_in1++;
        if (dofull) {
            if (hc_sat(n, g->out, d1, d2) != 1) A->hc_bad++;
            int cyc[MAXN], v = 0;
            for (int i = 0; i < n; i++) { cyc[i] = v; v = g->succ[v]; }
            if (A->nc == A->cap) { A->cap = A->cap ? 2 * A->cap : 1024; A->codes = realloc(A->codes, sizeof(u128) * A->cap); }
            A->codes[A->nc++] = rotcode(n, g->out, cyc);
        }
        return;
    }
    const char *Ls = "abcde";
    for (int i = 0; i < 5; i++) { GS h; if (gs_apply(g, Ls[i], &h)) gs_visit(&h, target, dofull, A, d1, d2); }
}

static double binom(int n, int k) { if (k < 0 || k > n || n < 0) return 0; double r = 1; for (int i = 1; i <= k; i++) r = r * (n - k + i) / i; return r; }
static long long Tformula(int n, int k) {
    int m = n - 4; long long s = 0;
    for (int d = 0; d <= m; d++) { int u = d + k - 3; if (u < 0) continue; s += (long long)(binom(m, u + d) * binom(u + 1, d) + 0.5); }
    return s;
}

static u128 *shelor_level(int n, int dofull, size_t *ncodes) {
    /* parallel over prefixes of length pl */
    int m = n - 4, pl = m < 6 ? m : 6;
    GS *pre = NULL; size_t np = 0, cp = 0;
    {   /* collect prefixes */
        GS stack[32]; int lvl = 0; (void)lvl;
        /* iterative DFS via recursion helper */
        GS seed; gs_seed(&seed);
        /* simple BFS by levels */
        GS *cur = malloc(sizeof(GS)); cur[0] = seed; size_t nc = 1;
        for (int l = 0; l < pl; l++) {
            GS *nx = malloc(sizeof(GS) * nc * 3); size_t nn = 0;
            for (size_t i = 0; i < nc; i++) for (int j = 0; j < 5; j++) { GS h; if (gs_apply(&cur[i], "abcde"[j], &h)) nx[nn++] = h; }
            free(cur); cur = nx; nc = nn;
        }
        pre = cur; np = nc; cp = nc; (void)cp; (void)stack;
    }
    int nth = omp_get_max_threads();
    GAcc *acc = calloc(nth, sizeof(GAcc));
#pragma omp parallel
    {
        int t = omp_get_thread_num();
        msk *d1 = dofull ? malloc(sizeof(msk) << n) : NULL, *d2 = dofull ? malloc(sizeof(msk) << n) : NULL;
#pragma omp for schedule(dynamic, 1)
        for (long i = 0; i < (long)np; i++) gs_visit(&pre[i], n, dofull, &acc[t], d1, d2);
        free(d1); free(d2);
    }
    GAcc T; memset(&T, 0, sizeof T);
    size_t tot = 0; for (int t = 0; t < nth; t++) tot += acc[t].nc;
    u128 *codes = dofull ? malloc(sizeof(u128) * (tot + 1)) : NULL; size_t p = 0;
    for (int t = 0; t < nth; t++) {
        T.words += acc[t].words; T.hc_bad += acc[t].hc_bad; T.cut_formula_bad += acc[t].cut_formula_bad;
        T.noin1 += acc[t].noin1; T.noout1 += acc[t].noout1; T.inserted_not_in1 += acc[t].inserted_not_in1;
        for (int k = 0; k <= MAXN; k++) T.hist[k] += acc[t].hist[k];
        if (dofull) { memcpy(codes + p, acc[t].codes, sizeof(u128) * acc[t].nc); p += acc[t].nc; }
        free(acc[t].codes);
    }
    size_t distinct = dofull ? uniq128(codes, tot) : 0;
    printf("SHELOR n=%d: words %lld (F(2n-6) = %lld)", n, T.words, fibo(2 * n - 6));
    if (dofull) printf(", hc != 1: %lld, distinct classes %zu", T.hc_bad, distinct);
    printf(", no in-degree-1 vertex: %lld, no out-degree-1 vertex: %lld, inserted vertex not of in-degree 1: %lld, cut-arc count != 3 + #U - #D: %lld\n",
           T.noin1, T.noout1, T.inserted_not_in1, T.cut_formula_bad);
    printf("   T(n,k) direct:");
    int fbad = 0;
    for (int k = 0; k <= n; k++) if (T.hist[k]) printf(" %d:%lld", k, T.hist[k]);
    for (int k = 0; k <= n; k++) if (T.hist[k] != Tformula(n, k)) fbad++;
    printf("   | closed form sum_d C(m,k-3+2d) C(k-2+d,d) mismatches: %d\n", fbad);
    printf("   checks: T(n,2) = 2^(n-5)? %s  T(n,3) = n 2^(n-7)? %s  T(n,4) = (n^2+n-8) 2^(n-10)? %s  top 1, n-4, C(n-3,2)? %s\n",
           (n >= 5 && T.hist[2] == (1LL << (n - 5))) ? "yes" : "no",
           (n >= 7 && T.hist[3] * 128 == (long long)n << n) ? "yes" : (n >= 6 && T.hist[3] * 128 == ((long long)n << n) ? "yes" : "no"),
           (n >= 7 && T.hist[4] * 1024 == ((long long)(n * n + n - 8)) << n) ? "yes" : "no",
           (n >= 5 && T.hist[n - 1] == 1 && T.hist[n - 2] == n - 4 && T.hist[n - 3] == (long long)binom(n - 3, 2)) ? "yes" : "no");
    fflush(stdout);
    free(pre);
    *ncodes = distinct;
    return codes;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage\n"); return 1; }
    printf("aud_uh threads=%d\n", omp_get_max_threads());
    if (!strcmp(argv[1], "all")) {
        int NB = atoi(argv[2]), NE = atoi(argv[3]), NS = atoi(argv[4]), NC = atoi(argv[5]);
        u128 *bset[MAXN + 1] = {0}; size_t bn[MAXN + 1] = {0};
        u128 *eset[MAXN + 1] = {0}; size_t en[MAXN + 1] = {0};
        for (int n = 3; n <= NB; n++) { double t0 = omp_get_wtime(); bset[n] = brute(n, &bn[n]); printf("   [%.1f s]\n", omp_get_wtime() - t0); fflush(stdout); }
        /* EXT from C3 */
        Level L; L.codes = malloc(sizeof(u128)); { msk o[3] = {2, 4, 1}; int c[3] = {0, 1, 2}; L.codes[0] = rotcode(3, o, c); } L.n = 1;
        eset[3] = malloc(sizeof(u128)); eset[3][0] = L.codes[0]; en[3] = 1;
        for (int n = 4; n <= NE; n++) {
            double t0 = omp_get_wtime();
            long long e1, e2, ni, cand;
            Level N = ext_level(n, &L, &e1, &e2, &ni, &cand);
            printf("EXT n=%d: out-sets tried %lld, hc=1 extensions with |O| or |I| = 1: %lld, with |O|,|I| >= 2: %lld, hc=1 with O not a cyclic interval: %lld, classes %zu (Garey %lld) [%.1f s]\n",
                   n, cand, e1, e2, ni, N.n, fibo(2 * n - 6), omp_get_wtime() - t0);
            fflush(stdout);
            eset[n] = N.codes; en[n] = N.n;
            L = N;
        }
        for (int n = 4; n <= NS || n <= NC; n++) {
            if (n < 4) continue;
            double t0 = omp_get_wtime();
            size_t nc;
            u128 *s = shelor_level(n, n <= NS, &nc);
            if (n <= NS) {
                if (n <= NE && en[n]) printf("   SHELOR set == EXT set at n=%d: %s\n", n, (nc == en[n] && !memcmp(s, eset[n], sizeof(u128) * nc)) ? "YES" : "NO");
                if (n <= NB && bn[n]) printf("   SHELOR set == BRUTE set at n=%d: %s\n", n, (nc == bn[n] && !memcmp(s, bset[n], sizeof(u128) * nc)) ? "YES" : "NO");
            }
            printf("   [%.1f s]\n", omp_get_wtime() - t0); fflush(stdout);
            free(s);
            if (n > NS && n >= NC) break;
        }
        for (int n = 3; n <= NB && n <= NE; n++) if (bn[n]) printf("BRUTE set == EXT set at n=%d: %s\n", n, (bn[n] == en[n] && !memcmp(bset[n], eset[n], sizeof(u128) * bn[n])) ? "YES" : "NO");
    }
    return 0;
}
