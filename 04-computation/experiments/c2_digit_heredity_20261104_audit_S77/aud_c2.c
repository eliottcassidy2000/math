/*
 * aud_c2.c -- independent audit program for THM-4873 / HYP-9466 (audit of opus-2026-11-04-S77).
 * Written from scratch for the audit; shares no code with c2_digit_heredity_20261104.c,
 * anymark_c2_digits_20261104.c or ternary_parent_digit_families_20261104.c.
 *
 * Conventions: a tournament on {0..n-1} is stored as out-neighbour bitmasks out[v].
 *   h(T) = # Hamiltonian paths (h(empty) = 1); s_v, e_v = # Hamiltonian paths starting / ending at v;
 *   B_v(T) = sum over W containing v of h(T[V \ W]) * s_v(T[W]);  digit sigma_v = B_v - h.
 *   C2(B, a): x = n with out(x) = {a}; y = n+1 with out(y) = {x}; every old vertex beats y; every old v != a beats x.
 *
 * Methods (deliberately different from the audited programs):
 *   - classes: canonical code = lexicographic minimum of the upper-triangle code over all relabelings that list the
 *     vertices by nondecreasing score (exhaustive inside equal-score cells, prefix pruning only).  Exact, because an
 *     isomorphism preserves scores, so isomorphic tournaments have the same set of score-sorted codes.
 *   - path counts: ONE "start" DP r[W][v] = # Hamiltonian paths of T[W] starting at v; it gives s_v(T[W]) and h(T[W])
 *     for every W at once, and B_v = sum_W h(V \ W) r[W][v].  End counts by a separate forward DP.
 *
 * Modes:
 *   classes NMAX                      class counts (expect 1 1 2 4 12 56 456 6880 191536)
 *   heredity NMIN NMAX                C2 child / grandchild digit census over all classes and all marks
 *   endcount NLAB NCLS                end-count lemma: all labelled tournaments of order <= NLAB, all classes <= NCLS
 *   sweep NMAX LIM                    any-mark one-step sweep h + 2 B_a and the any-mark ready spectra R_1, R_2
 *   triples N                         rows "h e B strong mult" for every (class, vertex) of order N (klein's format)
 *   family TMAX BOUND                 the S77 explicit digit-1 family (DP check of sigma_a = 1 and of H)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef long long ll;
typedef unsigned long long u64;
#define MAXN 14

/* ------------------------------------------------------------------ canonical form ---------------------------- */
typedef struct {
    int n, tot; const uint32_t *out; int req[MAXN]; int perm[MAXN]; int used[MAXN]; int sc[MAXN];
    u64 best; int have;
} Canon;

static void canon_rec(Canon *C, int q, u64 code, int nbits) {
    if (q == C->n) { if (!C->have || code < C->best) { C->best = code; C->have = 1; } return; }
    for (int v = 0; v < C->n; v++) {
        if (C->used[v] || C->sc[v] != C->req[q]) continue;
        u64 c = code;
        for (int i = 0; i < q; i++) c = (c << 1) | ((C->out[C->perm[i]] >> v) & 1u);
        int nb = nbits + q;
        if (C->have && nb > 0) { u64 bp = C->best >> (C->tot - nb); if (c > bp) continue; }
        C->used[v] = 1; C->perm[q] = v;
        canon_rec(C, q + 1, c, nb);
        C->used[v] = 0;
    }
}

static u64 canon(int n, const uint32_t *out) {
    Canon C; memset(&C, 0, sizeof C);
    C.n = n; C.tot = n * (n - 1) / 2; C.out = out;
    for (int v = 0; v < n; v++) C.sc[v] = __builtin_popcount(out[v]);
    for (int v = 0; v < n; v++) C.req[v] = C.sc[v];
    for (int i = 1; i < n; i++) { int x = C.req[i], j = i - 1; while (j >= 0 && C.req[j] > x) { C.req[j + 1] = C.req[j]; j--; } C.req[j + 1] = x; }
    canon_rec(&C, 0, 0, 0);
    return C.best;
}

static void decode(int n, u64 code, uint32_t *out) {
    int b = n * (n - 1) / 2 - 1;
    for (int i = 0; i < n; i++) out[i] = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) {
        int bit = (int)((code >> b) & 1); b--;
        if (bit) out[i] |= 1u << j; else out[j] |= 1u << i;
    }
}

/* simple open-addressing set of u64 codes */
typedef struct { u64 *k; char *u; size_t cap, cnt; } Set;
static void set_init(Set *s, size_t cap) { s->cap = cap; s->cnt = 0; s->k = malloc(cap * sizeof(u64)); s->u = calloc(cap, 1); }
static int set_add(Set *s, u64 key) {
    size_t i = (size_t)((key * 0x9E3779B97F4A7C15ULL) >> 17) % s->cap;
    while (s->u[i]) { if (s->k[i] == key) return 0; if (++i == s->cap) i = 0; }
    s->u[i] = 1; s->k[i] = key; s->cnt++; return 1;
}

static u64 *CL[MAXN]; static size_t NCL[MAXN];
static int cmp_u64(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }

static void build_classes(int nmax) {
    { uint32_t one[1] = {0}; CL[1] = malloc(sizeof(u64)); CL[1][0] = canon(1, one); NCL[1] = 1; }
    for (int n = 2; n <= nmax; n++) {
        size_t cap = NCL[n - 1] * ((size_t)1 << (n - 1)) * 2 + 1024;
        Set s; set_init(&s, cap);
        for (size_t i = 0; i < NCL[n - 1]; i++) {
            uint32_t o[MAXN]; decode(n - 1, CL[n - 1][i], o);
            for (uint32_t O = 0; O < (1u << (n - 1)); O++) {
                uint32_t t[MAXN];
                for (int v = 0; v < n - 1; v++) t[v] = o[v] | (((O >> v) & 1u) ? 0u : (1u << (n - 1)));
                t[n - 1] = O;
                set_add(&s, canon(n, t));
            }
        }
        CL[n] = malloc(s.cnt * sizeof(u64)); NCL[n] = 0;
        for (size_t i = 0; i < s.cap; i++) if (s.u[i]) CL[n][NCL[n]++] = s.k[i];
        qsort(CL[n], NCL[n], sizeof(u64), cmp_u64);
        free(s.k); free(s.u);
    }
    /* sanity: canon(decode(code)) == code */
    for (int n = 1; n <= nmax; n++) for (size_t i = 0; i < NCL[n]; i++) {
        uint32_t o[MAXN]; decode(n, CL[n][i], o);
        if (canon(n, o) != CL[n][i]) { fprintf(stderr, "canon round-trip failure n=%d\n", n); exit(1); }
    }
}

/* ------------------------------------------------------------------ path DPs ----------------------------------- */
/* r[W*n+v] = # Hamiltonian paths of T[W] starting at v (v in W); hW[W] = h(T[W]), hW[0] = 1 */
static void start_dp(int n, const uint32_t *out, ll *r, ll *hW) {
    int full = (1 << n) - 1;
    hW[0] = 1;
    for (int W = 1; W <= full; W++) {
        ll tot = 0;
        for (uint32_t m = (uint32_t)W; m; m &= m - 1) {
            int v = __builtin_ctz(m); int Wv = W & ~(1 << v); ll c;
            if (!Wv) c = 1;
            else { c = 0; for (uint32_t o = out[v] & (uint32_t)Wv; o; o &= o - 1) c += r[(size_t)Wv * n + __builtin_ctz(o)]; }
            r[(size_t)W * n + v] = c; tot += c;
        }
        hW[W] = tot;
    }
}
/* B[v] and sigma[v] for every v; returns H; also starts s[v] if s != NULL */
static ll all_digits(int n, const uint32_t *out, ll *r, ll *hW, ll *Bv, ll *sig, ll *s) {
    start_dp(n, out, r, hW);
    int full = (1 << n) - 1;
    for (int v = 0; v < n; v++) Bv[v] = 0;
    for (int W = 1; W <= full; W++) {
        ll hc = hW[full ^ W];
        for (uint32_t m = (uint32_t)W; m; m &= m - 1) { int v = __builtin_ctz(m); Bv[v] += hc * r[(size_t)W * n + v]; }
    }
    ll H = hW[full];
    for (int v = 0; v < n; v++) sig[v] = Bv[v] - H;
    if (s) for (int v = 0; v < n; v++) s[v] = r[(size_t)full * n + v];
    return H;
}
/* end counts by forward DP f[W*n+v] = # Hamiltonian paths of T[W] ending at v */
static ll end_counts(int n, const uint32_t *out, ll *f, ll *e) {
    int full = (1 << n) - 1; uint32_t in[MAXN];
    for (int v = 0; v < n; v++) in[v] = (~out[v]) & (uint32_t)full & ~(1u << v);
    ll H = 0;
    for (int W = 1; W <= full; W++)
        for (uint32_t m = (uint32_t)W; m; m &= m - 1) {
            int v = __builtin_ctz(m); int Wv = W & ~(1 << v); ll c;
            if (!Wv) c = 1;
            else { c = 0; for (uint32_t o = in[v] & (uint32_t)Wv; o; o &= o - 1) c += f[(size_t)Wv * n + __builtin_ctz(o)]; }
            f[(size_t)W * n + v] = c;
        }
    for (int v = 0; v < n; v++) { e[v] = f[(size_t)full * n + v]; H += e[v]; }
    return H;
}

/* reachability and strong components */
static void reach_sets(int n, const uint32_t *out, uint32_t *R) {
    for (int v = 0; v < n; v++) {
        uint32_t seen = 1u << v, fr = seen;
        while (fr) { uint32_t nf = 0; for (uint32_t m = fr; m; m &= m - 1) nf |= out[__builtin_ctz(m)]; nf &= ~seen; seen |= nf; fr = nf; }
        R[v] = seen;
    }
}
static int is_strong(int n, const uint32_t *out) {
    uint32_t R[MAXN]; reach_sets(n, out, R); uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    for (int v = 0; v < n; v++) if (R[v] != full) return 0;
    return 1;
}
/* component of v = {u : u in R[v] and v in R[u]}; returns mask of the last component and the number of components */
static uint32_t last_component(int n, const uint32_t *out, int *ncomp, int *all_other_singletons) {
    uint32_t R[MAXN]; reach_sets(n, out, R); uint32_t comp[MAXN]; uint32_t last = 0; int nc = 0, okk = 1; uint32_t seen = 0;
    for (int v = 0; v < n; v++) { comp[v] = 0; for (int u = 0; u < n; u++) if (((R[v] >> u) & 1) && ((R[u] >> v) & 1)) comp[v] |= 1u << u; }
    for (int v = 0; v < n; v++) { if (R[v] == comp[v]) last = comp[v]; if (!((seen >> v) & 1)) { nc++; seen |= comp[v]; } }
    for (int v = 0; v < n; v++) if (!((last >> v) & 1) && __builtin_popcount(comp[v]) != 1) okk = 0;
    *ncomp = nc; *all_other_singletons = okk; return last;
}

static void c2child(int n, const uint32_t *out, int a, uint32_t *t) {
    int x = n, y = n + 1;
    for (int v = 0; v < n; v++) { t[v] = out[v] | (1u << y); if (v != a) t[v] |= 1u << x; }
    t[x] = 1u << a; t[y] = 1u << x;
}

static void scoreseq(int n, const uint32_t *out, char *buf) {
    int s[MAXN]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(out[v]);
    for (int i = 1; i < n; i++) { int x = s[i], j = i - 1; while (j >= 0 && s[j] > x) { s[j + 1] = s[j]; j--; } s[j + 1] = x; }
    int p = 0; buf[p++] = '(';
    for (int i = 0; i < n; i++) p += sprintf(buf + p, "%s%d", i ? "," : "", s[i]);
    buf[p++] = ')'; buf[p] = 0;
}

/* ------------------------------------------------------------------ heredity ----------------------------------- */
static void mode_heredity(int nmin, int nmax) {
    build_classes(nmax);
    int NB = nmax + 4;
    ll *r = malloc(sizeof(ll) * ((size_t)1 << NB) * NB), *hW = malloc(sizeof(ll) << NB), *f = malloc(sizeof(ll) * ((size_t)1 << NB) * NB);
    long long tot_pairs = 0, tot_g = 0;
    for (int n = nmin; n <= nmax; n++) {
        long pairs = 0, lawfail1 = 0, lawfail2 = 0, natx = 0, naty = 0, natx2 = 0, naty2 = 0;
        long ch1 = 0, ch2 = 0, ch_old12 = 0, ready_viol = 0, ready_cases = 0, below_h = 0;
        long g = 0, g1 = 0, g2 = 0, g12 = 0, hist[5] = {0}, g2_oldmark = 0, g_old12_marks = 0, t1strong = 0;
        long dmin_viol = 0;
        printf("== |B| = %d: %zu classes\n", n, NCL[n]);
        for (size_t ci = 0; ci < NCL[n]; ci++) {
            uint32_t B[MAXN]; decode(n, CL[n][ci], B);
            ll BB[MAXN], sB[MAXN], eB[MAXN];
            ll h = all_digits(n, B, r, hW, BB, sB, NULL);
            ll h2 = end_counts(n, B, f, eB);
            if (h2 != h) { printf("H MISMATCH\n"); exit(1); }
            int strongB = is_strong(n, B);
            int transB = (h == 1);
            int sink = 0; for (int v = 0; v < n; v++) if (B[v] == 0) sink = 1;
            int ttc3 = (h == 3 && !sink);
            char sq[128]; scoreseq(n, B, sq);
            for (int a = 0; a < n; a++) {
                uint32_t T1[MAXN]; c2child(n, B, a, T1);
                int m = n + 2; ll B1[MAXN], s1[MAXN], e1[MAXN];
                ll H1 = all_digits(m, T1, r, hW, B1, s1, NULL);
                end_counts(m, T1, f, e1);
                pairs++;
                if (H1 != h + 2 * BB[a]) lawfail1++;
                if (s1[n] != 0) natx++;
                if (s1[n + 1] != h - eB[a]) naty++;
                /* out-degree bound sigma_v >= d(d-1) (dual of READY (7)) on the child */
                for (int u = 0; u < m; u++) { ll d = __builtin_popcount(T1[u]); if (s1[u] < d * (d - 1)) dmin_viol++; }
                int s1strong = is_strong(m, T1); t1strong += s1strong;
                int c1 = 0, c2 = 0; ll mo = -1;
                for (int u = 0; u < m; u++) { if (s1[u] == 1) c1 = 1; if (s1[u] == 2) c2 = 1; }
                for (int u = 0; u < n; u++) { if (s1[u] == 1 || s1[u] == 2) ch_old12++; if (mo < 0 || s1[u] < mo) mo = s1[u]; }
                if (n >= 1 && mo < h) below_h++;
                ch1 += c1; ch2 += c2;
                if (strongB && n >= 4) { ready_cases++; if (c1 || c2) ready_viol++; }
                if (c1 || c2) {
                    printf("  child: B#%zu score%s h=%lld %s%s%s a=%d (outdeg_B(a)=%d, e_a=%lld) H1=%lld small marks:", ci, sq, h,
                           transB ? "transitive" : "", ttc3 ? "TT=>C3" : "", strongB ? " strong" : "", a, __builtin_popcount(B[a]), eB[a], H1);
                    for (int u = 0; u < m; u++) if (s1[u] == 1 || s1[u] == 2)
                        printf(" %s(digit %lld)", u == n ? "x" : u == n + 1 ? "y" : u == a ? "a" : "old", s1[u]);
                    printf("\n");
                }
                /* grandchildren */
                for (int u = 0; u < m; u++) {
                    uint32_t T2[MAXN]; c2child(m, T1, u, T2);
                    int k = m + 2; ll B2[MAXN], s2[MAXN];
                    ll H2 = all_digits(k, T2, r, hW, B2, s2, NULL);
                    g++;
                    if (H2 != H1 + 2 * B1[u]) lawfail2++;
                    if (s2[m] != 0) natx2++;
                    if (s2[m + 1] != H1 - e1[u]) naty2++;
                    int d1 = 0, d2 = 0; ll sm = -1; int oldsmall = 0;
                    for (int w = 0; w < k; w++) {
                        if (s2[w] == 1) d1 = 1;
                        if (s2[w] == 2) d2 = 1;
                        if (s2[w] > 0 && (sm < 0 || s2[w] < sm)) sm = s2[w];
                        if (w < m && (s2[w] == 1 || s2[w] == 2)) { oldsmall = 1; g_old12_marks++; }
                    }
                    g1 += d1; g2 += d2; g12 += (d1 && d2); g2_oldmark += oldsmall;
                    int bucket = sm < 0 ? 0 : sm <= 2 ? 1 : sm <= 5 ? 2 : sm <= 20 ? 3 : 4; hist[bucket]++;
                    if (d1 || d2) {
                        printf("  grandchild: B#%zu score%s a=%d u=%s H1=%lld H2=%lld small marks:", ci, sq, a,
                               u == n ? "x1" : u == n + 1 ? "y1" : u == a ? "a" : "old", H1, H2);
                        for (int w = 0; w < k; w++) if (s2[w] == 1 || s2[w] == 2)
                            printf(" %s(digit %lld)", w == m ? "x2" : w == m + 1 ? "y2" : w == u ? "u[old mark of T1]" : w == a ? "a" : w == n ? "x1" : w == n + 1 ? "y1" : "old", s2[w]);
                        printf("\n");
                    }
                }
            }
        }
        tot_pairs += pairs; tot_g += g;
        printf("|B|=%d: pairs %ld | C2 law failures child %ld grandchild %ld | native sigma_x!=0: %ld, sigma_y!=h-e_a: %ld (grandchild level: %ld, %ld)\n",
               n, pairs, lawfail1, lawfail2, natx, naty, natx2, naty2);
        printf("        children with a digit-1 mark %ld, digit-2 mark %ld; old marks with digit 1 or 2: %ld; some old mark below h: %ld; strong children: %ld\n",
               ch1, ch2, ch_old12, below_h, t1strong);
        printf("        READY s4 dual check (strong B, |B|>=4): %ld cases, %ld with a digit-1/2 mark in the child; sigma_v < d(d-1) violations in children: %ld\n",
               ready_cases, ready_viol, dmin_viol);
        printf("        grandchildren %ld | with digit-1 mark %ld, digit-2 mark %ld, both %ld | grandchildren with an OLD mark of digit 1/2: %ld (old marks: %ld)\n",
               g, g1, g2, g12, g2_oldmark, g_old12_marks);
        printf("        smallest positive digit: none %ld, <=2 %ld, 3-5 %ld, 6-20 %ld, >20 %ld\n", hist[0], hist[1], hist[2], hist[3], hist[4]);
        fflush(stdout);
    }
    printf("TOTAL over |B| = %d..%d: pairs %lld, grandchildren %lld\n", nmin, nmax, tot_pairs, tot_g);
}

/* ------------------------------------------------------------------ end-count lemma ---------------------------- */
static void endcount_check(int n, const uint32_t *out, ll *r, ll *hW, ll *f, long *st) {
    ll e[MAXN], s[MAXN], Bv[MAXN], sig[MAXN];
    ll H = end_counts(n, out, f, e);
    ll H2 = all_digits(n, out, r, hW, Bv, sig, s);
    if (H != H2) { st[0]++; return; }
    st[1]++;
    /* components of T and of its converse */
    int nc, oth; uint32_t last = last_component(n, out, &nc, &oth);
    uint32_t conv[MAXN]; uint32_t full = (1u << n) - 1;
    for (int v = 0; v < n; v++) conv[v] = (~out[v]) & full & ~(1u << v);
    int ncc, othc; uint32_t first = last_component(n, conv, &ncc, &othc); /* last component of the converse = first of T */
    int ttc3 = (__builtin_popcount(last) == 3 && oth);   /* TT_j => C3 */
    int c3tt = (__builtin_popcount(first) == 3 && othc); /* C3 => TT_k */
    for (int v = 0; v < n; v++) {
        if (H >= 3 && e[v] == H - 1) st[2]++;
        if (H >= 3 && s[v] == H - 1) st[3]++;
        int lhs = (e[v] == H - 2), rhs = ttc3 && ((last >> v) & 1);
        if (lhs != rhs) st[4]++;
        if (lhs) st[5]++;
        int lhs2 = (s[v] == H - 2), rhs2 = c3tt && ((first >> v) & 1);
        if (lhs2 != rhs2) st[6]++;
        if (H == 1 && e[v] > 1) st[7]++;
        ll d = __builtin_popcount(out[v]);
        if (sig[v] < d * (d - 1)) st[8]++;
        if (sig[v] == 1 && d != 1) st[9]++;
        if (sig[v] == 2 && d > 2) st[10]++;
        if (sig[v] < 0) st[11]++;
        if (H >= 3) { ll gap = H - e[v]; if (gap >= 3 && gap <= 6) st[12 + (int)(gap - 3)]++; }
    }
}

static void mode_endcount(int nlab, int ncls) {
    int NB = (nlab > ncls ? nlab : ncls);
    ll *r = malloc(sizeof(ll) * ((size_t)1 << NB) * NB), *hW = malloc(sizeof(ll) << NB), *f = malloc(sizeof(ll) * ((size_t)1 << NB) * NB);
    for (int n = 1; n <= nlab; n++) {
        long st[20] = {0};
        int tot = n * (n - 1) / 2;
        for (u64 code = 0; code < (1ULL << tot); code++) { uint32_t o[MAXN]; decode(n, code, o); endcount_check(n, o, r, hW, f, st); }
        printf("labelled order %d: %ld tournaments (H mismatch %ld) | e_v=H-1 with H>=3: %ld | s_v=H-1 with H>=3: %ld | e_v=H-2 <-> (TT=>C3, v in C3) disagreements: %ld (e_v=H-2 occurrences %ld) | start version disagreements: %ld | H=1 with e_v>1: %ld | sigma_v<d(d-1): %ld, sigma_v=1 with outdeg!=1: %ld, sigma_v=2 with outdeg>2: %ld, sigma<0: %ld | e_v=H-3..H-6 occurrences: %ld %ld %ld %ld\n",
               n, st[1], st[0], st[2], st[3], st[4], st[5], st[6], st[7], st[8], st[9], st[10], st[11], st[12], st[13], st[14], st[15]);
        fflush(stdout);
    }
    if (ncls > 0) {
        build_classes(ncls);
        for (int n = 1; n <= ncls; n++) {
            long st[20] = {0};
            for (size_t i = 0; i < NCL[n]; i++) { uint32_t o[MAXN]; decode(n, CL[n][i], o); endcount_check(n, o, r, hW, f, st); }
            printf("classes order %d: %ld classes (H mismatch %ld) | e_v=H-1 with H>=3: %ld | s_v=H-1: %ld | e_v=H-2 disagreements: %ld (occurrences %ld) | start disagreements: %ld | H=1, e_v>1: %ld | sigma_v<d(d-1): %ld, sigma_v=1 & outdeg!=1: %ld, sigma_v=2 & outdeg>2: %ld, sigma<0: %ld | e_v=H-3..H-6: %ld %ld %ld %ld\n",
                   n, st[1], st[0], st[2], st[3], st[4], st[5], st[6], st[7], st[8], st[9], st[10], st[11], st[12], st[13], st[14], st[15]);
            fflush(stdout);
        }
    }
}

/* ------------------------------------------------------------------ sweep -------------------------------------- */
static void mode_sweep(int nmax, ll lim) {
    build_classes(nmax);
    ll VL = 3 * lim + 10;
    char *V = calloc(VL + 1, 1), *Hs = calloc(lim + 1, 1), *R0 = calloc(lim + 1, 1), *R1 = calloc(lim + 1, 1), *R2 = calloc(lim + 1, 1);
    char *R1s = calloc(lim + 1, 1), *R2s = calloc(lim + 1, 1);  /* restricted to strong B with a sole exit a */
    ll *r = malloc(sizeof(ll) * ((size_t)1 << nmax) * nmax), *hW = malloc(sizeof(ll) << nmax);
    size_t ncls = 0;
    for (int n = 1; n <= nmax; n++) {
        ncls += NCL[n];
        for (size_t i = 0; i < NCL[n]; i++) {
            uint32_t o[MAXN]; decode(n, CL[n][i], o);
            ll Bv[MAXN], sig[MAXN]; ll h = all_digits(n, o, r, hW, Bv, sig, NULL);
            int st = is_strong(n, o);
            if (h <= lim) Hs[h] = 1;
            for (int a = 0; a < n; a++) {
                ll v = h + 2 * Bv[a]; if (v <= VL) V[v] = 1;
                if (h <= lim) {
                    if (sig[a] == 0) R0[h] = 1;
                    if (sig[a] == 1) R1[h] = 1;
                    if (sig[a] == 2) R2[h] = 1;
                    if (st && __builtin_popcount(o[a]) == 1) { if (sig[a] == 1) R1s[h] = 1; if (sig[a] == 2) R2s[h] = 1; }
                }
            }
        }
    }
    printf("sweep over %zu classes of orders 1..%d (every mark)\n", ncls, nmax);
    printf("odd v <= %lld NOT of the form h + 2 B_a:", lim);
    long miss = 0; for (ll v = 1; v <= lim; v += 2) if (!V[v]) { miss++; printf(" %lld", v); }
    printf("  [%ld]\n", miss);
    ll first_gap = -1; for (ll v = 69; v <= VL; v += 2) if (!V[v]) { first_gap = v; break; }
    printf("every odd v in [69, %lld] attained (first gap above 69: %lld)\n", first_gap - 2, first_gap);
    long nh = 0, n1 = 0, n2 = 0, n12 = 0, n1s = 0, n2s = 0; ll hmax = 0;
    for (ll v = 1; v <= lim; v += 2) if (Hs[v]) { nh++; hmax = v; if (R1[v]) n1++; if (R2[v]) n2++; if (R1[v] && R2[v]) n12++; if (R1s[v]) n1s++; if (R2s[v]) n2s++; }
    printf("realised odd h <= %lld: %ld (largest %lld); with a digit-1 mark: %ld; with a digit-2 mark: %ld; both: %ld\n", lim, nh, hmax, n1, n2, n12);
    printf("restricted to strong B with a sole-exit mark: digit 1: %ld values, digit 2: %ld values\n", n1s, n2s);
    printf("h with a digit-1 mark:"); for (ll v = 1; v <= lim; v += 2) if (R1[v]) printf(" %lld", v); printf("\n");
    printf("h with a digit-1 mark but not via a strong sole exit:"); for (ll v = 1; v <= lim; v += 2) if (R1[v] && !R1s[v]) printf(" %lld", v); printf("\n");
    printf("h with a digit-2 mark but not via a strong sole exit:"); for (ll v = 1; v <= lim; v += 2) if (R2[v] && !R2s[v]) printf(" %lld", v); printf("\n");
    /* Ternary parents.  N odd, sigma = 2N mod 3, h = (N - 2 sigma)/3.  A digit-1/2 mark forces T = TT => S with S strong
     * and H(T) = H(S) (or T transitive), and Busch's bound H(S) >= 5^((|S|-1)/3) puts every strong S with H(S) <= 42 at
     * order <= 8, so for h <= 42 (N <= 130) the classes of order <= 8 decide whether h has a mark of digit sigma. */
    if (nmax >= 8) {
        printf("odd N <= 130 with NO ternary parent (unique odd predecessor h lacks the digit in every tournament):");
        for (ll N = 3; N <= 130; N += 2) {
            int sg = (int)((2 * N) % 3); ll h = (N - 2 * sg) / 3;
            int ok = (h >= 1 && h <= lim && Hs[h] && (sg == 0 || (sg == 1 && R1[h]) || (sg == 2 && R2[h])));
            if (!ok) printf(" %lld[h=%lld,digit %d%s]", N, h, sg, V[N] ? ",attained by h+2B_a" : "");
        }
        printf("\n");
    }
}

/* ------------------------------------------------------------------ triples ------------------------------------ */
typedef struct { ll h, e, B; int st; } Row;
static int cmp_row(const void *a, const void *b) {
    const Row *x = a, *y = b;
    if (x->h != y->h) return x->h < y->h ? -1 : 1;
    if (x->e != y->e) return x->e < y->e ? -1 : 1;
    if (x->B != y->B) return x->B < y->B ? -1 : 1;
    return x->st - y->st;
}
static void mode_triples(int N) {
    build_classes(N);
    ll *r = malloc(sizeof(ll) * ((size_t)1 << N) * N), *hW = malloc(sizeof(ll) << N), *f = malloc(sizeof(ll) * ((size_t)1 << N) * N);
    Row *rows = malloc(sizeof(Row) * NCL[N] * N); size_t nr = 0;
    for (size_t i = 0; i < NCL[N]; i++) {
        uint32_t o[MAXN]; decode(N, CL[N][i], o);
        ll Bv[MAXN], sig[MAXN], e[MAXN]; ll h = all_digits(N, o, r, hW, Bv, sig, NULL); end_counts(N, o, f, e);
        int st = is_strong(N, o);
        for (int a = 0; a < N; a++) { rows[nr].h = h; rows[nr].e = e[a]; rows[nr].B = Bv[a]; rows[nr].st = st; nr++; }
    }
    qsort(rows, nr, sizeof(Row), cmp_row);
    for (size_t i = 0; i < nr;) {
        size_t j = i; while (j < nr && cmp_row(&rows[i], &rows[j]) == 0) j++;
        printf("%lld %lld %lld %d %zu\n", rows[i].h, rows[i].e, rows[i].B, rows[i].st, j - i);
        i = j;
    }
}

/* ------------------------------------------------------------------ family (DP side) --------------------------- */
/* C(b): r=0, w=1, z_i = 1+i; r->w; w->z_i; z_i->z_j (i<j); r->z_i iff b_i=1.  B = C^op + a (a = t+2, a -> r only). */
static void mode_family(int tmax, ll bound) {
    int NB = tmax + 3;
    ll *r = malloc(sizeof(ll) * ((size_t)1 << NB) * NB), *hW = malloc(sizeof(ll) << NB);
    for (int t = 1; t <= tmax; t++) {
        int n = t + 2; long bad = 0;
        for (uint32_t b = 0; b < (1u << t); b++) {
            uint32_t C[MAXN] = {0};
            C[0] |= 1u << 1;
            for (int i = 1; i <= t; i++) {
                C[1] |= 1u << (1 + i);
                for (int j = i + 1; j <= t; j++) C[1 + i] |= 1u << (1 + j);
                if ((b >> (i - 1)) & 1u) C[0] |= 1u << (1 + i); else C[1 + i] |= 1u << 0;
            }
            uint32_t Bt[MAXN] = {0};
            for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (u != v && ((C[v] >> u) & 1u)) Bt[u] |= 1u << v;  /* converse */
            Bt[n] = 1u << 0; for (int v = 1; v < n; v++) Bt[v] |= 1u << n;
            ll Bv[MAXN], sig[MAXN]; ll H = all_digits(n + 1, Bt, r, hW, Bv, sig, NULL);
            if (sig[n] != 1) bad++;
            printf("F %d %u %lld\n", t, b, H);
        }
        fprintf(stderr, "t=%d: digit-1 failures %ld\n", t, bad);
    }
    (void)bound;
}

/* ------------------------------------------------------------------ proofcheck --------------------------------- */
/* Checks, over all classes B of order 2..NMAX and every mark a, the identities used in the audit's all-order proof:
 *   (R12) sigma_a(C2(B,a)) = 2 B_a(B) - h(B) - h(B - a)            [multiport RESULT.md (12): S_a' = 4S - H(B-a)]
 *   (Z)   z = sink of B, z != a:  H(T1 - z) - e_z(T1) = h(B - z) + B_a(B - z)  and  sigma_z(T1) >= 4
 *   (W)   outdeg_B(a) = 1, a -> w:  sigma_a(T1) = e_w(B - a) + B_w(B - a);  outdeg_B(a) = 0: sigma_a(T1) = 0
 *   (D)   every vertex of T1 of out-degree >= 3 has digit >= 6 (consequence of sigma_v >= d(d-1))
 */
static void induced(int n, const uint32_t *out, uint32_t keep, uint32_t *o2, int *map) {
    int k = 0; int idx[MAXN];
    for (int v = 0; v < n; v++) if ((keep >> v) & 1u) { idx[v] = k; map[k] = v; k++; } else idx[v] = -1;
    for (int v = 0; v < n; v++) if ((keep >> v) & 1u) {
        uint32_t m = 0;
        for (int u = 0; u < n; u++) if (((keep >> u) & 1u) && ((out[v] >> u) & 1u)) m |= 1u << idx[u];
        o2[idx[v]] = m;
    }
}
static void mode_proofcheck(int nmax) {
    build_classes(nmax);
    int NB = nmax + 2;
    ll *r = malloc(sizeof(ll) * ((size_t)1 << NB) * NB), *hW = malloc(sizeof(ll) << NB), *f = malloc(sizeof(ll) * ((size_t)1 << NB) * NB);
    for (int n = 2; n <= nmax; n++) {
        long cases = 0, bad12 = 0, zc = 0, badZ = 0, badZ4 = 0, wc = 0, badW = 0, badW0 = 0, badD = 0, small_other = 0;
        for (size_t ci = 0; ci < NCL[n]; ci++) {
            uint32_t B[MAXN]; decode(n, CL[n][ci], B);
            ll BB[MAXN], sB[MAXN]; ll h = all_digits(n, B, r, hW, BB, sB, NULL);
            int sink = -1; for (int v = 0; v < n; v++) if (B[v] == 0) sink = v;
            for (int a = 0; a < n; a++) {
                uint32_t T1[MAXN]; c2child(n, B, a, T1); int m = n + 2;
                ll B1[MAXN], s1[MAXN], e1[MAXN]; ll H1 = all_digits(m, T1, r, hW, B1, s1, NULL); end_counts(m, T1, f, e1);
                cases++;
                /* (R12) */
                uint32_t Ba[MAXN]; int mp[MAXN]; induced(n, B, ((1u << n) - 1) & ~(1u << a), Ba, mp);
                ll tmpB[MAXN], tmps[MAXN]; ll hBa = (n - 1 >= 1) ? all_digits(n - 1, Ba, r, hW, tmpB, tmps, NULL) : 1;
                if (s1[a] != 2 * BB[a] - h - hBa) bad12++;
                /* (Z) */
                if (sink >= 0 && sink != a) {
                    zc++;
                    uint32_t Bz[MAXN]; int mz[MAXN]; induced(n, B, ((1u << n) - 1) & ~(1u << sink), Bz, mz);
                    int az = -1; for (int k = 0; k < n - 1; k++) if (mz[k] == a) az = k;
                    ll bz[MAXN], sz[MAXN]; ll hz = all_digits(n - 1, Bz, r, hW, bz, sz, NULL);
                    uint32_t T1z[MAXN]; int mz1[MAXN]; induced(m, T1, ((1u << m) - 1) & ~(1u << sink), T1z, mz1);
                    ll b1z[MAXN], s1z[MAXN]; ll H1z = all_digits(m - 1, T1z, r, hW, b1z, s1z, NULL);
                    if (H1z - e1[sink] != hz + bz[az]) badZ++;
                    if (s1[sink] < 4) badZ4++;
                }
                /* (W) */
                int da = __builtin_popcount(B[a]);
                if (da == 0 && s1[a] != 0) badW0++;
                if (da == 1) {
                    wc++;
                    int w = __builtin_ctz(B[a]); int wk = -1; for (int k = 0; k < n - 1; k++) if (mp[k] == w) wk = k;
                    ll bw[MAXN], sw[MAXN], ew[MAXN]; all_digits(n - 1, Ba, r, hW, bw, sw, NULL); end_counts(n - 1, Ba, f, ew);
                    if (s1[a] != ew[wk] + bw[wk]) badW++;
                }
                /* (D) and: small digits only at y or a */
                for (int u = 0; u < m; u++) {
                    int d = __builtin_popcount(T1[u]);
                    if (d >= 3 && s1[u] < 6) badD++;
                    if ((s1[u] == 1 || s1[u] == 2) && u != n + 1 && u != a) small_other++;
                }
                (void)H1;
            }
        }
        printf("|B|=%d: %ld (B,a) | (R12) failures %ld | sink cases %ld: identity failures %ld, sigma_z<4: %ld | outdeg-1 marks %ld: (W) failures %ld; outdeg-0 mark sigma!=0: %ld | outdeg>=3 with digit<6: %ld | small digit at a vertex other than y, a: %ld\n",
               n, cases, bad12, zc, badZ, badZ4, wc, badW, badW0, badD, small_other);
        fflush(stdout);
    }
    /* (L) component reduction used in the audit's M4: for every tournament T (classes of order <= nmax) and mark u:
     *   u outside the last strong component S  =>  sigma_u(T) >= H(T);
     *   u in S  =>  sigma_u(T) = h(T - S) * sigma_u(S). */
    for (int n = 1; n <= nmax; n++) {
        long marks = 0, badOut = 0, badIn = 0;
        for (size_t ci = 0; ci < NCL[n]; ci++) {
            uint32_t T[MAXN]; decode(n, CL[n][ci], T);
            ll Bv[MAXN], sg[MAXN]; ll H = all_digits(n, T, r, hW, Bv, sg, NULL);
            int nc, oth; uint32_t L = last_component(n, T, &nc, &oth);
            uint32_t S[MAXN], P[MAXN]; int ms[MAXN], mpp[MAXN];
            int ns = __builtin_popcount(L), np = n - ns;
            induced(n, T, L, S, ms);
            ll bS[MAXN], sS[MAXN]; all_digits(ns, S, r, hW, bS, sS, NULL);
            ll hP = 1;
            if (np > 0) { induced(n, T, ((1u << n) - 1) & ~L, P, mpp); ll bP[MAXN], sP[MAXN]; hP = all_digits(np, P, r, hW, bP, sP, NULL); }
            for (int u = 0; u < n; u++) {
                marks++;
                if (!((L >> u) & 1u)) { if (sg[u] < H) badOut++; }
                else { int k = -1; for (int j = 0; j < ns; j++) if (ms[j] == u) k = j; if (sg[u] != hP * sS[k]) badIn++; }
            }
        }
        printf("(L) order %d: %ld marks | outside last component with sigma < H: %ld | last-component scaling failures: %ld\n", n, marks, badOut, badIn);
    }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: see header\n"); return 1; }
    if (!strcmp(argv[1], "proofcheck")) { mode_proofcheck(atoi(argv[2])); return 0; }
    const char *m = argv[1];
    if (!strcmp(m, "classes")) {
        int nmax = atoi(argv[2]); build_classes(nmax);
        for (int n = 1; n <= nmax; n++) printf("order %d: %zu classes\n", n, NCL[n]);
    } else if (!strcmp(m, "heredity")) mode_heredity(atoi(argv[2]), atoi(argv[3]));
    else if (!strcmp(m, "endcount")) mode_endcount(atoi(argv[2]), atoi(argv[3]));
    else if (!strcmp(m, "sweep")) mode_sweep(atoi(argv[2]), atoll(argv[3]));
    else if (!strcmp(m, "triples")) mode_triples(atoi(argv[2]));
    else if (!strcmp(m, "family")) mode_family(atoi(argv[2]), atoll(argv[3]));
    else { fprintf(stderr, "unknown mode\n"); return 1; }
    return 0;
}
