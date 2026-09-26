/*
 * procgen_seven3_20260926_markov.c -- exact rho_max of a variable-depth sign rule of q n +- 1 by MARKOV REFINEMENT
 * (lane "seven3", session collatz-procgen-20260922, 2026-09-26).  C port of Markov/howard/certify_upper in
 * procgen_seven3_20260926_lib.py (same method, separate implementation).
 *
 * Input (stdin): "q" then lines "c d s" (odd residue c mod 2^d, d <= 62, sign s = +1/-1): a partition of the odd
 * 2-adic integers into classes with constant sign (the rule).  The even integers form the leaf (0 mod 2).
 * Refinement: every leaf L = (c, d) has the image class T(L) = (T(c) mod 2^(d-1), d-1) (T(c) = c/2 or (q c + s)/2);
 * leaves are split (children inherit the sign) until every image is a union of leaves.  Leaf graph: L -> every leaf
 * inside T(L); weight 1 on odd leaves.  Lemma M: closed walks <-> periodic orbits, same density.
 * Max mean cycle: Howard policy iteration (doubles) gives a cycle of density a/p (exact integers); then the least
 * fixed point psi(L) = max(0, e(L) + max_{L'} psi(L')), e = p - a (odd), -a (even), is computed with exact int64
 * and checked on EVERY edge (Lemma G1 on the leaf graph).  If it does not converge (Howard missed a denser cycle),
 * Dinkelbach steps on the pointer graph take over (a pointer cycle created by strict increases is denser than F).
 * Output: "LEAVES n", "RHO a p" (reduced), "CERT ok maxpsi", "CYCLE c_0 d_0 c_1 d_1 ..." (witness walk of leaves).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef int64_t i64;
typedef uint64_t u64;

static i64 Q;
/* trie */
static int32_t *ch0, *ch1; static int8_t *sg; static int8_t *dep; static u64 *res; static int32_t nn, cap_n;
static int32_t newnode(u64 c, int d, int s) {
    if (nn == cap_n) {
        cap_n = cap_n ? 2 * cap_n : 1 << 16;
        ch0 = realloc(ch0, cap_n * 4); ch1 = realloc(ch1, cap_n * 4); sg = realloc(sg, cap_n); dep = realloc(dep, cap_n); res = realloc(res, cap_n * 8);
        if (!ch0 || !ch1 || !sg || !dep || !res) { fprintf(stderr, "oom\n"); exit(2); }
    }
    ch0[nn] = ch1[nn] = -1; sg[nn] = (int8_t)s; dep[nn] = (int8_t)d; res[nn] = c; return nn++;
}
static inline u64 mask(int d) { return d >= 64 ? ~0ULL : ((1ULL << d) - 1); }
/* worklist of leaves */
static int32_t *wl; static i64 wn, wcap;
static void wpush(int32_t x) { if (wn == wcap) { wcap = wcap ? 2 * wcap : 1 << 16; wl = realloc(wl, wcap * 4); } wl[wn++] = x; }

static void split(int32_t x) {
    u64 c = res[x]; int d = dep[x];
    int s = sg[x];
    int32_t a = newnode(c, d + 1, s), b = newnode(c | (1ULL << d), d + 1, s);
    ch0[x] = a; ch1[x] = b;   /* x now internal; sign kept only for leaves */
    wpush(a); wpush(b);
}
static int32_t descend(u64 c, int d) { /* node on the path to (c, d) that is a leaf or the node of depth d */
    int32_t x = 0;
    while (dep[x] < d && ch0[x] >= 0) x = ((c >> dep[x]) & 1) ? ch1[x] : ch0[x];
    return x;
}
static void ensure_union(u64 c, int d) {
    int32_t x = descend(c, d);
    while (dep[x] < d) { split(x); x = ((c >> dep[x]) & 1) ? ch1[x] : ch0[x]; }
}
static inline u64 image_res(int32_t x) {
    u64 c = res[x]; int d = dep[x];
    u64 t;
    if (!(c & 1)) t = c >> 1;
    else {  /* (q c + s)/2 mod 2^(d-1): compute q*c + s mod 2^d with wraparound arithmetic */
        u64 m = (u64)Q * c + (u64)(i64)sg[x];
        t = (m & mask(d)) >> 1;
    }
    return t & mask(d - 1);
}

/* leaf enumeration */
static int32_t *lid; /* node -> leaf index or -1 */
static int32_t *leaves; static int32_t nl;
static int32_t *eoff; static int32_t *edst; static i64 ne;

static void collect_under(int32_t x, int32_t **buf, i64 *bn, i64 *bc) {
    /* iterative DFS */
    static int32_t *st = 0; static i64 sc = 0;
    i64 sp = 0;
    if (!st) { sc = 1 << 16; st = malloc(sc * 4); }
    st[sp++] = x;
    while (sp) {
        int32_t y = st[--sp];
        if (ch0[y] < 0) { if (*bn == *bc) { *bc = *bc ? 2 * *bc : 1 << 20; *buf = realloc(*buf, *bc * 4); } (*buf)[(*bn)++] = lid[y]; }
        else { if (sp + 2 > sc) { sc *= 2; st = realloc(st, sc * 4); } st[sp++] = ch0[y]; st[sp++] = ch1[y]; }
    }
}

int main(void) {
    if (scanf("%lld", (long long *)&Q) != 1) return 2;
    int32_t root = newnode(0, 0, 0);
    (void)root;
    unsigned long long c; int d, s;
    while (scanf("%llu %d %d", &c, &d, &s) == 3) {
        if (d < 1 || d > 62 || !(c & 1)) { fprintf(stderr, "bad leaf\n"); return 2; }
        int32_t x = 0;
        while (dep[x] < d) {
            if (ch0[x] < 0) { int32_t a = newnode(res[x], dep[x] + 1, 0), b = newnode(res[x] | (1ULL << dep[x]), dep[x] + 1, 0); ch0[x] = a; ch1[x] = b; }
            x = ((c >> dep[x]) & 1) ? ch1[x] : ch0[x];
        }
        if (ch0[x] >= 0) { fprintf(stderr, "overlapping leaves\n"); return 2; }
        sg[x] = (int8_t)s;
    }
    /* every odd leaf must have a sign; even leaves have sign 0 */
    for (int32_t x = 0; x < nn; x++) if (ch0[x] < 0) {
        if ((res[x] & 1) && dep[x] >= 1 && sg[x] == 0) { fprintf(stderr, "uncovered odd class %llu mod 2^%d\n", (unsigned long long)res[x], dep[x]); return 2; }
        if (!(res[x] & 1)) sg[x] = 0;
        wpush(x);
    }
    while (wn) {
        int32_t x = wl[--wn];
        if (ch0[x] >= 0 || dep[x] == 0) continue;
        ensure_union(image_res(x), dep[x] - 1);
    }
    lid = malloc((size_t)nn * 4);
    nl = 0;
    for (int32_t x = 0; x < nn; x++) lid[x] = (ch0[x] < 0) ? nl++ : -1;
    leaves = malloc((size_t)nl * 4);
    for (int32_t x = 0; x < nn; x++) if (lid[x] >= 0) leaves[lid[x]] = x;
    eoff = malloc(((size_t)nl + 1) * 4);
    int32_t *buf = 0; i64 bn = 0, bc = 0;
    for (int32_t i = 0; i < nl; i++) {
        int32_t x = leaves[i];
        eoff[i] = (int32_t)bn;
        int32_t y = descend(image_res(x), dep[x] - 1);
        if (dep[y] != dep[x] - 1 && ch0[y] >= 0) { fprintf(stderr, "internal: image not resolved\n"); return 3; }
        if (ch0[y] < 0 && dep[y] < dep[x] - 1) { fprintf(stderr, "internal: image inside a coarser leaf\n"); return 3; }
        collect_under(y, &buf, &bn, &bc);
    }
    eoff[nl] = (int32_t)bn; ne = bn; edst = buf;
    printf("LEAVES %d EDGES %lld\n", nl, (long long)ne);
    int8_t *w = malloc(nl);
    for (int32_t i = 0; i < nl; i++) w[i] = (int8_t)(res[leaves[i]] & 1);

    /* Howard */
    int32_t *pol = malloc((size_t)nl * 4); double *eta = malloc((size_t)nl * 8), *h = malloc((size_t)nl * 8);
    int8_t *state = malloc(nl); int32_t *path = malloc((size_t)nl * 4);
    for (int32_t i = 0; i < nl; i++) { pol[i] = edst[eoff[i]]; for (int32_t j = eoff[i]; j < eoff[i + 1]; j++) if (w[edst[j]]) { pol[i] = edst[j]; break; } }
    for (int it = 0; it < 100000; it++) {
        memset(state, 0, nl);
        for (int32_t s0 = 0; s0 < nl; s0++) {
            if (state[s0]) continue;
            i64 pn = 0; int32_t u = s0;
            while (!state[u]) { state[u] = 1; path[pn++] = u; u = pol[u]; }
            if (state[u] == 1) {
                i64 j = 0; while (path[j] != u) j++;
                i64 L = pn - j, a = 0; for (i64 t = j; t < pn; t++) a += w[path[t]];
                double m = (double)a / (double)L;
                h[u] = 0; eta[u] = m;
                for (i64 t = pn - 1; t > j; t--) { int32_t v = path[t]; eta[v] = m; h[v] = w[v] - m + h[pol[v]]; }
                for (i64 t = j; t < pn; t++) state[path[t]] = 2;
                pn = j;
            }
            for (i64 t = pn - 1; t >= 0; t--) { int32_t v = path[t]; eta[v] = eta[pol[v]]; h[v] = w[v] - eta[v] + h[pol[v]]; state[v] = 2; }
        }
        int changed = 0; const double eps = 1e-10;
        for (int32_t v = 0; v < nl; v++) {
            double be = eta[pol[v]]; int32_t bt = pol[v];
            for (int32_t j = eoff[v]; j < eoff[v + 1]; j++) { int32_t t = edst[j]; if (eta[t] > be + eps) { be = eta[t]; bt = t; } }
            if (bt != pol[v]) { pol[v] = bt; changed = 1; }
        }
        if (changed) continue;
        for (int32_t v = 0; v < nl; v++) {
            double cur = h[pol[v]]; int32_t bt = pol[v];
            for (int32_t j = eoff[v]; j < eoff[v + 1]; j++) { int32_t t = edst[j]; if (fabs(eta[t] - eta[v]) < eps && h[t] > cur + eps) { cur = h[t]; bt = t; } }
            if (bt != pol[v]) { pol[v] = bt; changed = 1; }
        }
        if (!changed) break;
    }
    /* best cycle of the final policy */
    i64 ba = -1, bp = 1; int32_t bstart = -1;
    memset(state, 0, nl);
    for (int32_t s0 = 0; s0 < nl; s0++) {
        if (state[s0]) continue;
        i64 pn = 0; int32_t u = s0;
        while (!state[u]) { state[u] = 1; path[pn++] = u; u = pol[u]; }
        if (state[u] == 1) {
            i64 j = 0; while (path[j] != u) j++;
            i64 L = pn - j, a = 0; for (i64 t = j; t < pn; t++) a += w[path[t]];
            if (a * bp > ba * L) { ba = a; bp = L; bstart = u; }
        }
        for (i64 t = 0; t < pn; t++) state[path[t]] = 2;
    }
    /* exact certification with Dinkelbach fallback */
    i64 *psi = malloc((size_t)nl * 8); int32_t *ptr = malloc((size_t)nl * 4);
    int32_t *roff = calloc((size_t)nl + 1, 4), *rdst = malloc((size_t)ne * 4);
    for (i64 j = 0; j < ne; j++) roff[edst[j] + 1]++;
    for (int32_t i = 0; i < nl; i++) roff[i + 1] += roff[i];
    { int32_t *fill = calloc(nl, 4); for (int32_t v = 0; v < nl; v++) for (int32_t j = eoff[v]; j < eoff[v + 1]; j++) { int32_t t = edst[j]; rdst[roff[t] + fill[t]++] = v; } free(fill); }
    int32_t *q = malloc(((size_t)nl + 1) * 4); int8_t *inq = malloc(nl);
    i64 A = ba, P = bp;
    { i64 g = A, b = P; while (b) { i64 t = g % b; g = b; b = t; } A /= g; P /= g; }
    int32_t *wit = 0; i64 wlen = 0;
    /* witness from Howard */
    { wit = malloc((size_t)bp * 4); int32_t u = bstart; for (i64 t = 0; t < bp; t++) { wit[t] = u; u = pol[u]; } wlen = bp; }
    int rounds = 0;
    const i64 CHECK = 1 << 21;
    while (1) {                        /* one Dinkelbach round per threshold A/P */
        rounds++;
        if (rounds > 10000) { fprintf(stderr, "no convergence\n"); return 3; }
        for (int32_t i = 0; i < nl; i++) { psi[i] = 0; ptr[i] = -1; inq[i] = 1; q[i] = i; }
        i64 qh = 0, qt = nl, qn = (i64)nl + 1, pops = 0;
        int restart = 0;
        while (qh != qt) {
            int32_t v = q[qh]; qh = (qh + 1) % qn; inq[v] = 0; pops++;
            i64 e = w[v] ? P - A : -A;
            i64 best = -1; int32_t arg = -1;
            for (int32_t j = eoff[v]; j < eoff[v + 1]; j++) { int32_t t = edst[j]; if (psi[t] > best) { best = psi[t]; arg = t; } }
            i64 val = e + best;
            if (val > psi[v]) {
                psi[v] = val; ptr[v] = arg;
                for (int32_t j = roff[v]; j < roff[v + 1]; j++) { int32_t r = rdst[j]; if (!inq[r]) { inq[r] = 1; q[qt] = r; qt = (qt + 1) % qn; } }
            }
            if (pops % CHECK == 0) {   /* pointer-cycle check: a pointer cycle is denser than A/P */
                memset(state, 0, nl);
                for (int32_t s0 = 0; s0 < nl && !restart; s0++) {
                    if (state[s0] || ptr[s0] < 0) continue;
                    i64 pn = 0; int32_t u = s0;
                    while (u >= 0 && !state[u]) { state[u] = 1; path[pn++] = u; u = ptr[u]; }
                    if (u >= 0 && state[u] == 1) {
                        i64 j = 0; while (path[j] != u) j++;
                        i64 L = pn - j, a = 0; for (i64 t = j; t < pn; t++) a += w[path[t]];
                        if (!(a * P > A * L)) { fprintf(stderr, "internal: pointer cycle not denser\n"); return 3; }
                        free(wit); wit = malloc((size_t)L * 4); for (i64 t = 0; t < L; t++) wit[t] = path[j + t]; wlen = L;
                        A = a; P = L; { i64 g = A, b = P; while (b) { i64 t = g % b; g = b; b = t; } A /= g; P /= g; }
                        restart = 1;
                    }
                    for (i64 t = 0; t < pn; t++) state[path[t]] = 2;
                }
                if (restart) break;
            }
        }
        if (!restart) break;
    }
    /* exact check of the potential on every edge */
    i64 mx = 0;
    for (int32_t v = 0; v < nl; v++) {
        i64 e = w[v] ? P - A : -A;
        if (psi[v] < 0) { printf("CERT FAIL\n"); return 1; }
        if (psi[v] > mx) mx = psi[v];
        for (int32_t j = eoff[v]; j < eoff[v + 1]; j++) if (psi[edst[j]] + e > psi[v]) { printf("CERT FAIL\n"); return 1; }
    }
    /* witness density must equal A/P */
    { i64 a = 0; for (i64 t = 0; t < wlen; t++) a += w[wit[t]]; if (a * P != A * wlen) { printf("WITNESS MISMATCH\n"); return 1; } }
    /* critical set: nodes in nontrivial SCCs of the tight subgraph (edges with psi(t) + e(v) == psi(v)); these are
       exactly the nodes on cycles of density A/P.  Iterative Tarjan. */
    {
        int32_t *idx = malloc((size_t)nl * 4), *low = malloc((size_t)nl * 4), *stk = malloc((size_t)nl * 4);
        int32_t *cs = malloc((size_t)nl * 4), *ci = malloc((size_t)nl * 4);
        int8_t *onst = calloc(nl, 1);
        for (int32_t i = 0; i < nl; i++) idx[i] = -1;
        int32_t counter = 0, sp = 0; i64 ncrit = 0; double mcrit = 0;
        for (int32_t r0 = 0; r0 < nl; r0++) {
            if (idx[r0] >= 0) continue;
            int32_t csp = 0;
            cs[csp] = r0; ci[csp] = eoff[r0]; csp++;
            idx[r0] = low[r0] = counter++; stk[sp++] = r0; onst[r0] = 1;
            while (csp) {
                int32_t v = cs[csp - 1];
                i64 e = w[v] ? P - A : -A;
                int pushed = 0;
                while (ci[csp - 1] < eoff[v + 1]) {
                    int32_t t = edst[ci[csp - 1]++];
                    if (psi[t] + e != psi[v]) continue;
                    if (idx[t] < 0) { idx[t] = low[t] = counter++; stk[sp++] = t; onst[t] = 1; cs[csp] = t; ci[csp] = eoff[t]; csp++; pushed = 1; break; }
                    else if (onst[t] && idx[t] < low[v]) low[v] = idx[t];
                }
                if (pushed) continue;
                if (low[v] == idx[v]) {
                    int32_t cnt = 0, top = sp;
                    while (1) { int32_t x = stk[--sp]; onst[x] = 0; cnt++; if (x == v) break; }
                    int self = 0;
                    if (cnt == 1) { for (int32_t j = eoff[v]; j < eoff[v + 1]; j++) if (edst[j] == v && psi[v] + e == psi[v]) self = 1; }
                    if (cnt > 1 || self) { ncrit += cnt; for (int32_t z = sp; z < top; z++) mcrit += ldexp(1.0, -dep[leaves[stk[z]]]); }
                }
                csp--;
                if (csp) { int32_t u = cs[csp - 1]; if (low[v] < low[u]) low[u] = low[v]; }
            }
        }
        printf("CRIT %lld %.6g\n", (long long)ncrit, mcrit);
    }
    printf("RHO %lld %lld\n", (long long)A, (long long)P);
    printf("CERT ok %lld rounds %d\n", (long long)mx, rounds);
    printf("CYCLE");
    for (i64 t = 0; t < wlen; t++) printf(" %llu %d %d", (unsigned long long)res[leaves[wit[t]]], dep[leaves[wit[t]]], sg[leaves[wit[t]]]);
    printf("\n");
    return 0;
}
