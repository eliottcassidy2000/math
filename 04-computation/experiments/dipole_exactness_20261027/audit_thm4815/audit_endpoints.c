/* audit_endpoints.c -- independent audit of THM-4815 (B3, B4, B5, B7, B8 descent, Lu-type arc counts).
 * Written from scratch for the audit; does not reuse producer code.
 * Input: gentourng ascii (upper triangle row by row, '1' at pair (i,j), i<j, means i->j).
 * Pi(T) = {(a,b): some Hamiltonian path runs a -> ... -> b}, computed by a subset DP:
 *   E[S][v] = set of starts a such that some path with vertex set S starts at a and ends at v.
 * Checks per class (all classes, strong or not):
 *   B7 : some x in X (max score), y in Y (min score) with (x,y) in Pi.
 *   B7(i)  : every pair x in X, y in Y with y -> x has (x,y) in Pi.
 *   B7(ii) : if X => Y (every max vertex beats every min vertex), every pair (x,y) in X x Y is in Pi.
 *   Lu arcs: x -> y with (x,y) not in Pi; classes; tie status; whether X => Y fails there.
 *   B3 : non-strong T: Pi = C_1 x C_k exactly; X in C_1; Y in C_k; min score on C_1 > max score on C_k.
 *   B4 : delta0 = k >= 1: (some (a,b) in Pi with delta0(new) = k-1) <=> strong.  Lipschitz |change| <= 1 on Pi.
 *   B5 : non-strong non-transitive: L-descent <=> no source and no sink.  Strong: (min, max) pairs lie in D_L.
 *   strong L-local minima: count, (lambda, kappa) histogram, kappa >= 2 count.
 *   B8 : odd n, non-regular: some (a,b) in Pi with R(new) = R - 1; and D_R = {s_a > m, s_b < m} exactly.
 *   note sec 2: strong non-regular with no Hamiltonian path from a min-score to a max-score vertex, by kappa;
 *               the 2-connected ones: do they still have an L-descent?
 * New scores after the move along a path from a to b: s'_v = n-1-s_v + [v=a] - [v=b] (computed explicitly).
 * Build: cc -O3 -march=native audit_endpoints.c -o audit_endpoints ; Run: gentourng -q n | ./audit_endpoints n
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n;
static uint32_t out[32], full;
static uint32_t *E;     /* E[S*n+v] */
static uint32_t Pi[32]; /* Pi[a] = set of ends b */

static void endpoint_relation(void) {
    size_t N = (size_t)1 << n;
    memset(E, 0, sizeof(uint32_t) * N * n);
    for (int v = 0; v < n; v++) E[((size_t)1 << v) * n + v] = 1u << v;
    for (size_t S = 1; S < N; S++) {
        for (int v = 0; v < n; v++) {
            uint32_t st = E[S * n + v];
            if (!st) continue;
            uint32_t nx = out[v] & ~(uint32_t)S & full;
            while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; E[(S | ((size_t)1 << w)) * n + w] |= st; }
        }
    }
    for (int a = 0; a < n; a++) Pi[a] = 0;
    for (int b = 0; b < n; b++) { uint32_t st = E[(N - 1) * n + b]; while (st) { int a = __builtin_ctz(st); st &= st - 1; Pi[a] |= 1u << b; } }
}

/* reachability inside vertex set vs from vertex v */
static uint32_t reach_in(uint32_t vs, int v) {
    uint32_t seen = 1u << v, fr = seen;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= out[__builtin_ctz(f)] & vs; nx &= ~seen; seen |= nx; fr = nx; }
    return seen;
}
static int strong_in(uint32_t vs) {   /* induced subtournament on vs strongly connected (|vs| = 1 counts as strong) */
    if (!vs) return 1;
    int pc = __builtin_popcount(vs); if (pc == 1) return 1; if (pc == 2) return 0;
    for (uint32_t w = vs; w; w &= w - 1) if (reach_in(vs, __builtin_ctz(w)) != vs) return 0;
    return 1;
}
static int kappa_of(void) {    /* least |S| with T - S not strong (|T-S| >= 2); T assumed strong */
    for (int k = 1; k <= n - 2; k++) {
        for (uint32_t S = 1; S < full; S++) if (__builtin_popcount(S) == k) if (!strong_in(full & ~S)) return k;
    }
    return n - 1;
}
static int lambda_of(void) {
    int best = 1 << 20;
    for (uint32_t X = 1; X < full; X++) { int f = 0; for (int x = 0; x < n; x++) if (X >> x & 1) f += __builtin_popcount(out[x] & ~X & full); if (f < best) best = f; }
    return best;
}
static int cmpint(const void *a, const void *b) { return *(const int *)a - *(const int *)b; }
static int Lpot2(const int *s) { int t[32]; memcpy(t, s, sizeof(int) * n); qsort(t, n, sizeof(int), cmpint); int l = 0; for (int i = 0; i < n; i++) l += abs(t[i] - i); return l; } /* = 2L */
static int d0pot(const int *s) { int d = 1 << 20; for (int v = 0; v < n; v++) { int x = s[v] < n - 1 - s[v] ? s[v] : n - 1 - s[v]; if (x < d) d = x; } return d; }
static int Rpot4(const int *s) { int r = 0; for (int v = 0; v < n; v++) r += abs(2 * s[v] - (n - 1)); return r; } /* = 4R at odd n */
static void newscores(const int *s, int a, int b, int *t) { for (int v = 0; v < n; v++) t[v] = n - 1 - s[v]; t[a]++; t[b]--; }

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n\n", argv[0]); return 1; }
    n = atoi(argv[1]); full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    E = malloc(sizeof(uint32_t) * ((size_t)1 << n) * n);
    char line[1024];
    long classes = 0, strongc = 0, regular = 0;
    long b7_fail = 0, i_arcs = 0, i_fail = 0, i_classes = 0, ii_classes = 0, ii_fail = 0;
    long i_arcs_strong = 0;
    long lu_arcs = 0, lu_classes = 0, lu_untied = 0, lu_XY = 0, lu_arcs_strong = 0, lu_classes_strong = 0;
    long b3_fail = 0, b3_score_fail = 0, b4_fail = 0, b4_lmin = 0, lip_fail = 0, b5_fail = 0, b5_lmin = 0;
    long b5_strong_minmax_fail = 0, sLmin = 0, sLmin_k2 = 0, hist[40][40]; memset(hist, 0, sizeof hist);
    long b8_nonreg = 0, b8_fail = 0, b8_char_fail = 0;
    long mm_fail = 0, mm_fail_k2 = 0, mm_fail_k2_nodesc = 0, mm_fail_k1 = 0;
    long nonstrong_nontrans = 0, ns_srcsink = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        classes++;
        for (int i = 0; i < n; i++) out[i] = 0;
        int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') out[i] |= 1u << j; else out[j] |= 1u << i; k++; }
        int s[32], Dl = -1, dl = 1 << 20;
        for (int v = 0; v < n; v++) { s[v] = __builtin_popcount(out[v]); if (s[v] > Dl) Dl = s[v]; if (s[v] < dl) dl = s[v]; }
        uint32_t X = 0, Y = 0;
        for (int v = 0; v < n; v++) { if (s[v] == Dl) X |= 1u << v; if (s[v] == dl) Y |= 1u << v; }
        endpoint_relation();
        /* strong components: top = vertices reaching all; bottom = vertices reached by all */
        uint32_t R[32], top = 0, bot = full;
        for (int v = 0; v < n; v++) { R[v] = reach_in(full, v); if (R[v] == full) top |= 1u << v; bot &= R[v]; }
        int strong = (top == full);
        if (strong) strongc++;
        int srcsink = 0; for (int v = 0; v < n; v++) if (s[v] == 0 || s[v] == n - 1) srcsink = 1;
        /* ---- B7 ---- */
        if (Dl == dl) { regular++; }
        else {
            int ok = 0; for (uint32_t x = X; x; x &= x - 1) if (Pi[__builtin_ctz(x)] & Y) ok = 1;
            if (!ok) { b7_fail++; printf("B7 FAIL: %s", line); }
            int anyyx = 0, clsfail = 0, clslu = 0;
            for (uint32_t xs = X; xs; xs &= xs - 1) { int x = __builtin_ctz(xs);
                for (uint32_t ys = Y; ys; ys &= ys - 1) { int y = __builtin_ctz(ys);
                    if (out[y] >> x & 1) { anyyx = 1; i_arcs++; if (strong) i_arcs_strong++; if (!(Pi[x] >> y & 1)) { i_fail++; clsfail = 1; } }
                    else if (!(Pi[x] >> y & 1)) { lu_arcs++; if (strong) lu_arcs_strong++; clslu = 1; }
                } }
            if (anyyx) i_classes++;
            if (clsfail) printf("B7(i) FAIL: %s", line);
            if (!anyyx) { ii_classes++; int f = 0; for (uint32_t xs = X; xs; xs &= xs - 1) if ((Pi[__builtin_ctz(xs)] & Y) != Y) f = 1; if (f) { ii_fail++; printf("B7(ii) FAIL: %s", line); } }
            if (clslu) { lu_classes++; if (strong) lu_classes_strong++; if (__builtin_popcount(X) == 1 && __builtin_popcount(Y) == 1) lu_untied++; if (!anyyx) lu_XY++; }
        }
        /* ---- B3 ---- */
        if (!strong) {
            for (int a = 0; a < n; a++) { uint32_t expect = (top >> a & 1) ? bot : 0; if (Pi[a] != expect) { b3_fail++; printf("B3 FAIL: %s", line); break; } }
            int mnTop = 1 << 20, mxBot = -1; for (int v = 0; v < n; v++) { if ((top >> v & 1) && s[v] < mnTop) mnTop = s[v]; if ((bot >> v & 1) && s[v] > mxBot) mxBot = s[v]; }
            if ((X & ~top) || (Y & ~bot) || !(mnTop > mxBot)) { b3_score_fail++; printf("B3 score FAIL: %s", line); }
        }
        /* ---- potentials on Pi ---- */
        int L2 = Lpot2(s), d0 = d0pot(s), R4 = Rpot4(s);
        int descL = 0, descD = 0, descR = 0, t[32];
        for (int a = 0; a < n; a++) for (uint32_t e = Pi[a]; e; e &= e - 1) { int b = __builtin_ctz(e);
            newscores(s, a, b, t);
            int l2 = Lpot2(t), dd = d0pot(t), r4 = Rpot4(t);
            if (abs(l2 - L2) > 2 || abs(dd - d0) > 1 || ((n & 1) && abs(r4 - R4) > 4)) lip_fail++;
            if (l2 == L2 - 2) descL = 1;
            if (dd == d0 - 1) descD = 1;
            if ((n & 1) && r4 == R4 - 4) descR = 1;
            if (n & 1) { int m = (n - 1) / 2; int pred = (s[a] > m) && (s[b] < m); if (pred != (r4 == R4 - 4)) b8_char_fail++; }
        }
        /* ---- B4 ---- */
        if (d0 >= 1) { if (descD != strong) { b4_fail++; printf("B4 FAIL: %s", line); } if (!descD) b4_lmin++; }
        /* ---- B5 ---- */
        if (L2 > 0 && !strong) {
            nonstrong_nontrans++; if (srcsink) ns_srcsink++;
            if (descL == srcsink) { b5_fail++; printf("B5 FAIL: %s", line); }
            if (!descL) b5_lmin++;
        }
        if (strong && n >= 3) {   /* (min, max) pairs in D_L */
            for (uint32_t ys = Y; ys; ys &= ys - 1) for (uint32_t xs = X; xs; xs &= xs - 1) { int a = __builtin_ctz(ys), b = __builtin_ctz(xs);
                if (a == b) continue; newscores(s, a, b, t); if (Lpot2(t) != L2 - 2) b5_strong_minmax_fail++; }
        }
        if (strong && L2 > 0 && !descL) {
            sLmin++; int kp = kappa_of(), lm = lambda_of(); hist[lm][kp]++; if (kp >= 2) { sLmin_k2++; printf("STRONG 2-CONNECTED L-LOCAL-MIN: %s", line); }
        }
        /* ---- B8 ---- */
        if ((n & 1) && R4 > 0) { b8_nonreg++; if (!descR) { b8_fail++; printf("B8 FAIL: %s", line); } }
        /* ---- note sec 2: min -> max Hamiltonian paths in strong non-regular classes ---- */
        if (strong && Dl != dl) {
            int ok = 0; for (uint32_t ys = Y; ys; ys &= ys - 1) if (Pi[__builtin_ctz(ys)] & X) ok = 1;
            if (!ok) { mm_fail++; int cut = 0; for (int v = 0; v < n; v++) if (!strong_in(full & ~(1u << v))) cut = 1;
                if (cut) mm_fail_k1++; else { mm_fail_k2++; if (!descL) mm_fail_k2_nodesc++; } }
        }
    }
    printf("n=%d classes=%ld strong=%ld regular=%ld\n", n, classes, strongc, regular);
    printf("B7: non-regular classes without a max->min Hamiltonian path: %ld\n", b7_fail);
    printf("B7(i): classes with some arc y->x (y min, x max): %ld; such arcs %ld (in strong classes %ld); arcs y->x with no HP x..y: %ld\n", i_classes, i_arcs, i_arcs_strong, i_fail);
    printf("B7(ii): classes with X => Y: %ld; of which some pair (x,y) not joined x..y: %ld\n", ii_classes, ii_fail);
    printf("Lu-type: arcs x->y (x max, y min) with no HP x..y: %ld arcs in %ld classes (strong only: %ld arcs in %ld classes); untied: %ld; with X=>Y: %ld\n", lu_arcs, lu_classes, lu_arcs_strong, lu_classes_strong, lu_untied, lu_XY);
    printf("B3: non-strong classes with Pi != C_1 x C_k: %ld; score-placement failures: %ld\n", b3_fail, b3_score_fail);
    printf("Lipschitz failures (L, delta0, R on Pi): %ld\n", lip_fail);
    printf("B4: delta0>=1 classes with (descent != strong): %ld; delta0 local minima: %ld\n", b4_fail, b4_lmin);
    printf("B5: non-strong non-transitive %ld (with source/sink %ld): (descent == source/sink) failures %ld; non-strong L-local minima %ld\n", nonstrong_nontrans, ns_srcsink, b5_fail, b5_lmin);
    printf("B5: strong classes, (min,max) pairs not in D_L: %ld\n", b5_strong_minmax_fail);
    printf("strong L-local minima: %ld; with kappa>=2: %ld; by (lambda,kappa):", sLmin, sLmin_k2);
    for (int a = 0; a < 40; a++) for (int b = 0; b < 40; b++) if (hist[a][b]) printf(" (%d,%d):%ld", a, b, hist[a][b]);
    printf("\n");
    if (n & 1) printf("B8: non-regular classes %ld, without an R-descent in Pi: %ld; D_R characterisation failures: %ld\n", b8_nonreg, b8_fail, b8_char_fail);
    printf("min->max: strong non-regular classes with no HP from a min-score to a max-score vertex: %ld (cut vertex %ld, 2-connected %ld; 2-connected without L-descent %ld)\n", mm_fail, mm_fail_k1, mm_fail_k2, mm_fail_k2_nodesc);
    return 0;
}
