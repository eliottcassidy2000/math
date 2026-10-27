/* endpoint_checks.c -- mac-mini-2026-10-27-dipole
 * The endpoint principle: for score and cut potentials the value after the move along a Hamiltonian path Q depends
 * only on the endpoints (a, b) of Q (dipole law d -> -d + e_a - e_b).  So T has a descending move iff the endpoint
 * relation Pi(T) = {(a,b): some Hamiltonian path runs from a to b} meets the descending set D_Phi(T).  This program
 * computes Pi(T) by a subset DP for every class read from gentourng and checks, without any move graph:
 *   (C1) non-strong T with components C_1 => ... => C_k:  Pi(T) = C_1 x C_k;
 *   (C2) delta0 (target: a source or sink):  T with delta0 >= 1 has a descending move iff T is strong;
 *   (C3) L (target: transitive):  a non-strong non-transitive T has a descending move iff it has no source and no sink;
 *   (C4) R (odd n, target: regular):  every non-regular T has a descending move  [conjecture tested];
 *   (C5) strong L-local minima, by (lambda, kappa) and by score sequence;  whether a Hamiltonian path joins some
 *        minimum-score vertex to some maximum-score vertex;
 *   (C6) Thomassen-type statement: kappa >= 2 implies every pair is joined by a Hamiltonian path in some direction.
 * Build: cc -O3 -march=native endpoint_checks.c -o endpoint_checks ;  Run: gentourng -q n | ./endpoint_checks n
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n;
static uint32_t out_[16];
static uint16_t pi_[16];          /* pi_[a] = bitmask of ends b of Hamiltonian paths starting at a */
static uint16_t dp[1 << 13][13];  /* dp[S][v] = bitmask of starts a of paths covering S ending at v */
static void endpoints(void) {
    int full = (1 << n) - 1;
    for (int S = 0; S <= full; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1 << v;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) { uint16_t st = dp[S][v]; if (!st) continue;
        uint32_t nx = out_[v] & ~S; while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dp[S | 1 << w][w] |= st; } }
    for (int a = 0; a < n; a++) pi_[a] = 0;
    for (int b = 0; b < n; b++) { uint16_t st = dp[full][b]; while (st) { int a = __builtin_ctz(st); st &= st - 1; pi_[a] |= 1 << b; } }
}
static int strongMask(uint32_t vs) {
    if (!vs) return 1; int start = __builtin_ctz(vs);
    uint32_t seen = 1u << start, fr = seen;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= out_[__builtin_ctz(f)] & vs; nx &= ~seen; seen |= nx; fr = nx; }
    if (seen != vs) return 0;
    seen = 1u << start; fr = seen;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) { int v = __builtin_ctz(f); for (uint32_t w = vs; w; w &= w - 1) { int u = __builtin_ctz(w); if (out_[u] >> v & 1) nx |= 1u << u; } } nx &= ~seen; seen |= nx; fr = nx; }
    return seen == vs;
}
static int Lof(const int *s) {   /* (1/2) sum |s_(i) - i| */
    int t[16]; memcpy(t, s, n * sizeof(int));
    for (int i = 1; i < n; i++) { int v = t[i], j = i - 1; while (j >= 0 && t[j] > v) { t[j + 1] = t[j]; j--; } t[j + 1] = v; }
    int l = 0; for (int i = 0; i < n; i++) l += abs(t[i] - i); return l / 2;
}
static int d0of(const int *s) { int d = 99; for (int v = 0; v < n; v++) { int x = s[v] < n - 1 - s[v] ? s[v] : n - 1 - s[v]; if (x < d) d = x; } return d; }
static int Rof(const int *s) { int r = 0; for (int v = 0; v < n; v++) r += abs(2 * s[v] - (n - 1)); return r / 4; }   /* odd n */
typedef int (*potf)(const int *);
static int descends(potf P, const int *s) {   /* exists (a,b) in Pi with P(n-1-s + e_a - e_b) = P(s) - 1 */
    int p0 = P(s), t[16];
    for (int a = 0; a < n; a++) for (uint32_t e = pi_[a]; e; e &= e - 1) { int b = __builtin_ctz(e);
        for (int v = 0; v < n; v++) t[v] = n - 1 - s[v]; t[a]++; t[b]--;
        if (P(t) == p0 - 1) return 1; }
    return 0;
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); uint32_t full = (1u << n) - 1;
    char line[256]; long classes = 0, strong = 0;
    long c1_bad = 0, c2_bad = 0, c2_lmin = 0, c3_bad = 0, c3_lmin_ns = 0, c4_bad = 0, c4_nonreg = 0, c6_bad = 0, c6_tested = 0;
    long lmin_strong_hist[16][16] = {{0}}; long lmin_strong = 0, lmin_strong_minmax = 0, strong_minmax_fail = 0;
    long c5_by_kappa[16] = {0}; long lmin_kappa2plus = 0, c7_bad = 0, c7_n = 0, c8_bad = 0, c9_bad = 0, c7a_bad = 0, c7b_bad = 0;
    long lu_minmax_pairs = 0, lu_minmax_bad = 0, lu_maxmin_pairs = 0, lu_maxmin_bad = 0, lu_maxmin_bad_classes = 0, lmin_with_minmax_arc = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        classes++;
        for (int i = 0; i < n; i++) out_[i] = 0;
        int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') out_[i] |= 1u << j; else out_[j] |= 1u << i; k++; }
        int s[16]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(out_[v]);
        endpoints();
        /* strong components via scores (Landau): sort by score, cut where partial sums are binomial */
        int ord[16]; for (int i = 0; i < n; i++) ord[i] = i;
        for (int i = 1; i < n; i++) { int v = ord[i], j = i - 1; while (j >= 0 && s[ord[j]] > s[v]) { ord[j + 1] = ord[j]; j--; } ord[j + 1] = v; }
        uint32_t comp[16]; int nc = 0, cum = 0; uint32_t cur = 0;
        for (int i = 0; i < n; i++) { cum += s[ord[i]]; cur |= 1u << ord[i]; if (cum == (i + 1) * i / 2) { comp[nc++] = cur; cur = 0; } }
        /* comp[0] = bottom component (lowest scores), comp[nc-1] = top */
        int isStrong = (nc == 1);
        int lam = 0, kap = 0;
        if (isStrong) {
            strong++;
            lam = 99; for (uint32_t X = 1; X < full; X++) { int f = 0; for (int x = 0; x < n; x++) if (X >> x & 1) f += __builtin_popcount(out_[x] & ~X & full); if (f < lam) lam = f; }
            kap = n - 1; for (uint32_t S = 1; S < full; S++) { int sz = __builtin_popcount(S); if (sz >= kap || sz > n - 2) continue; if (!strongMask(full & ~S)) kap = sz; }
        }
        /* (C1) */
        if (!isStrong) {
            uint32_t top = comp[nc - 1], bot = comp[0];
            for (int a = 0; a < n; a++) { uint32_t expect = (top >> a & 1) ? bot : 0; if (pi_[a] != expect) { c1_bad++; break; } }
        }
        /* (C2) delta0 */
        int d0 = d0of(s);
        if (d0 >= 1) { int de = descends(d0of, s); if (de != isStrong) c2_bad++; if (!de) c2_lmin++; }
        /* (C3) L */
        int L = Lof(s);
        if (L > 0) {
            int de = descends(Lof, s);
            if (!isStrong) { int srcsink = 0; for (int v = 0; v < n; v++) if (s[v] == 0 || s[v] == n - 1) srcsink = 1; if (de == srcsink) c3_bad++; if (!de) c3_lmin_ns++; }
            else if (!de) { lmin_strong++; lmin_strong_hist[lam][kap]++; if (kap >= 2) lmin_kappa2plus++; }
        }
        /* (C5) strong: Hamiltonian path from some min-score vertex to some max-score vertex */
        if (isStrong) {
            int mn = 99, mx = -1; for (int v = 0; v < n; v++) { if (s[v] < mn) mn = s[v]; if (s[v] > mx) mx = s[v]; }
            uint32_t MX = 0; for (int v = 0; v < n; v++) if (s[v] == mx) MX |= 1u << v;
            int ok = 0; for (int a = 0; a < n; a++) if (s[a] == mn && (pi_[a] & MX)) ok = 1;
            if (!ok && mn != mx) { strong_minmax_fail++; c5_by_kappa[kap]++; if (L > 0 && !descends(Lof, s)) lmin_strong_minmax++; }
        }
        /* (C7) strong nonregular: a Hamiltonian path from some max-score vertex to some min-score vertex;
           (C8) from EVERY max-score vertex to some vertex of score < (n-1)/2;  (C9) from some vertex of score > (n-1)/2
           to some vertex of score < (n-1)/2 (any n; for odd n this is the R-descent) */
        if (isStrong) {
            int mn = 99, mx = -1; for (int v = 0; v < n; v++) { if (s[v] < mn) mn = s[v]; if (s[v] > mx) mx = s[v]; }
            if (mn != mx) {
                c7_n++;
                uint32_t MN = 0, LO = 0, HI = 0; for (int v = 0; v < n; v++) { if (s[v] == mn) MN |= 1u << v; if (2 * s[v] < n - 1) LO |= 1u << v; if (2 * s[v] > n - 1) HI |= 1u << v; }
                int ok7 = 0, ok8 = 1, ok9 = 0;
                for (int a = 0; a < n; a++) { if (s[a] == mx && (pi_[a] & MN)) ok7 = 1; if (s[a] == mx && !(pi_[a] & LO)) ok8 = 0; if ((HI >> a & 1) && (pi_[a] & LO)) ok9 = 1; }
                if (!ok7) c7_bad++; if (!ok8 && LO) c8_bad++; if (!ok9 && LO && HI) c9_bad++;
                int ok7a = 1, ok7b = 1;   /* (C7a) every max-score vertex reaches some min-score vertex; (C7b) every max to every min */
                for (int a = 0; a < n; a++) if (s[a] == mx) { if (!(pi_[a] & MN)) ok7a = 0; if ((pi_[a] & MN) != MN) ok7b = 0; }
                if (!ok7a) c7a_bad++; if (!ok7b) c7b_bad++;
                /* Lu (2018), as read from the abstract: for a min-max pair joined by an arc u -> v, a Hamiltonian path runs from u to v
                   (min-max arc: always; max-min arc: with a few exceptions).  Checked over all pairs of strong nonregular classes. */
                int badc = 0, hasMinMaxArc = 0;
                for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (u != v && (out_[u] >> v & 1)) {
                    if (s[u] == mn && s[v] == mx) { lu_minmax_pairs++; hasMinMaxArc = 1; if (!(pi_[v] >> u & 1)) lu_minmax_bad++; }
                    if (s[u] == mx && s[v] == mn) { lu_maxmin_pairs++; if (!(pi_[u] >> v & 1)) { lu_maxmin_bad++; badc = 1; } } }
                lu_maxmin_bad_classes += badc;
                if (hasMinMaxArc && L > 0 && !descends(Lof, s)) lmin_with_minmax_arc++;
            }
        }
        /* (C4) R, odd n */
        if (n % 2 == 1) { int R = Rof(s); if (R > 0) { c4_nonreg++; if (!descends(Rof, s)) c4_bad++; } }
        /* (C6) */
        if (isStrong && kap >= 2) { c6_tested++; for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (!((pi_[a] >> b & 1) || (pi_[b] >> a & 1))) { c6_bad++; a = n; break; } }
    }
    printf("n=%d classes=%ld strong=%ld\n", n, classes, strong);
    printf("(C1) non-strong: Pi(T) = C_1 x C_k fails for %ld classes\n", c1_bad);
    printf("(C2) delta0 >= 1: 'descending move <=> strong' fails for %ld; delta0 local minima %ld\n", c2_bad, c2_lmin);
    printf("(C3) L, non-strong non-transitive: 'descending move <=> no source and no sink' fails for %ld; non-strong L-local minima %ld\n", c3_bad, c3_lmin_ns);
    printf("(C3') strong L-local minima %ld, by (lambda,kappa):", lmin_strong);
    for (int l = 0; l < 16; l++) for (int q = 0; q < 16; q++) if (lmin_strong_hist[l][q]) printf(" (%d,%d):%ld", l, q, lmin_strong_hist[l][q]);
    printf("; with kappa >= 2: %ld\n", lmin_kappa2plus);
    printf("(C5) strong nonregular classes with no Hamiltonian path from a min-score to a max-score vertex: %ld (of which L-local minima %ld); by kappa:", strong_minmax_fail, lmin_strong_minmax);
    for (int q = 0; q < 16; q++) if (c5_by_kappa[q]) printf(" %d:%ld", q, c5_by_kappa[q]); printf("\n");
    if (n % 2 == 1) printf("(C4) R (regular target): non-regular classes %ld, without a descending move %ld\n", c4_nonreg, c4_bad);
    printf("(C6) kappa >= 2: classes %ld, with a pair joined in neither direction %ld\n", c6_tested, c6_bad);
    printf("(C7) strong nonregular %ld: no Hamiltonian path from a max-score to a min-score vertex: %ld\n", c7_n, c7_bad);
    printf("(C7a) some max-score vertex reaches no min-score vertex: %ld;  (C7b) some (max, min) pair not joined max -> min: %ld\n", c7a_bad, c7b_bad);
    printf("(Lu) min->max arcs y->x %ld, without a Hamiltonian path x -> y (max to min) %ld; max->min arcs %ld, not extendable %ld (in %ld classes)\n", lu_minmax_pairs, lu_minmax_bad, lu_maxmin_pairs, lu_maxmin_bad, lu_maxmin_bad_classes);
    printf("(Lu') strong L-local minima having a min->max arc: %ld\n", lmin_with_minmax_arc);
    printf("(C8) some max-score vertex starts no Hamiltonian path ending below (n-1)/2: %ld\n", c8_bad);
    printf("(C9) no Hamiltonian path from a vertex above (n-1)/2 to a vertex below it: %ld\n", c9_bad);
    return 0;
}
