/* audit_b7_random.c -- random tests of THM-4815 B7 (and its cases (i), (ii)) beyond the exhaustive orders.
 * Written from scratch for the audit.  For each sample: X = max-score vertices, Y = min-score vertices; for each x in X
 * a subset DP gives the set of ends of Hamiltonian paths starting at x.  Checks:
 *   B7     some x in X, y in Y with a Hamiltonian path x..y;
 *   B7(i)  every pair with y -> x has such a path;  B7(ii) if X => Y, every pair has such a path;
 *   B8     (odd n) some Hamiltonian path from a vertex above (n-1)/2 to a vertex below it (implied by B7).
 * Samplers: 0 uniform; 1 near-regular (rotational base, many score-preserving 3-cycle reversals, then 1..4 random arc
 * flips); 2 tied extremes (near-regular, then extra flips toward a planted set of max and min vertices); 3 near-
 * transitive (transitive + k random flips); 4 Thomassen-type plant (random tournament on V - {x}, made so that T - x
 * ends with ... -> y -> z, z -> x, x -> y, then random flips elsewhere).
 * Build: cc -O3 -march=native audit_b7_random.c -o audit_b7_random ; Run: ./audit_b7_random n samples mode seed
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n; static uint32_t out[32], full; static uint32_t *reach;
static uint64_t rs;
static uint64_t rnd(void) { rs += 0x9E3779B97F4A7C15ULL; uint64_t z = rs; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }
static void setarc(int u, int v) { out[u] |= 1u << v; out[v] &= ~(1u << u); }
static void flip(int u, int v) { if (out[u] >> v & 1) setarc(v, u); else setarc(u, v); }
static uint32_t ends_from(int x) {     /* set of b such that some Hamiltonian path runs x .. b */
    size_t N = (size_t)1 << n;
    memset(reach, 0, sizeof(uint32_t) * N);
    reach[(size_t)1 << x] = 1u << x;
    for (size_t S = 1; S < N; S++) {
        uint32_t r = reach[S]; if (!r) continue;
        for (uint32_t e = r; e; e &= e - 1) { int v = __builtin_ctz(e); uint32_t nx = out[v] & ~(uint32_t)S & full;
            for (; nx; nx &= nx - 1) { int w = __builtin_ctz(nx); reach[S | ((size_t)1 << w)] |= 1u << w; } }
    }
    return reach[N - 1];
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); long samples = atol(argv[2]); int mode = atoi(argv[3]); rs = 0x1234567ULL ^ ((uint64_t)atol(argv[4]) << 20);
    full = (1u << n) - 1; reach = malloc(sizeof(uint32_t) * ((size_t)1 << n));
    long tested = 0, regular = 0, b7f = 0, if_ = 0, iif = 0, i_cls = 0, ii_cls = 0, tied = 0, b8f = 0, lu = 0;
    for (long it = 0; it < samples; it++) {
        for (int i = 0; i < n; i++) out[i] = 0;
        if (mode == 0 || mode == 4) { for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (rnd() & 1) setarc(i, j); else setarc(j, i); } }
        if (mode == 3) { for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) setarc(i, j); int k = 2 + rnd() % (2 * n); for (int r = 0; r < k; r++) { int a = rnd() % n, b = rnd() % n; if (a != b) flip(a, b); } }
        if (mode == 1 || mode == 2) {
            int h = n / 2;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { int d = j - i; if (d <= h && !(n % 2 == 0 && d == h && i >= h)) setarc(i, j); else setarc(j, i); }
            for (int r = 0; r < 30 * n; r++) { int a = rnd() % n, b = rnd() % n, c = rnd() % n; if (a == b || b == c || a == c) continue;
                if ((out[a] >> b & 1) && (out[b] >> c & 1) && (out[c] >> a & 1)) { setarc(b, a); setarc(c, b); setarc(a, c); } }
            int f = 1 + rnd() % 4; for (int r = 0; r < f; r++) { int a = rnd() % n, b = rnd() % n; if (a != b) flip(a, b); }
            if (mode == 2) {   /* plant ties: pick sets P (to raise) and M (to lower), then push arcs from M to P */
                int p = 2 + rnd() % 3, q = 2 + rnd() % 3; int P[8], M[8]; uint32_t used = 0;
                for (int i = 0; i < p; i++) { int v; do v = rnd() % n; while (used >> v & 1); used |= 1u << v; P[i] = v; }
                for (int i = 0; i < q; i++) { int v; do v = rnd() % n; while (used >> v & 1); used |= 1u << v; M[i] = v; }
                for (int i = 0; i < p; i++) for (int j = 0; j < q; j++) if (rnd() % 3) setarc(P[i], M[j]); else setarc(M[j], P[i]);
            }
        }
        if (mode == 4) {   /* Thomassen-type plant: x -> y, y -> z, z -> x, z beats nobody else, y beats only z outside x */
            int x = rnd() % n, y, z; do y = rnd() % n; while (y == x); do z = rnd() % n; while (z == x || z == y);
            setarc(x, y); setarc(y, z); setarc(z, x);
            for (int w = 0; w < n; w++) if (w != x && w != y && w != z) { setarc(w, z); if (rnd() % 4) setarc(w, y); }
        }
        int s[32], D = -1, d = 99; for (int v = 0; v < n; v++) { s[v] = __builtin_popcount(out[v]); if (s[v] > D) D = s[v]; if (s[v] < d) d = s[v]; }
        tested++;
        if (D == d) { regular++; continue; }
        uint32_t X = 0, Y = 0; for (int v = 0; v < n; v++) { if (s[v] == D) X |= 1u << v; if (s[v] == d) Y |= 1u << v; }
        if (__builtin_popcount(X) > 1 || __builtin_popcount(Y) > 1) tied++;
        int ok = 0, anyyx = 0, fi = 0, fii = 0, lucls = 0; int m2 = n - 1;
        uint32_t HI = 0, LO = 0; for (int v = 0; v < n; v++) { if (2 * s[v] > m2) HI |= 1u << v; if (2 * s[v] < m2) LO |= 1u << v; }
        for (uint32_t xs = X; xs; xs &= xs - 1) { int x = __builtin_ctz(xs); uint32_t e = ends_from(x);
            if (e & Y) ok = 1;
            for (uint32_t ys = Y; ys; ys &= ys - 1) { int y = __builtin_ctz(ys);
                if (out[y] >> x & 1) { anyyx = 1; if (!(e >> y & 1)) fi = 1; } else if (!(e >> y & 1)) lucls = 1; } }
        if (!anyyx) { ii_cls++; if (lucls) fii = 1; } else i_cls++;
        lu += lucls;
        if (!ok) { b7f++; printf("B7 FAIL n=%d mode=%d:", n, mode); for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) putchar((out[i] >> j & 1) ? '1' : '0'); putchar('\n'); }
        if (fi) { if_++; printf("B7(i) FAIL\n"); }
        if (fii) { iif++; printf("B7(ii) FAIL\n"); }
        if ((n & 1) && !ok) b8f++;   /* (B8 descent follows from B7 at odd n) */
    }
    printf("n=%d mode=%d samples=%ld regular %ld tied-extremes %ld | B7 failures %ld | case (i) classes %ld failures %ld | case (ii) classes %ld failures %ld | Lu-type pairs present in %ld samples\n",
           n, mode, tested, regular, tied, b7f, i_cls, if_, ii_cls, iif, lu);
    return 0;
}
