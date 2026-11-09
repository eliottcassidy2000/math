/* dom_cover.c - constructive coverage of odd numbers by dominating-set counts of even (von Bromssen; RPGFD) graphs
 * (paleyeven lane, mac-mini-2026-11-09).
 *
 * Base: all connected even (von Bromssen; RPGFD) graphs from an exhaustive dom_even run (lines "n m d E|O ..."), giving (value, size).
 * Operations with proved evenness (PROOF.md, Part B, Lemmas B2-B4); every operation strictly increases the value, so
 * values are finalised in increasing order (least realisation size is exact for this generating system):
 *   (U) G u H                         d = d(G) d(H)                          size g + h
 *   (C) K_1 join (G u E_j), j >= 1,   d = 2^t + d(G),  t = g + j >= 2       size t + 1
 *   (J) E_r join (G u E_j), r >= 2,   d = (2^r - 1)(2^t - 1) + 1 + d(G)     size t + r
 *       j >= 1, excluding E_r join E_r (= K_{r,r}, t = r, G empty)
 * The empty graph on 0 vertices has d = 1 (unit of (U)); G = K_1 is the base value 1 on 1 vertex.
 *   (G) (A u E_i) join (G u E_j), i, j >= 1, d = (2^a - 1)(2^t - 1) + d(A) + d(G), size a + t  [optional, AMAX]
 * usage: dom_cover basefile LOG2B [AMAX]   -> missing odd values <= 2^LOG2B and the worst size / log2(value)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

static int GENJ = 0; static uint64_t AMAX = 255;
int main(int argc, char **argv) {
    if (argc > 3) { GENJ = 1; AMAX = strtoull(argv[3], 0, 10); }
    if (argc < 3) { fprintf(stderr, "usage: dom_cover basefile LOG2B\n"); return 1; }
    int L = atoi(argv[2]); uint64_t B = 1ULL << L;
    uint8_t *sz = malloc(B + 1); memset(sz, 255, B + 1);     /* 255 = unreached */
    sz[1] = 0;                                               /* empty graph */
    FILE *f = fopen(argv[1], "r"); char line[4096];
    while (fgets(line, sizeof line, f)) {
        int n, m; unsigned long long d; char par;
        if (sscanf(line, "%d %d %llu %c", &n, &m, &d, &par) != 4) continue;
        if (par == 'E' && d <= B && n < sz[d]) sz[d] = (uint8_t)n;
    }
    fclose(f);
    for (uint64_t v = 1; v <= B; v += 2) {
        if (sz[v] == 255) continue;
        int s = sz[v];
        /* cones and E_r joins over G u E_j */
        for (int t = s + 1; t < 64 && (1ULL << t) <= B; t++) {
            uint64_t c = (1ULL << t) + v;
            if (t >= 2 && c <= B && sz[c] > t + 1) sz[c] = (uint8_t)(t + 1);   /* t = 1 would be K_2 (odd) */
            for (int r = 2; r < 64; r++) {
                uint64_t a = ((1ULL << r) - 1);
                if (a > B / ((1ULL << t) - 1)) break;
                uint64_t val = a * ((1ULL << t) - 1) + 1 + v;
                if (val > B) break;
                if (v == 1 && s == 0 && r == t) continue;        /* E_r join E_r = K_{r,r} */
                if (sz[val] > t + r) sz[val] = (uint8_t)(t + r);
            }
        }
        /* (G) general join (A u E_i) join (G u E_j), i, j >= 1: d = (2^a - 1)(2^t - 1) + d(A) + d(G), size a + t;
         *     A ranges over final values alpha <= min(v, AMAX); excluded when (alpha, a) == (v, t) (possible swap). */
        if (GENJ) for (uint64_t al = 1; al <= v && al <= AMAX; al += 2) {
            if (sz[al] == 255) continue;
            if (al + v > B) break;
            int sa = sz[al];
            for (int a = sa + 1; a < 40; a++) {
                uint64_t A2 = (1ULL << a) - 1;
                if (A2 > B) break;
                for (int t = s + 1; t < 40; t++) {
                    uint64_t T2 = (1ULL << t) - 1;
                    if (A2 > (B - al - v) / T2) break;
                    uint64_t val = A2 * T2 + al + v;
                    if (val > B) break;
                    if (al == v && a == t) continue;
                    if (a + t < 255 && sz[val] > a + t) sz[val] = (uint8_t)(a + t);
                }
            }
        }
        /* products with already final smaller values w <= v (w >= 5) */
        if (v >= 5) for (uint64_t w = 5; w <= v && w <= B / v; w += 2) {
            if (sz[w] == 255) continue;
            uint64_t p = v * w; int ss = s + sz[w];
            if (ss < 255 && sz[p] > ss) sz[p] = (uint8_t)ss;
        }
    }
    long miss = 0; double worst = 0; uint64_t wv = 0;
    printf("missing odd values <= 2^%d:", L);
    for (uint64_t v = 1; v <= B; v += 2) {
        if (sz[v] == 255) { if (miss < 60) printf(" %llu", (unsigned long long)v); miss++; continue; }
        if (v > 15) { double r = sz[v] / log2((double)v); if (r > worst) { worst = r; wv = v; } }
    }
    printf("\ntotal missing %ld; worst size/log2(value) for values > 15: %.4f at %llu (size %d)\n", miss, worst,
           (unsigned long long)wv, wv ? sz[wv] : 0);
    /* size profile at the top */
    int maxs = 0; for (uint64_t v = B / 2 + 1; v <= B; v += 2) if (sz[v] != 255 && sz[v] > maxs) maxs = sz[v];
    printf("largest least-size among values in (2^%d, 2^%d]: %d\n", L - 1, L, maxs);
    return 0;
}
