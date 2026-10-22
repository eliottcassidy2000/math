/* a3_hcount.c -- AUDIT (independent of circ_counts.c and ie_ham.c): number H of Hamiltonian paths of a
 * circulant tournament Cay(Z_n, S), S = connection set (exactly one of k, n-k for k = 1..(n-1)/2).
 *
 *   mode dp : Held-Karp over subsets of Z_n \ {0}, paths starting at 0; H = n * #(paths from 0)  (n <= 23, uint64).
 *   mode ie : Karp inclusion-exclusion  H = sum_{X subset Z_n} (-1)^(n-|X|) W(X),  W(X) = #walks on n vertices inside X,
 *             summed over binary necklaces (rotation orbits; iterative FKM, weight = period).  Arithmetic in
 *             unsigned __int128, i.e. modulo 2^128; exact because 0 <= H < 2^128 for n <= 32.
 *
 * stdin : lines "n s_1 ... s_m";  stdout : "n s_1 ... s_m : H".   Usage: ./a3_hcount dp|ie < sets.txt
 * Audit of THM-4751 (2026-10-21).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned __int128 u128;

static void print_u128(u128 x) {
    char buf[64]; int i = 63; buf[i] = 0;
    if (x == 0) { printf("0"); return; }
    while (x) { buf[--i] = (char)('0' + (int)(x % 10)); x /= 10; }
    printf("%s", buf + i);
}

static int n;
static uint32_t outm[33];   /* outm[i] = bitmask of out-neighbours of i */

static void build(const int *S, int ms) {
    int in[64] = {0};
    for (int t = 0; t < ms; t++) in[((S[t] % n) + n) % n] = 1;
    for (int i = 0; i < n; i++) {
        outm[i] = 0;
        for (int d = 1; d < n; d++) if (in[d]) outm[i] |= 1u << ((i + d) % n);
    }
    /* tournament check */
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        int a = (outm[i] >> j) & 1, b = (outm[j] >> i) & 1;
        if (a + b != 1) { fprintf(stderr, "not a tournament: n=%d pair %d %d\n", n, i, j); exit(2); }
    }
}

/* ---------------- dp ---------------- */
static u128 h_dp(void) {
    int N1 = n - 1;                 /* vertices 1..n-1 -> bits 0..n-2 */
    size_t full = (size_t)1 << N1;
    uint64_t *f = calloc(full * (size_t)N1, sizeof(uint64_t));
    if (!f) { fprintf(stderr, "dp: out of memory\n"); exit(3); }
    uint32_t o[32];
    for (int v = 1; v < n; v++) o[v - 1] = (outm[v] >> 1);          /* drop vertex 0, shift */
    uint32_t o0 = outm[0] >> 1;
    for (int b = 0; b < N1; b++) if ((o0 >> b) & 1) f[((size_t)1 << b) * N1 + b] = 1;
    for (size_t mask = 1; mask < full; mask++) {
        uint64_t *row = f + mask * N1;
        uint32_t mm = (uint32_t)mask;
        while (mm) {
            int v = __builtin_ctz(mm); mm &= mm - 1;
            uint64_t c = row[v]; if (!c) continue;
            uint32_t nx = o[v] & ~(uint32_t)mask;
            while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; f[(mask | ((size_t)1 << w)) * N1 + w] += c; }
        }
    }
    u128 tot = 0; uint64_t *last = f + (full - 1) * N1;
    for (int v = 0; v < N1; v++) tot += last[v];
    free(f);
    return tot * (u128)n;
}

/* ---------------- ie ---------------- */
static u128 walks(uint32_t X) {
    u128 v[32], w[32]; int idx[32], k = 0;
    for (int i = 0; i < n; i++) if ((X >> i) & 1) idx[k++] = i;
    if (k == 0) return 0;
    for (int a = 0; a < k; a++) v[idx[a]] = 1;
    for (int step = 1; step < n; step++) {
        for (int a = 0; a < k; a++) {
            int i = idx[a]; uint32_t mm = outm[i] & X; u128 s = 0;
            while (mm) { int j = __builtin_ctz(mm); mm &= mm - 1; s += v[j]; }
            w[i] = s;
        }
        for (int a = 0; a < k; a++) v[idx[a]] = w[idx[a]];
    }
    u128 t = 0; for (int a = 0; a < k; a++) t += v[idx[a]];
    return t;
}

static u128 h_ie(long *nneck) {
    int a[34]; memset(a, 0, sizeof a);   /* a[1..n] */
    u128 acc = 0; long cnt = 0;
    /* visit all-zero necklace: X = empty, W = 0 */
    cnt++;
    for (;;) {
        int i = n;
        while (i > 0 && a[i] == 1) i--;
        if (i == 0) break;
        a[i] = 1;
        for (int j = i + 1; j <= n; j++) a[j] = a[j - i];
        if (n % i == 0) {
            uint32_t X = 0; int sz = 0;
            for (int t = 1; t <= n; t++) if (a[t]) { X |= 1u << (t - 1); sz++; }
            u128 wv = walks(X) * (u128)i;     /* orbit size = period i */
            if ((n - sz) & 1) acc -= wv; else acc += wv;
            cnt++;
        }
    }
    *nneck = cnt;
    return acc;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s dp|ie < sets\n", argv[0]); return 1; }
    int mode_dp = !strcmp(argv[1], "dp");
    char line[4096];
    while (fgets(line, sizeof line, stdin)) {
        int vals[64], k = 0; char *p = line, *e;
        for (;;) { long x = strtol(p, &e, 10); if (e == p) break; vals[k++] = (int)x; p = e; }
        if (k < 1) continue;
        n = vals[0];
        if (n > 32 || (mode_dp && n > 23)) { fprintf(stderr, "n too large for mode\n"); return 1; }
        build(vals + 1, k - 1);
        u128 H; long nn = 0;
        if (mode_dp) H = h_dp(); else H = h_ie(&nn);
        for (int t = 0; t < k; t++) printf(t ? " %d" : "%d", vals[t]);
        printf(" : "); print_u128(H);
        if (!mode_dp) printf("  [necklaces=%ld]", nn);
        printf("\n"); fflush(stdout);
    }
    return 0;
}
