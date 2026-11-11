/* aud_reach.c - auditor's independent recomputation of the pre-registered reach runs E1 and E3 (pure towers).
 * Written from Theorem B / Corollary B of THM-5130 (level map on states), not from the session's tower_reach.c.
 *
 * State (r, M, h, E).  A level is (delta in {1,2}, x in [0, XMAX], m in [0, MMAX]): x = (sinks added to the base)
 * + (tail vertices below the pair); with s sinks and i = x - s + 1 the level map of Corollary B gives
 *     r+ = -r,  h+ = r + 2^(x+1+delta+m) M,  M+ = alpha 2^m (h + r) + 2 gamma 2^x 3^m M - [x = 0] E,
 *     E+ = alpha 2^x M + [m = 0] (r - h)        (alpha = 2^delta + 1, gamma = 3^delta + 2^delta)
 * (a sink doubles M and kills E, so the dependence is on x only).
 * Start states: (TT_(k+1), {1}) = (1, 2^k, 1, [k = 0]), k = 0..KMAX, and (TT, {}) = (1, 0, 1, 0).
 * Output: for each depth d, the number of odd c <= CMAX occurring as odd part of M with offset +1 / -1 at depth <= d;
 * number of offset -1 states whose M is a power of two; least depth and 2-adic valuation for chosen c.
 * usage: aud_reach XMAX MMAX KMAX DEPTH CMAX
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef __int128 i128;
static int XMAX, MMAX, KMAX, DEPTH, CMAX;
static unsigned char *minDepthPlus, *minDepthMinus;   /* least depth at which c occurs (255 = never) */
static unsigned char *v2Plus, *v2Minus;               /* a valuation at the least depth (least one) */
static unsigned long long pow2minus = 0, states = 0;
static i128 P3[8];

static void record(int r, i128 M, int depth) {
    if (M <= 0) return;
    int v = 0;
    /* count trailing zeros of a 128-bit value */
    unsigned long long lo = (unsigned long long)M, hi = (unsigned long long)(M >> 64);
    if (lo) v = __builtin_ctzll(lo); else v = 64 + __builtin_ctzll(hi);
    i128 c = M >> v;
    if (r == -1 && c == 1) pow2minus++;
    if (c > CMAX) return;
    int ci = (int)c;
    unsigned char *md = (r == 1) ? minDepthPlus : minDepthMinus;
    unsigned char *vv = (r == 1) ? v2Plus : v2Minus;
    if (depth < md[ci]) { md[ci] = (unsigned char)depth; vv[ci] = (unsigned char)(v > 250 ? 250 : v); }
    else if (depth == md[ci] && v < vv[ci]) vv[ci] = (unsigned char)v;
}

static void go(int r, i128 M, i128 h, i128 E, int depth) {
    states++;
    if (r == 1 || r == -1) record(r, M, depth);
    if (depth == DEPTH) return;
    for (int delta = 1; delta <= 2; delta++) {
        i128 alpha = (delta == 1) ? 3 : 5, gamma = (delta == 1) ? 5 : 13;
        for (int x = 0; x <= XMAX; x++) for (int m = 0; m <= MMAX; m++) {
            i128 Mp = alpha * ((i128)1 << m) * (h + r) + 2 * gamma * ((i128)1 << x) * P3[m] * M - (x == 0 ? E : 0);
            i128 hp = r + ((i128)1 << (x + 1 + delta + m)) * M;
            i128 Ep = alpha * ((i128)1 << x) * M + (m == 0 ? (r - h) : 0);
            go(-r, Mp, hp, Ep, depth + 1);
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: aud_reach XMAX MMAX KMAX DEPTH CMAX\n"); return 1; }
    XMAX = atoi(argv[1]); MMAX = atoi(argv[2]); KMAX = atoi(argv[3]); DEPTH = atoi(argv[4]); CMAX = atoi(argv[5]);
    P3[0] = 1; for (int i = 1; i < 8; i++) P3[i] = 3 * P3[i - 1];
    minDepthPlus = malloc(CMAX + 2); minDepthMinus = malloc(CMAX + 2); v2Plus = malloc(CMAX + 2); v2Minus = malloc(CMAX + 2);
    memset(minDepthPlus, 255, CMAX + 2); memset(minDepthMinus, 255, CMAX + 2); memset(v2Plus, 255, CMAX + 2); memset(v2Minus, 255, CMAX + 2);
    /* optional 6th argument "pairs": start instead from the two transitive pairs (TT_2,{1,2}) = (-1, 3, 1, 1) and
       (TT_3,{1,3}) = (-1, 5, 1, 1) (the exploratory run e1x); optional 7th argument: file receiving the list
       "sign c leastdepth" */
    int pairs = (argc > 6 && strcmp(argv[6], "pairs") == 0);
    if (!pairs) {
        for (int k = 0; k <= KMAX; k++) go(1, (i128)1 << k, 1, k == 0 ? 1 : 0, 0);
        go(1, 0, 1, 0, 0);
    } else {
        go(-1, 3, 1, 1, 0);
        go(-1, 5, 1, 1, 0);
    }
    if (argc > 7) {
        FILE *f = fopen(argv[7], "w");
        for (int c = 1; c <= CMAX; c += 2) {
            if (minDepthPlus[c] != 255) fprintf(f, "+1 %d %d\n", c, minDepthPlus[c]);
            if (minDepthMinus[c] != 255) fprintf(f, "-1 %d %d\n", c, minDepthMinus[c]);
        }
        fclose(f);
    }
    printf("box: delta in {1,2}, 0 <= x <= %d, 0 <= m <= %d, start k <= %d (+ empty out-set), depth <= %d, odd c <= %d\n", XMAX, MMAX, KMAX, DEPTH, CMAX);
    for (int d = 0; d <= DEPTH; d++) {
        int np = 0, nm = 0;
        for (int c = 1; c <= CMAX; c += 2) { if (minDepthPlus[c] <= d) np++; if (minDepthMinus[c] <= d) nm++; }
        printf("depth <= %d: offset +1: %d odd c ; offset -1: %d\n", d, np, nm);
    }
    printf("states visited: %llu ; offset -1 states with M a power of two: %llu\n", states, pow2minus);
    printf("odd c in [3,129] never reached with offset +1:");
    for (int c = 3; c <= 129; c += 2) if (minDepthPlus[c] == 255) printf(" %d", c);
    printf("\n");
    int least = -1; for (int c = 3; c <= CMAX; c += 2) if (minDepthPlus[c] == 255) { least = c; break; }
    printf("least odd c >= 3 not reached with offset +1: %d\n", least);
    int watch[] = {1, 5, 11, 19, 23, 29, 73, 85, 91, 95, 125, 127};
    for (unsigned i = 0; i < sizeof watch / sizeof watch[0]; i++) {
        int c = watch[i]; if (c > CMAX) continue;
        printf("  c = %d: +1 at depth %d (v2 %d) ; -1 at depth %d (v2 %d)   [255 = never]\n", c, minDepthPlus[c], v2Plus[c], minDepthMinus[c], v2Minus[c]);
    }
    return 0;
}
