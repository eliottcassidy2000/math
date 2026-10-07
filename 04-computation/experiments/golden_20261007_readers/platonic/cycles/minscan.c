/* minscan.c -- all integer cycles of the Terras map T(x) = x/2 (x even), (3x+1)/2 (x odd)
 * with period <= LMAX, found by scanning the candidate minimal-|.| element of each cycle.
 * A cycle of period L <= LMAX with k odd steps has min element x_min <= XB (positive side)
 * or |x|_min <= YB (negative side); XB, YB come from bounds.py (product identity
 * prod_{odd}(3 + 1/x_i) = 2^L).  For each candidate c we iterate until the orbit
 *   - returns to c          -> c is the min-|.| element of a cycle (period = return time),
 *   - drops below c in |.|   -> c is not the min element of a cycle,
 *   - runs LMAX steps        -> c is not on a cycle of period <= LMAX (counted).
 * Negative side: x = -y, y -> y/2 (y even), (3y-1)/2 (y odd).
 * usage: minscan XB YB LMAX
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;

static void report(int sign, uint64_t c, int L) {
    /* recompute word and k */
    char w[256]; int k = 0; u128 v = c; int t;
    for (t = 0; t < L && t < 255; t++) {
        int b = (int)(v & 1); w[t] = b ? '1' : '0'; k += b;
        if (sign > 0) v = b ? (3 * v + 1) >> 1 : v >> 1;
        else          v = b ? (3 * v - 1) >> 1 : v >> 1;
    }
    w[t] = 0;
    printf("CYCLE sign=%c min|x|=%llu period L=%d odd steps k=%d word(from min)=%s\n",
           sign > 0 ? '+' : '-', (unsigned long long)c, L, k, L < 255 ? w : "(long)");
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: minscan XB YB LMAX\n"); return 1; }
    uint64_t XB = strtoull(argv[1], 0, 10), YB = strtoull(argv[2], 0, 10);
    int LMAX = atoi(argv[3]);
    uint64_t capped = 0, ovf = 0, maxsteps_pos = 0, maxsteps_neg = 0; u128 maxval = 0;
    for (uint64_t c = 1; c <= XB; c++) {
        u128 v = c; int t;
        for (t = 1; t <= LMAX; t++) {
            v = (v & 1) ? (3 * v + 1) >> 1 : v >> 1;
            if (v > maxval) maxval = v;
            if (v == c) { report(+1, c, t); break; }
            if (v < c) break;
            if (v >> 120) { ovf++; break; }
        }
        if (t > LMAX) capped++;
        if ((uint64_t)t > maxsteps_pos && t <= LMAX) maxsteps_pos = t;
    }
    for (uint64_t c = 1; c <= YB; c++) {
        u128 v = c; int t;
        for (t = 1; t <= LMAX; t++) {
            v = (v & 1) ? (3 * v - 1) >> 1 : v >> 1;
            if (v > maxval) maxval = v;
            if (v == c) { report(-1, c, t); break; }
            if (v < c) break;
            if (v >> 120) { ovf++; break; }
        }
        if (t > LMAX) capped++;
        if ((uint64_t)t > maxsteps_neg && t <= LMAX) maxsteps_neg = t;
    }
    printf("done XB=%llu YB=%llu LMAX=%d capped=%llu overflow=%llu max_stop_pos=%llu max_stop_neg=%llu maxval~%.3e\n",
           (unsigned long long)XB, (unsigned long long)YB, LMAX, (unsigned long long)capped,
           (unsigned long long)ovf, (unsigned long long)maxsteps_pos, (unsigned long long)maxsteps_neg,
           (double)maxval);
    return 0;
}
