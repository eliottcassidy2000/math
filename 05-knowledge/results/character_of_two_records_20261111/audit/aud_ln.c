/* AUDIT program (independent): l_n = number of labelled tournaments on {1..n} in which the largest vertex
 * reaches every vertex, by brute force over all 2^C(n,2) tournaments (Gray code), n <= 8.
 * Also counts the strong ones. Compare with A001187 (connected labelled graphs): 1, 1, 4, 38, 728, 26704,
 * 1866256, 251548592.
 * usage: aud_ln n
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char **argv) {
    int n = argc > 1 ? atoi(argv[1]) : 5;
    if (n < 1 || n > 8) return 1;
    int pi[28], pj[28], np = 0, i, j;
    for (i = 0; i < n; i++) for (j = i + 1; j < n; j++) { pi[np] = i; pj[np] = j; np++; }
    unsigned out[8];
    /* start: bits = 0: j -> i for all i < j */
    for (i = 0; i < n; i++) out[i] = 0;
    for (i = 0; i < np; i++) out[pj[i]] |= 1u << pi[i];
    uint64_t total = 1ULL << np, g, ln = 0, sn = 0;
    unsigned full = (1u << n) - 1;
    for (g = 0; g < total; g++) {
        if (g) {
            int b = __builtin_ctzll(g);
            out[pi[b]] ^= 1u << pj[b];
            out[pj[b]] ^= 1u << pi[b];
        }
        /* reach set of the largest vertex n-1 */
        unsigned r = 1u << (n - 1), fr = r;
        while (fr) {
            unsigned nx = 0, t = fr;
            while (t) { int v = __builtin_ctz(t); t &= t - 1; nx |= out[v]; }
            fr = nx & ~r;
            r |= fr;
        }
        if (r == full) {
            ln++;
            /* strong iff also every vertex reaches n-1: reach in the converse */
            unsigned in[8];
            for (i = 0; i < n; i++) in[i] = ~out[i] & full & ~(1u << i);
            unsigned r2 = 1u << (n - 1), f2 = r2;
            while (f2) {
                unsigned nx = 0, t = f2;
                while (t) { int v = __builtin_ctz(t); t &= t - 1; nx |= in[v]; }
                f2 = nx & ~r2;
                r2 |= f2;
            }
            if (r2 == full) sn++;
        }
    }
    printf("n = %d: tournaments %llu; largest vertex reaches every vertex: %llu; strong: %llu\n",
           n, (unsigned long long)total, (unsigned long long)ln, (unsigned long long)sn);
    return 0;
}
