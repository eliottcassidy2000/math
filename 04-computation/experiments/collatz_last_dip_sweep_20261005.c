/* collatz_last_dip_sweep_20261005.c
 *
 * Exhaustive test of the LAST-DIP LEMMA below X (default 2^28):
 *   for odd m* < n with Syracuse orbit m* = n_0 -> ... -> n_l = n and every n_j > n (1 <= j < l),
 *   the valuation word is rising: 2^A < 3^l.
 * For every odd m < X the Syracuse orbit is followed until it drops below m (the early break is sound:
 * if the orbit dips to n' < m and later visits n'' > m, the last dip into n'' is found from n' < m).
 * For every odd orbit value n with m < n < X, the last orbit value v < n before n is located and the
 * segment word (l, A) from v to n is tested.  Also records the maximum of 2^A/3^l (as a log2 excess
 * A - l*log2(3) <= 0 means rising) and the count of segments.
 *
 * Build: cc -O2 -o collatz_last_dip_sweep collatz_last_dip_sweep_20261005.c
 * Run:   ./collatz_last_dip_sweep [X]           (X <= 2^40; memory O(orbit length))
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

typedef unsigned __int128 u128;

int main(int argc, char **argv) {
    uint64_t X = (argc > 1) ? strtoull(argv[1], 0, 10) : (1ULL << 28);
    const double LOG2_3 = 1.5849625007211563;
    uint64_t *orb = malloc(sizeof(uint64_t) * 100000);
    uint32_t *orbl = malloc(sizeof(uint32_t) * 100000);
    uint32_t *orbA = malloc(sizeof(uint32_t) * 100000);
    unsigned long long segments = 0, violations = 0;
    double best_excess = -1e9; uint64_t best_v = 0, best_n = 0; uint32_t best_l = 0, best_A = 0;
    uint64_t maxval = 0;
    int longest = 0;
    for (uint64_t m = 3; m < X; m += 2) {
        int len = 0;
        orb[0] = m; orbl[0] = 0; orbA[0] = 0; len = 1;
        uint64_t n = m; uint32_t l = 0, A = 0;
        for (;;) {
            u128 y = (u128)3 * n + 1;
            uint32_t a = 0;
            while ((y & 1) == 0) { y >>= 1; a++; }
            if (y > (u128)0xFFFFFFFFFFFFFFFFULL) { fprintf(stderr, "overflow at m=%llu\n", (unsigned long long)m); return 2; }
            n = (uint64_t)y; l++; A += a;
            if (n > maxval) maxval = n;
            if (n < m) break;
            if (n == 1) break;
            if (n < X) {
                /* last orbit value below n */
                for (int k = len - 1; k >= 0; k--) {
                    if (orb[k] < n) {
                        uint32_t dl = l - orbl[k], dA = A - orbA[k];
                        double excess = (double)dA - (double)dl * LOG2_3;   /* > 0 means 2^A > 3^l (NON-rising) */
                        segments++;
                        if (excess > best_excess) { best_excess = excess; best_v = orb[k]; best_n = n; best_l = dl; best_A = dA; }
                        if (excess > 0) {
                            /* exact check with integers to avoid floating error: 2^dA > 3^dl ? */
                            u128 p2 = (u128)1 << (dA > 126 ? 126 : dA);
                            u128 p3 = 1; for (uint32_t t = 0; t < dl && p3 < ((u128)1 << 126); t++) p3 *= 3;
                            if (dA > 126 || p2 > p3) {
                                violations++;
                                printf("VIOLATION: m*=%llu n=%llu l=%u A=%u\n", (unsigned long long)orb[k], (unsigned long long)n, dl, dA);
                            }
                        }
                        break;
                    }
                }
            }
            if (len >= 100000) { fprintf(stderr, "orbit too long at m=%llu\n", (unsigned long long)m); return 3; }
            orb[len] = n; orbl[len] = l; orbA[len] = A; len++;
            if (len > longest) longest = len;
        }
    }
    printf("X=%llu last-dip segments=%llu violations=%llu\n", (unsigned long long)X, segments, violations);
    printf("largest A - l*log2(3) = %.9f at m*=%llu n=%llu l=%u A=%u  (2^A/3^l = %.9f)\n",
           best_excess, (unsigned long long)best_v, (unsigned long long)best_n, best_l, best_A, pow(2.0, best_excess));
    printf("longest orbit prefix stored=%d, max orbit value=%llu\n", longest, (unsigned long long)maxval);
    return violations ? 1 : 0;
}
