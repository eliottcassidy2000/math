/* Transitive-plus-one strong tournaments (w_1 = 1, w_m = 0): check 1/3 <= P(u_2 before z) <= 2/3 for every word, m <= MMAX,
   by direct enumeration of all valid A (sets of u's before z): A valid iff (A empty or w_maxA = 0) and (B empty or w_minB = 1).
   Independent of the closed-form decomposition used in the proof.  Usage: ./tpo_u2_check MMAX  (cost about 8^m/4; m <= 20 takes minutes) */
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv) {
    int MMAX = atoi(argv[1]);
    for (int m = 3; m <= MMAX; m++) {
        long long words = 0, bad = 0, eq_hi = 0, eq_lo = 0;
        double pmin = 1, pmax = 0;
        for (long long mid = 0; mid < (1LL << (m - 2)); mid++) {
            int w = 1 | (int)(mid << 1);            /* bit i-1 = w_i ; w_1 = 1, w_m = 0 */
            long long H = 0, in2 = 0;
            for (long long A = 0; A < (1LL << m); A++) {
                long long B = ((1LL << m) - 1) ^ A;
                if (A) { int mx = 63 - __builtin_clzll(A); if (w >> mx & 1) continue; }
                if (B) { int mn = __builtin_ctzll(B); if (!(w >> mn & 1)) continue; }
                H++; if (A >> 1 & 1) in2++;
            }
            words++;
            if (3 * in2 < H || 3 * in2 > 2 * H) { bad++; if (bad < 4) printf("  VIOLATION m=%d w=%x H=%lld in2=%lld\n", m, w, H, in2); }
            if (3 * in2 == 2 * H) eq_hi++;
            if (3 * in2 == H) eq_lo++;
            double p = (double)in2 / H; if (p < pmin) pmin = p; if (p > pmax) pmax = p;
        }
        printf("m=%2d (n=%2d): %lld words, violations %lld, P(u_2 before z) in [%.5f, %.5f], equality 2/3: %lld, 1/3: %lld\n",
               m, m + 1, words, bad, pmin, pmax, eq_hi, eq_lo);
        fflush(stdout);
    }
    return 0;
}
