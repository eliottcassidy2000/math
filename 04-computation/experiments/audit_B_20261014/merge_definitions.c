/* Audit B (2026-10-14): two merge definitions for the orbits of n and n + h under T(n) = n/2, (3n+1)/2.
 *   COINC : the first equal-time coincidence T^k(n) = T^k(n+h) happens before either orbit reaches 1
 *           (THM-4627 statement 3 / collatz_merge_lags_20261013.py Part B: no odd-count condition);
 *   EQODD : the same, and the numbers of odd steps up to k agree (THM-4627 (2)(a)'s merge event; the condition
 *           used by prime_shift_merging.py).  After the first coincidence the orbits are identical, so EQODD holds
 *           iff the odd counts agree at the first coincidence.
 * Reported for odd n <= N, primes p <= N, lesser twin primes p <= N (p, p+2 prime), and also for ALL n <= N.
 * 128-bit integers are used for the orbit values.
 * Usage: merge_definitions N
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;

static int first_coinc(uint64_t n, uint64_t h, int *eqodd) {
    u128 x = n, y = n + h; long ax = 0, ay = 0;
    for (long k = 0; k < 100000; k++) {
        if (x == y) { *eqodd = (ax == ay); return 1; }
        if (x == 1 || y == 1) { *eqodd = 0; return 0; }
        if (x & 1) { x = (3 * x + 1) >> 1; ax++; } else x >>= 1;
        if (y & 1) { y = (3 * y + 1) >> 1; ay++; } else y >>= 1;
    }
    *eqodd = 0; return -1;
}

int main(int argc, char **argv) {
    uint64_t N = argc > 1 ? strtoull(argv[1], 0, 10) : 10000000ULL;
    unsigned char *comp = calloc(N + 4, 1);
    comp[0] = comp[1] = 1;
    for (uint64_t i = 2; i * i <= N + 3; i++) if (!comp[i]) for (uint64_t j = i * i; j <= N + 3; j += i) comp[j] = 1;
    for (int h = 1; h <= 2; h++) {
        /* classes: 0 all n, 1 odd n, 2 odd primes, 3 lesser twin primes */
        double cnt[4] = {0}, co[4] = {0}, eo[4] = {0}; long unf = 0;
        for (uint64_t n = 1; n <= N; n++) {
            int e, c = first_coinc(n, h, &e);
            if (c < 0) unf++;
            int cls[4] = {1, (int)(n & 1), (n > 2 && !comp[n]), (n > 2 && !comp[n] && !comp[n + 2])};
            for (int t = 0; t < 4; t++) if (cls[t]) { cnt[t]++; if (c == 1) co[t]++; if (c == 1 && e) eo[t]++; }
        }
        const char *nm[4] = {"all n", "odd n", "odd primes", "lesser twin primes"};
        for (int t = 0; t < 4; t++)
            printf("N = %llu h = %d %-20s count %10.0f  COINC %.5f  EQODD %.5f  (coincidences with unequal odd counts: %.0f)\n",
                   (unsigned long long)N, h, nm[t], cnt[t], co[t] / cnt[t], eo[t] / cnt[t], co[t] - eo[t]);
        printf("unfinished: %ld\n", unf);
        fflush(stdout);
    }
    return 0;
}
