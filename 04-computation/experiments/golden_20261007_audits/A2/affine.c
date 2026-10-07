#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main() {
    int KM = 26; long N = 1L << KM; uint8_t *lab = calloc(N + 2, 1);
    for (long m = 1; m <= N; m++) {
        if (m == 1) { lab[1] = 1; continue; } if (m == 5) { lab[5] = 2; continue; } if (m == 17) { lab[17] = 3; continue; }
        unsigned __int128 x = m;
        while (1) { x = (x & 1) ? (3 * x - 1) / 2 : x / 2; if (x < (unsigned __int128)m) { lab[m] = lab[(long)x]; break; } }
    }
    printf(" k   P[lab(m)!=lab(m+1)]  P[lab(m)!=lab(3m)]  P[lab(m)!=lab(m-7)]  (random-pair baseline 1-sum d^2 = 0.666)\n");
    for (int k = 10; k <= 24; k += 2) {
        long lo = 1L << (k - 1), hi = 1L << k, t = 0, c1 = 0, c3 = 0, c7 = 0;
        for (long m = lo; m < hi; m++) { t++; c1 += lab[m] != lab[m + 1]; c3 += lab[m] != lab[3 * m]; c7 += lab[m] != lab[m - 7]; }
        printf("%2d   %.4f   %.4f   %.4f\n", k, (double)c1 / t, (double)c3 / t, (double)c7 / t);
    }
    return 0;
}
