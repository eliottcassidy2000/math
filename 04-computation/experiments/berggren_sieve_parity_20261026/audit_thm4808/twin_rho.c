// Independent check of HYP-9418 (3): pi2(X) = #{p <= X : p, p+2 prime} by an odd-only segmented sieve, and
// rho(X) = (X/2) prod_{3 <= q <= sqrt X} (1 - 2/q) / pi2(X), at X = 10^3 .. 10^LMAX (own code).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv) {
    int LMAX = argc > 1 ? atoi(argv[1]) : 10;
    int64_t XMAX = 1; for (int i = 0; i < LMAX; i++) XMAX *= 10;
    int64_t LIM = XMAX + 2000;          // a prime beyond XMAX + 2 closes the last checkpoint
    int64_t r = (int64_t)sqrt((double)LIM) + 2;
    char *small = calloc(r + 1, 1);
    int64_t *pr = malloc(sizeof(int64_t) * (r + 1)); int np = 0;
    for (int64_t i = 2; i <= r; i++) if (!small[i]) { pr[np++] = i; for (int64_t j = i * i; j <= r; j += i) small[j] = 1; }
    const int64_t SEG = 1 << 22;                 // odd numbers per segment
    char *seg = malloc(SEG);
    int64_t next_cp = 1000, count = 0;
    int64_t prev_odd_prime = -10;                // last odd prime seen (for twin detection)
    int cpk = 3;
    // odd numbers 2i+1, i = 0.. ; handle 3,5,7 naturally (1 is not prime)
    for (int64_t lo = 0; 2 * lo + 1 <= LIM; lo += SEG) {
        int64_t hi = lo + SEG;                    // indices [lo, hi)
        memset(seg, 1, SEG);
        for (int k = 1; k < np; k++) {           // odd primes
            int64_t p = pr[k];
            int64_t start = p * p; if (start > 2 * (hi - 1) + 1) break;
            // first odd multiple >= max(p*p, 2lo+1)
            int64_t m0 = 2 * lo + 1; if (start < m0) { start = ((m0 + p - 1) / p) * p; if (start % 2 == 0) start += p; }
            for (int64_t j = (start - 1) / 2; j < hi; j += p) seg[j - lo] = 0;
        }
        if (lo == 0) seg[0] = 0;                  // 1 is not prime
        for (int64_t i = lo; i < hi; i++) {
            int64_t v = 2 * i + 1; if (v > LIM) break;
            if (!seg[i - lo]) continue;
            // v is an odd prime; the pair (v-2, v) is twin if prev == v - 2; count by lower member p = v - 2
            while (v - 2 > next_cp && next_cp <= XMAX) {
                // all lower members <= next_cp are counted once v - 2 > next_cp
                double prod = 1.0; int64_t s = (int64_t)floor(sqrt((double)next_cp));
                while ((s + 1) * (s + 1) <= next_cp) s++; while (s * s > next_cp) s--;
                for (int k = 1; k < np && pr[k] <= s; k++) prod *= 1.0 - 2.0 / (double)pr[k];
                double rho = (next_cp / 2.0) * prod / (double)count;
                printf("10^%d pi2=%lld prod=%.10f rho=%.5f\n", cpk, (long long)count, prod, rho);
                fflush(stdout);
                next_cp *= 10; cpk++;
            }
            if (prev_odd_prime == v - 2 && v - 2 <= XMAX) count++;
            prev_odd_prime = v;
        }
    }
    return 0;
}
