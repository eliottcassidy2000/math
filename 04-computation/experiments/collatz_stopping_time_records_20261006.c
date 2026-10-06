/* collatz_stopping_time_records_20261006.c
 *
 * "Wanted entries only" (Alman--Vassilevska Williams 2026): a last-dip violation in the residual window
 * (l, A, [X, N_l]) needs a source m* < n <= N_l whose T-stopping time exceeds A (it stays above n > m* for A
 * T-steps).  So instead of sweeping all last-dip segments it suffices to know  max { sigma_T(m) : m < N }.
 * This program computes the stopping-time records (first m attaining a new maximum of sigma_T) for odd m < X,
 * with sigma_T(m) = least j with T^j(m) < m, T(x) = (3x+1)/2 (odd) / x/2 (even), and prints max sigma_T below
 * each power of two.  Compare with the residual A's: 4701 (X = 2^20 level), 24727 (2^28), 125743 (2^32).
 *
 * Build: cc -O2 -o stoprec collatz_stopping_time_records_20261006.c ;  Run: ./stoprec [X]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
int main(int argc, char **argv) {
    uint64_t X = (argc > 1) ? strtoull(argv[1], 0, 10) : (1ULL << 32);
    uint64_t best = 0; uint64_t next_pow = 2;
    for (uint64_t m = 3; m < X; m += 2) {
        if (m >= next_pow) {
            printf("max sigma_T below 2^%d (%llu): %llu\n", (int)__builtin_ctzll(next_pow), (unsigned long long)next_pow, (unsigned long long)best);
            fflush(stdout);
            next_pow <<= 1;
        }
        u128 x = m; uint64_t j = 0;
        while (x >= m) {
            if (x & 1) x = (3 * x + 1) >> 1; else x >>= 1;
            j++;
            if (j > 100000000ULL) { printf("m=%llu exceeded 1e8 steps\n", (unsigned long long)m); return 2; }
        }
        if (j > best) { best = j; printf("record: sigma_T(%llu) = %llu\n", (unsigned long long)m, (unsigned long long)j); fflush(stdout); }
    }
    printf("max sigma_T below X=%llu: %llu\n", (unsigned long long)X, (unsigned long long)best);
    return 0;
}
