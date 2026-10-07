/* Independent zeroless-powers-of-2 verifier (mac-mini-2026-10-07-oaimath3).
   Keeps 2^n mod 10^(9*NL) in base-10^9 limbs (little-endian), doubles each step, and finds the position of the
   rightmost decimal 0 (1 = units digit).  Prints every new record of that position, and fails loudly if a window of
   9*NL digits has no 0 (then the exponent must be checked in full).  Start state: 2^N0 mod 10^(9*NL) from stdin. */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#define NL 32
int main(int argc, char **argv) {
    unsigned long long n0 = strtoull(argv[1], 0, 10), n1 = strtoull(argv[2], 0, 10);
    uint32_t L[NL];
    for (int i = 0; i < NL; i++) if (scanf("%u", &L[i]) != 1) return 2;   /* limbs low to high */
    int rec = 0; unsigned long long n;
    for (n = n0; n < n1; n++) {
        /* rightmost zero of current value 2^n */
        int pos = -1;
        for (int i = 0; i < NL && pos < 0; i++) {
            uint32_t v = L[i];
            for (int d = 0; d < 9; d++) { if (v % 10 == 0) { pos = 9 * i + d + 1; break; } v /= 10; }
        }
        if (pos < 0) { printf("NO ZERO IN WINDOW at n=%llu\n", n); fflush(stdout); return 1; }
        if (pos > rec) { rec = pos; printf("n=%llu rightmost_zero_position=%d\n", n, pos); fflush(stdout); }
        /* double mod 10^(9*NL) */
        uint32_t c = 0;
        for (int i = 0; i < NL; i++) { uint64_t t = (uint64_t)L[i] * 2 + c; L[i] = (uint32_t)(t % 1000000000u); c = (uint32_t)(t / 1000000000u); }
    }
    printf("done n in [%llu, %llu): record rightmost-zero position %d; final limbs:", n0, n1, rec);
    for (int i = NL - 1; i >= 0; i--) printf(" %09u", L[i]);
    printf("\n");
    return 0;
}
