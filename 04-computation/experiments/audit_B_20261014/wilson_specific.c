/* Audit B: Wilson quotient w_p and the Giuga sum S(p) = sum_{k<p} k^(p-1) modulo p^2 for specific primes
 * (the further "Giuga primes" 780887 and 56151923 reported on primepuzzles.net puzzle 998 from the
 * Costa-Gerbicz-Harvey near-Wilson data).  (p-1)! by a sequential product mod p^2; S(p) directly for p < 10^6.
 * Usage: wilson_specific p [p ...] */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128; typedef uint64_t u64;
static u64 mm(u64 a, u64 b, u64 m) { return (u64)(((u128)a * b) % m); }
static u64 pw(u64 b, u64 e, u64 m) { u64 r = 1; b %= m; while (e) { if (e & 1) r = mm(r, b, m); b = mm(b, b, m); e >>= 1; } return r; }
int main(int argc, char **argv) {
  for (int i = 1; i < argc; i++) {
    u64 p = strtoull(argv[i], 0, 10), p2 = p * p, f = 1;
    for (u64 k = 2; k < p; k++) f = mm(f, k, p2);
    if ((f + 1) % p) { printf("%llu: not prime (Wilson fails)\n", (unsigned long long)p); continue; }
    u64 w = ((f + 1) / p) % p;
    printf("p = %llu: w_p mod p = %llu (p - 1 = %llu): w_p == -1: %d", (unsigned long long)p, (unsigned long long)w,
           (unsigned long long)(p - 1), w == p - 1);
    if (p < 1000000) {
      u64 S = 0; for (u64 k = 1; k < p; k++) { S += pw(k, p - 1, p2); if (S >= p2) S -= p2; }
      printf("; S(p) mod p^2 == -1: %d", S == p2 - 1);
    }
    printf("\n"); fflush(stdout);
  }
  return 0;
}
