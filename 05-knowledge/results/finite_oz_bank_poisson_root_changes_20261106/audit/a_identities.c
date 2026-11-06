/* a_identities.c - independent audit of THM-4888 (A) and (B) (auditor's own code).
 * Q is built as a full array by the defining recurrence, in increasing m:
 *     Q(0) = 0, Q(2k) = Q(k), Q(2k+1) = Q(k) + k + (k mod 2).
 * Checks, for every m < 2^NB:
 *   (i)   Session 206 digit identity  Q(m) = sum_j b_j floor(m / 2^(j+1)) + #{j : b_j = b_(j+1) = 1}   (exact)
 *   (ii)  THM-4888 (A)  Q(m) = m^2 - r(m) (mod 3), r(m) = number of maximal runs of 1-bits         (exact)
 *   (iii) Session 205 bounds  Q even,  0 <= Q(m) <= m - 1 (m >= 1)
 *   (iv)  Session 206 involution T(m) = 2^n - 1 - rev_(n-1)(m - 2^(n-1)) (n = bit length):
 *         T(T(m)) = m, T preserves bit length, and Q(T(m)) - Q(m) = 2 (T(m) - m)
 * Prints T(9), T(4) and counts of failures.  usage: a_identities NB   (NB <= 30)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned long long u64;
static int bitlen(u64 m){ int n = 0; while (m){ n++; m >>= 1; } return n; }
static u64 rev(u64 a, int k){ u64 r = 0; for (int i = 0; i < k; i++){ r = (r << 1) | (a & 1); a >>= 1; } return r; }
static u64 T(u64 m){ if (!m) return 0; int n = bitlen(m); u64 a = m - (1ULL << (n - 1)); return (1ULL << n) - 1 - rev(a, n - 1); }
static int runs(u64 m){ int r = 0, prev = 0; while (m){ int b = m & 1; if (b && !prev) r++; prev = b; m >>= 1; } return r; }
int main(int argc, char **argv){
  int NB = argc > 1 ? atoi(argv[1]) : 24; u64 N = 1ULL << NB;
  uint32_t *Q = malloc(N * sizeof(uint32_t)); if (!Q){ fprintf(stderr, "oom\n"); return 1; }
  Q[0] = 0;
  for (u64 m = 1; m < N; m++){ u64 k = m >> 1; Q[m] = (m & 1) ? (uint32_t)(Q[k] + k + (k & 1)) : Q[k]; }
  printf("T(9) = %llu, T(4) = %llu, T(11) = %llu, T(7) = %llu\n", T(9), T(4), T(11), T(7));
  u64 bad_digit = 0, bad_mod3 = 0, bad_even = 0, bad_bound = 0, bad_T = 0, bad_inv = 0, fixed = 0;
  for (u64 m = 0; m < N; m++){
    u64 q = Q[m];
    /* (i) digit identity */
    u64 s = 0; for (int j = 0; (m >> j); j++) if ((m >> j) & 1) s += m >> (j + 1);
    s += (u64)__builtin_popcountll(m & (m >> 1));
    if (s != q) bad_digit++;
    /* (ii) mod-3 law */
    long long lhs = (long long)(q % 3), rhs = (long long)(((m % 3) * (m % 3)) % 3) - runs(m) % 3;
    if (((lhs - rhs) % 3 + 3) % 3) bad_mod3++;
    /* (iii) */
    if (q & 1) bad_even++;
    if (m >= 1 && q > m - 1) bad_bound++;
    /* (iv) */
    if (m >= 1){
      u64 t = T(m);
      if (T(t) != m || bitlen(t) != bitlen(m)) bad_inv++;
      if (t < N && (long long)Q[t] - (long long)q != 2 * ((long long)t - (long long)m)) bad_T++;
      if (t == m) fixed++;
    }
  }
  printf("NB=%d checked m < %llu\n", NB, N);
  printf("digit identity failures: %llu\nmod-3 law failures: %llu\nQ odd: %llu\nQ > m-1: %llu\n", bad_digit, bad_mod3, bad_even, bad_bound);
  printf("T not an involution of the dyadic block: %llu\nQ(T m) - Q(m) != 2(T m - m): %llu\nfixed points of T (1 <= m < 2^NB): %llu\n", bad_inv, bad_T, fixed);
  /* sanity: the law must FAIL for a wrong sign, so the check is not vacuous */
  u64 wrong = 0; for (u64 m = 0; m < 4096; m++){ long long l = Q[m] % 3, r = ((m % 3) * (m % 3) + runs(m)) % 3; if (((l - r) % 3 + 3) % 3) wrong++; }
  printf("control (wrong law Q = m^2 + r mod 3) failures below 4096: %llu (must be > 0)\n", wrong);
  free(Q); return 0;
}
