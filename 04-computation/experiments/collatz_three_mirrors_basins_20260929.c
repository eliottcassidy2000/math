/* Basins of the saddle chain of 1 + 4^6 (4097, 3073, 2305, 1729, 1297, 973, 730), of 27 and of 137 under 3x+1
   (T(n) = n/2 even, (3n+1)/2 odd), S23 2026-09-29.
   For every n <= N: code[n] = the first target hit (a small index into the target list) or 0 if the orbit falls
   below n without hitting a target (then code[n] = code[v] of that smaller v; memoised).  Because the targets lie on
   two orbits (1729's orbit runs 4097 -> ... -> 730 -> 365 -> ... -> 137 -> ... and 27's orbit joins at 137), we
   record the FIRST target hit; the basin of a target a is then {n : the orbit of n passes through a} = the set of n
   whose first target is a or an ancestor of a on the chain.  Output per dyadic range: the counts, the cumulative
   densities, and the count against N^0.84 (Krasikov-Lagarias) for the root 1729.
   Build: gcc -O3 -o basins3 collatz_three_mirrors_basins_20260929.c ; run: ./basins3 1073741824 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128;
/* targets in orbit order along the chain: index 1..7 = 4097, 3073, 2305, 1729, 1297, 973, 730; 8 = 27; 9 = 137 */
static const uint64_t targets[10] = {0, 4097, 3073, 2305, 1729, 1297, 973, 730, 27, 137};
int main(int argc, char **argv) {
  uint64_t N = argc > 1 ? strtoull(argv[1], 0, 10) : (1ULL << 28);
  uint8_t *code = malloc(N + 1);
  if (!code) { fprintf(stderr, "alloc fail\n"); return 1; }
  memset(code, 0, N + 1);
  code[1] = 0;  /* the orbit of 1 is the cycle 1 -> 2 -> 1: never below 1, never a target */
  for (uint64_t n = 2; n <= N; n++) {
    u128 v = n;
    for (;;) {
      if (v <= 4097) {
        int hit = 0;
        for (int i = 1; i <= 9; i++) if (v == targets[i]) { code[n] = (uint8_t)i; hit = 1; break; }
        if (hit) break;
      }
      if (v < n) { code[n] = code[(uint64_t)v]; break; }
      if (v & 1) v = (3 * v + 1) / 2; else v >>= 1;
    }
  }
  int S = 0; while ((1ULL << (S + 1)) <= N) S++;
  /* basin membership: orbit passes through target a  <=>  first target hit is a or an earlier chain point (index <= a's index) for a in the chain 1..7;
     for 137 (index 9): first target in {1..7, 8, 9} all pass through 137 (the whole chain and 27 reach 137); for 27: first target 8 only. */
  printf("N = %llu; cumulative counts and densities of the basins B(a) = {n <= N : orbit passes through a}\n", (unsigned long long)N);
  printf("range          B(4097)   B(3073)   B(2305)   B(1729)   B(1297)   B(973)    B(730)    B(27)     B(137)\n");
  uint64_t cum[10] = {0};
  for (int s = 0; s <= S; s++) {
    uint64_t lo = 1ULL << s, hi = (1ULL << (s + 1)) - 1; if (hi > N) hi = N;
    uint64_t cnt[10] = {0};
    for (uint64_t n = lo; n <= hi; n++) cnt[code[n]]++;
    for (int i = 0; i < 10; i++) cum[i] += cnt[i];
    if (s >= 20 || s == S) {
      /* chain basins are nested: B(chain_i) = first target in {1..i} */
      uint64_t b[10]; double d[10];
      uint64_t tot = 0;
      for (int i = 1; i <= 7; i++) { tot += cum[i]; b[i] = tot; }
      b[8] = cum[8];
      b[9] = tot + cum[8] + cum[9];
      for (int i = 1; i <= 9; i++) d[i] = (double)b[i] / (double)hi;
      printf("[1,2^%2d]  ", s + 1);
      for (int i = 1; i <= 9; i++) printf(" %.6f", d[i]);
      printf("   |B(1729)| = %llu vs N^0.84 = %.3e (ratio %.3f)\n", (unsigned long long)b[4], pow((double)hi, 0.84), (double)b[4] / pow((double)hi, 0.84));
    }
  }
  /* first-target shares on the top range */
  {
    uint64_t lo = 1ULL << S, hi = N; uint64_t cnt[10] = {0};
    for (uint64_t n = lo; n <= hi; n++) cnt[code[n]]++;
    printf("top range [2^%d, N]: first-target shares:", S);
    for (int i = 1; i <= 9; i++) printf(" %llu:%.6f", (unsigned long long)targets[i], (double)cnt[i] / (double)(hi - lo + 1));
    printf("  none:%.6f\n", (double)cnt[0] / (double)(hi - lo + 1));
  }
  free(code);
  return 0;
}
