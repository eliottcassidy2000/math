/* General basin sieve for up to 8 targets under 3x+1 (T(n) = n/2 even, (3n+1)/2 odd), S23 2026-09-29.
   code[n] = bitmask of the targets that the orbit of n passes through: the targets met before the orbit first falls
   below n, OR the mask of that smaller value (memoised).  No orbit relation between the targets is assumed.
   Output: per dyadic range the cumulative density of each basin B(a) = {n <= N : orbit passes through a}.
   Build: gcc -O3 -o basins_mask collatz_three_mirrors_basins_mask_20260929.c ; run: ./basins_mask N t1 t2 ... t8 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: basins_mask N t1 [t2 ... t8]\n"); return 1; }
  uint64_t N = strtoull(argv[1], 0, 10);
  int T = argc - 2; if (T > 8) T = 8;
  uint64_t targets[8]; uint64_t tmax = 0;
  for (int i = 0; i < T; i++) { targets[i] = strtoull(argv[2 + i], 0, 10); if (targets[i] > tmax) tmax = targets[i]; }
  uint8_t *code = malloc(N + 1);
  if (!code) { fprintf(stderr, "alloc fail\n"); return 1; }
  memset(code, 0, N + 1);
  code[1] = 0;
  for (uint64_t n = 2; n <= N; n++) {
    u128 v = n; uint8_t mask = 0;
    for (;;) {
      if (v <= tmax) for (int i = 0; i < T; i++) if (v == targets[i]) mask |= (uint8_t)(1u << i);
      if (v < n) { mask |= code[(uint64_t)v]; break; }
      if (v & 1) v = (3 * v + 1) / 2; else v >>= 1;
    }
    code[n] = mask;
  }
  int S = 0; while ((1ULL << (S + 1)) <= N) S++;
  printf("N = %llu; targets:", (unsigned long long)N);
  for (int i = 0; i < T; i++) printf(" %llu", (unsigned long long)targets[i]);
  printf("\nrange     "); for (int i = 0; i < T; i++) printf(" B(%llu)", (unsigned long long)targets[i]); printf("\n");
  uint64_t cum[8] = {0};
  for (int s = 0; s <= S; s++) {
    uint64_t lo = 1ULL << s, hi = (1ULL << (s + 1)) - 1; if (hi > N) hi = N;
    for (uint64_t n = lo; n <= hi; n++) { uint8_t m = code[n]; for (int i = 0; i < T; i++) if (m & (1u << i)) cum[i]++; }
    if (s >= 20 || s == S) {
      printf("[1,2^%2d] ", s + 1);
      for (int i = 0; i < T; i++) printf(" %.6f", (double)cum[i] / (double)hi);
      printf("\n");
    }
  }
  free(code);
  return 0;
}
