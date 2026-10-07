/* B_L[t] = #(zeroless L-digit y with y == t mod 2^L); print min/max for each L (auditor's code). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
int main(int argc, char **argv) {
  int LMAX = atoi(argv[1]);
  uint64_t *B = calloc((size_t)1 << LMAX, 8), *Bn = calloc((size_t)1 << LMAX, 8);
  B[0] = 1;
  for (int L = 1; L <= LMAX; L++) {
    size_t sz = (size_t)1 << L; memset(Bn, 0, sz * 8);
    for (size_t y = 0; y < (sz >> 1); y++) { uint64_t v = B[y]; if (!v) continue; for (uint64_t a = 1; a <= 9; a++) Bn[(10 * y + a) & (sz - 1)] += v; }
    uint64_t *t = B; B = Bn; Bn = t;
    uint64_t mx = 0, mn = UINT64_MAX; size_t amx = 0, amn = 0;
    for (size_t s = 0; s < sz; s++) { if (B[s] > mx) { mx = B[s]; amx = s; } if (B[s] < mn) { mn = B[s]; amn = s; } }
    printf("L=%d B[0]=%llu max=%llu (s=%zu) min=%llu (s=%zu)\n", L, (unsigned long long)B[0], (unsigned long long)mx, amx, (unsigned long long)mn, amn);
    fflush(stdout);
  }
  return 0;
}
