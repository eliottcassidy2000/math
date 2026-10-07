/* verify.c -- check that 2^n contains a decimal 0 for n in [n0, n1).
   Trailing window: 2^n mod 10^(18*L) kept EXACTLY in L limbs of base 1e18, doubled each step.
   Carry-chain-free doubling: since B=1e18 is even, carry out of limb j is (2*l_j >= B),
   independent of the carry into limb j (2*l_j + 1 >= B with 2*l_j < B would force 2*l_j = B-1, odd).
   A zero inside the trailing window is a genuine digit of 2^n whenever 2^n has more than 18*L digits
   (true for n >= 957 when L = 16), hence certifies that 2^n contains a 0.
   Leading window (statistics only): first 54 significant digits of 2^n, truncated, 3 limbs;
   relative error <= (#divisions by 10)*1e-53, so >= 36 leading digits are reliable for n <= 1e12.
   Output: suffix histogram (position of first zero from the right), joint (prefix,suffix) histogram,
   record holders, survivors.
   usage: verify n0 n1 "<288-digit decimal of 2^n0 mod 10^288>" "<54 leading digits of 2^n0>" tag
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#define L 16
#define W (18*L)
#define B 1000000000000000000ULL
#define PMAX 36
static uint8_t fz6[1000000];   /* position (1..6) of first zero from the right in 6-digit zero-padded v, 0 if none */
static uint8_t lz6[1000000];   /* position (1..6) of first zero from the LEFT in 6-digit zero-padded v, 0 if none */
static uint64_t hist[W+2];
static uint64_t joint[PMAX+1][W+2];
static inline int suffix_pos(const uint64_t *l) {
  for (int j = 0; j < L; j++) {
    uint64_t x = l[j];
    for (int c = 0; c < 3; c++) {
      uint32_t ch = (uint32_t)(x % 1000000ULL); x /= 1000000ULL;
      if (fz6[ch]) return 18*j + 6*c + fz6[ch];
    }
  }
  return 0;
}
/* leading digits: h[0] (18 digits, in [1e17,1e18)), h[1], h[2] (18 digits each) */
static inline int prefix_zeroless_len(const uint64_t *h) {
  /* returns number of leading zeroless digits, capped at PMAX */
  int cnt = 0;
  for (int j = 0; j < 3; j++) {
    uint64_t x = h[j];
    uint32_t ch[3];
    ch[2] = (uint32_t)(x % 1000000ULL); x /= 1000000ULL;
    ch[1] = (uint32_t)(x % 1000000ULL); x /= 1000000ULL;
    ch[0] = (uint32_t)x;
    for (int c = 0; c < 3; c++) {
      if (lz6[ch[c]]) { cnt += lz6[ch[c]] - 1; return cnt > PMAX ? PMAX : cnt; }
      cnt += 6;
      if (cnt >= PMAX) return PMAX;
    }
  }
  return cnt > PMAX ? PMAX : cnt;
}
int main(int argc, char **argv) {
  if (argc < 6) { fprintf(stderr, "usage\n"); return 1; }
  uint64_t n0 = strtoull(argv[1], 0, 10), n1 = strtoull(argv[2], 0, 10);
  const char *ts = argv[3], *ls = argv[4], *tag = argv[5];
  for (int v = 0; v < 1000000; v++) {
    int d[6], x = v;
    for (int i = 0; i < 6; i++) { d[i] = x % 10; x /= 10; } /* d[0] least significant */
    fz6[v] = 0; for (int i = 0; i < 6; i++) if (d[i] == 0) { fz6[v] = i + 1; break; }
    lz6[v] = 0; for (int i = 5; i >= 0; i--) if (d[i] == 0) { lz6[v] = 6 - i; break; }
  }
  uint64_t l[L], h[3];
  if ((int)strlen(ts) != W) { fprintf(stderr, "bad trailing string len %zu\n", strlen(ts)); return 1; }
  for (int j = 0; j < L; j++) { /* limb j = digits [18j, 18j+18) from the right */
    uint64_t x = 0; const char *p = ts + W - 18*(j+1);
    for (int i = 0; i < 18; i++) x = 10*x + (p[i]-'0');
    l[j] = x;
  }
  if ((int)strlen(ls) != 54) { fprintf(stderr, "bad leading string\n"); return 1; }
  for (int j = 0; j < 3; j++) { uint64_t x = 0; for (int i = 0; i < 18; i++) x = 10*x + (ls[18*j+i]-'0'); h[j] = x; }
  int record = 0; uint64_t survivors = 0;
  int bestcov_n = -1; double bestcov = 0;
  clock_t t0 = clock();
  for (uint64_t n = n0; n < n1; n++) {
    int pos = suffix_pos(l);
    int pre = prefix_zeroless_len(h);
    if (pos == 0) { survivors++; printf("%s SURVIVOR n=%llu (no zero in last %d digits)\n", tag, (unsigned long long)n, W); fflush(stdout); pos = W+1; }
    hist[pos]++;
    joint[pre][pos]++;
    if (pos > record) { record = pos; printf("%s RECORD n=%llu first_zero_from_right=%d (zeroless suffix %d) prefix=%d\n", tag, (unsigned long long)n, pos, pos-1, pre); fflush(stdout); }
    /* trailing doubling, chain-free */
    uint64_t b[L];
    for (int j = 0; j < L; j++) { uint64_t d = l[j] << 1; b[j] = d >= B; l[j] = d - (b[j] ? B : 0); }
    for (int j = L-1; j >= 1; j--) l[j] += b[j-1];
    /* leading doubling with truncation */
    {
      uint64_t d2 = h[2] << 1, c2 = d2 >= B; h[2] = d2 - (c2 ? B : 0);
      uint64_t d1 = (h[1] << 1) + c2, c1 = d1 >= B; h[1] = d1 - (c1 ? B : 0);
      h[0] = (h[0] << 1) + c1;
      if (h[0] >= B) { /* divide 54-digit number by 10, truncating */
        uint64_t r0 = h[0] % 10; h[0] /= 10;
        uint64_t r1 = h[1] % 10; h[1] = h[1] / 10 + r0 * 100000000000000000ULL;
        h[2] = h[2] / 10 + r1 * 100000000000000000ULL;
      }
    }
    if (((n - n0) & ((1ULL<<32)-1)) == 0 && n > n0) {
      double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
      printf("%s progress n=%llu elapsed=%.1fs\n", tag, (unsigned long long)n, el); fflush(stdout);
    }
  }
  double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
  printf("%s DONE range=[%llu,%llu) count=%llu survivors=%llu record=%d elapsed=%.1fs\n", tag,
         (unsigned long long)n0, (unsigned long long)n1, (unsigned long long)(n1-n0), (unsigned long long)survivors, record, el);
  /* final state check values */
  printf("%s FINAL_TRAIL_LOW18=%018llu FINAL_LEAD_HI18=%llu\n", tag, (unsigned long long)l[0], (unsigned long long)h[0]);
  printf("%s HIST", tag); for (int i = 1; i <= W+1; i++) printf(" %llu", (unsigned long long)hist[i]); printf("\n");
  printf("%s JOINT_BEGIN\n", tag);
  for (int p = 0; p <= PMAX; p++) { printf("%s J %d", tag, p); for (int i = 1; i <= W+1; i++) printf(" %llu", (unsigned long long)joint[p][i]); printf("\n"); }
  printf("%s JOINT_END\n", tag);
  return 0;
}
