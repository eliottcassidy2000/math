/* verify2.c -- fast check that 2^n contains a decimal 0 for n in [n0,n1), trailing window only (288 digits).
   Tier A: limbs 0,1 (36 digits) of 2^n mod 10^288 updated every step (carry-chain-free doubling, base B=1e18).
   Tier B: limbs 2..15 kept lazily: H = floor(2^ns / B^2) mod B^14 at sync time ns, plus the carry accumulator
           C = sum_{t<g} 2^(g-1-t) b1(ns+t) of tier-A carry-outs, so that H(ns+g) = 2^g H(ns) + C (mod B^14).
   Multiplication by 2^r (r<=18) is chain-free: since B = 2^18 5^18,
           (2^r l) mod B = 2^r (l mod 2^(18-r) 5^18),  floor(2^r l / B) = floor(l / (2^(18-r) 5^18)).
   Sync happens when tier A has no zero (then tier B digits are needed) or when g reaches 56 (so C < 2^56 < B).
   Output: histogram of the position of the first zero from the right, records, survivors.
   usage: verify2 n0 n1 "<288-digit trailing window of 2^n0>" tag            */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#define NL 16
#define W (18*NL)
#define B 1000000000000000000ULL
static uint8_t fz6[1000000];
static uint64_t hist[W+2];
static inline int fz18(uint64_t x) { /* first zero from right within 18-digit limb, 1..18, or 0 */
  uint32_t c = (uint32_t)(x % 1000000ULL); if (fz6[c]) return fz6[c];
  x /= 1000000ULL; c = (uint32_t)(x % 1000000ULL); if (fz6[c]) return 6 + fz6[c];
  c = (uint32_t)(x / 1000000ULL); if (fz6[c]) return 12 + fz6[c];
  return 0;
}
static uint64_t P5_18;
static inline void mulpow2(uint64_t *h, int nl, int r) { /* h (nl limbs, little-endian) *= 2^r mod B^nl, r <= 18 */
  uint64_t d = P5_18 << (18 - r);   /* 2^(18-r) 5^18 */
  uint64_t prevq = 0;
  for (int j = 0; j < nl; j++) {
    uint64_t q = h[j] / d, rem = h[j] - q * d;
    h[j] = (rem << r) + prevq;
    prevq = q;
  }
}
int main(int argc, char **argv) {
  uint64_t n0 = strtoull(argv[1], 0, 10), n1 = strtoull(argv[2], 0, 10);
  const char *ts = argv[3], *tag = argv[4];
  P5_18 = 1; for (int i = 0; i < 18; i++) P5_18 *= 5;
  for (int v = 0; v < 1000000; v++) { int x = v; fz6[v] = 0; for (int i = 0; i < 6; i++) { if (x % 10 == 0) { fz6[v] = i + 1; break; } x /= 10; } }
  uint64_t l[NL];
  if ((int)strlen(ts) != W) { fprintf(stderr, "bad len\n"); return 1; }
  for (int j = 0; j < NL; j++) { uint64_t x = 0; const char *p = ts + W - 18*(j+1); for (int i = 0; i < 18; i++) x = 10*x + (p[i]-'0'); l[j] = x; }
  uint64_t a0 = l[0], a1 = l[1];
  uint64_t *h = l + 2; const int nh = NL - 2;
  uint64_t C = 0; int g = 0;
  int record = 0; uint64_t survivors = 0, syncs = 0;
  clock_t t0 = clock();
  for (uint64_t n = n0; n < n1; n++) {
    int pos = fz18(a0);
    if (!pos) { int p1 = fz18(a1); if (p1) pos = 18 + p1; }
    if (!pos || g >= 56) {  /* C < 2^56 < B keeps the ripple add single-step */
      /* sync tier B to time n */
      int gg = g;
      while (gg > 0) { int r = gg > 18 ? 18 : gg; mulpow2(h, nh, r); gg -= r; }
      /* add C with ripple carry */
      uint64_t cc = C;
      for (int j = 0; j < nh && cc; j++) { uint64_t s = h[j] + cc; if (s >= B) { h[j] = s - B; cc = 1; } else { h[j] = s; cc = 0; } }
      C = 0; g = 0; syncs++;
      if (!pos) {
        for (int j = 0; j < nh; j++) { int p = fz18(h[j]); if (p) { pos = 36 + 18*j + p; break; } }
      }
    }
    if (!pos) { survivors++; printf("%s SURVIVOR n=%llu\n", tag, (unsigned long long)n); fflush(stdout); pos = W + 1; }
    hist[pos]++;
    if (pos > record) { record = pos; printf("%s RECORD n=%llu first_zero_from_right=%d (zeroless suffix %d)\n", tag, (unsigned long long)n, pos, pos-1); fflush(stdout); }
    /* step tier A */
    uint64_t d0 = a0 << 1, b0 = d0 >= B; a0 = d0 - (b0 ? B : 0);
    uint64_t d1 = a1 << 1, b1 = d1 >= B; a1 = d1 - (b1 ? B : 0) + b0;
    C = (C << 1) | b1; g++;
    if (((n - n0) & ((1ULL<<33)-1)) == 0 && n > n0) { printf("%s progress n=%llu elapsed=%.1fs\n", tag, (unsigned long long)n, (double)(clock()-t0)/CLOCKS_PER_SEC); fflush(stdout); }
  }
  /* final sync for state check */
  { int gg = g; while (gg > 0) { int r = gg > 18 ? 18 : gg; mulpow2(h, nh, r); gg -= r; }
    uint64_t cc = C; for (int j = 0; j < nh && cc; j++) { uint64_t s = h[j] + cc; if (s >= B) { h[j] = s - B; cc = 1; } else { h[j] = s; cc = 0; } } }
  printf("%s DONE range=[%llu,%llu) count=%llu survivors=%llu record=%d syncs=%llu elapsed=%.1fs\n", tag, (unsigned long long)n0, (unsigned long long)n1,
         (unsigned long long)(n1-n0), (unsigned long long)survivors, record, (unsigned long long)syncs, (double)(clock()-t0)/CLOCKS_PER_SEC);
  printf("%s FINAL_STATE ", tag); printf("%018llu", (unsigned long long)l[NL-1]);
  for (int j = NL-2; j >= 2; j--) printf("%018llu", (unsigned long long)l[j]);
  printf("%018llu%018llu\n", (unsigned long long)a1, (unsigned long long)a0);
  printf("%s HIST", tag); for (int i = 1; i <= W+1; i++) printf(" %llu", (unsigned long long)hist[i]); printf("\n");
  return 0;
}
