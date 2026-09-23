// collatz_procgen_20260922_ladder_f2x.c -- the F_2[x] sibling of 3x+1 (Hicks-Mullen-Yucas-Zavislak 2008).
//
// Bit i of an unsigned integer is the coefficient of x^i.  "Odd" means f(0)=1.
// Shortcut map (as in the ladder note):  T(f) = f/x (f even),  ((x+1)f+1)/x = f ^ ((f^1)>>1) (f odd).
// HMYZ/ABP (non-shortcut) map: C(f) = (1+x)f+1 (f odd), f/x (f even); t_HMYZ(f) = min{k>=0: C^k(f)=1}
// = (shortcut steps) + (number of odd shortcut steps), since C(C(f)) = T(f) for odd f.
//
// usage: f2x [D=24] [K=20]
//  (1) identities: x*(T(f)+1) = (x+1)(f+1) for odd f, deg f <= 16;
//      T^k(x^k u + 1) = (x+1)^k u + 1 for odd u, k <= 12, deg u <= 12;
//  (2) for every nonzero f with deg f <= D: deg law, run-length law r(f) = v_x(f+1) (f != 1),
//      and the orbit reaches 1 (memoised census with in-progress marks = cycle detector);
//      per-degree max/mean shortcut and HMYZ stopping times vs the bound d^2+2d;
//      identity t_HMYZ(f) = 2 t_short(f) - deg f (Inselmann 2024) checked for every f;
//  (3) for k <= K: f mod x^k -> first k parities is a bijection; the all-odd word has the
//      single preimage 1 mod x^k (the F_2 exceptional count is 1 at every level).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static inline int deg32(uint32_t f) { return 31 - __builtin_clz(f); }  // f != 0
static inline uint32_t T(uint32_t f) { return (f & 1u) ? (f ^ ((f ^ 1u) >> 1)) : (f >> 1); }
static uint64_t clmul(uint64_t a, uint64_t b) { uint64_t r = 0; while (b) { if (b & 1) r ^= a; a <<= 1; b >>= 1; } return r; }
static void polystr(uint32_t f, char *buf) {  // e.g. x^5+x^2+1
  buf[0] = 0; int first = 1;
  for (int i = 31; i >= 0; i--) if (f >> i & 1u) {
    char t[16];
    if (i == 0) sprintf(t, "1"); else if (i == 1) sprintf(t, "x"); else sprintf(t, "x^%d", i);
    if (!first) strcat(buf, "+");
    strcat(buf, t); first = 0;
  }
}

int main(int argc, char **argv) {
  int D = argc > 1 ? atoi(argv[1]) : 24;
  int K = argc > 2 ? atoi(argv[2]) : 20;
  if (D < 1 || D > 28 || K < 1 || K > 26) { fprintf(stderr, "bad args\n"); return 1; }
  long fails = 0;

  // (1) identities
  long nid = 0;
  for (uint32_t f = 1; f < (1u << 17); f += 2) {  // odd, deg <= 16
    uint64_t lhs = (uint64_t)(T(f) ^ 1u) << 1, rhs = clmul(f ^ 1u, 3);
    if (lhs != rhs) fails++;
    nid++;
  }
  printf("(1a) x*(T(f)+1) == (x+1)*(f+1) for all %ld odd f with deg<=16: %s\n", nid, fails ? "FAIL" : "PASS");
  long nid2 = 0, f2 = 0;
  for (int k = 1; k <= 12; k++)
    for (uint32_t u = 1; u < (1u << 13); u += 2) {
      uint32_t f = (u << k) ^ 1u, g = f;
      for (int j = 0; j < k; j++) { if (!(g & 1u)) { f2++; break; } g = T(g); }
      uint64_t want = 1;
      for (int j = 0; j < k; j++) want = clmul(want, 3);
      want = clmul(want, u) ^ 1u;
      if ((uint64_t)g != want) f2++;
      nid2++;
    }
  fails += f2;
  printf("(1b) T^k(x^k u+1) == (x+1)^k u+1 with all k iterates odd, k<=12, odd u deg<=12 (%ld cases): %s\n",
         nid2, f2 ? "FAIL" : "PASS");

  // (2) census
  uint32_t SZ = 1u << (D + 1);
  uint16_t *ts = malloc(sizeof(uint16_t) * SZ), *th = malloc(sizeof(uint16_t) * SZ);
  uint32_t *stk = malloc(sizeof(uint32_t) * 65536);
  if (!ts || !th || !stk) { fprintf(stderr, "oom\n"); return 1; }
  const uint16_t UNK = 0xFFFF, INP = 0xFFFE;
  for (uint32_t f = 0; f < SZ; f++) ts[f] = th[f] = UNK;
  ts[1] = th[1] = 0;
  long degfail = 0, runfail = 0, cyclefail = 0;
  for (uint32_t f = 2; f < SZ; f++) {
    // local laws
    uint32_t g1 = T(f);
    int d = deg32(f);
    if (g1 == 0 || deg32(g1) != ((f & 1u) ? d : d - 1)) degfail++;
    int r = 0; uint32_t g = f;
    while (g & 1u) { r++; g = T(g); if (r > 64) break; }
    if (r != __builtin_ctz(f ^ 1u)) runfail++;
    if (ts[f] != UNK) continue;
    int sp = 0; g = f;
    while (ts[g] == UNK) {
      if (sp >= 65536) { fprintf(stderr, "stack\n"); return 1; }
      stk[sp++] = g; ts[g] = INP; g = T(g);
      if (g >= SZ) { fprintf(stderr, "degree increased?\n"); return 1; }
    }
    if (ts[g] == INP) { cyclefail++; for (int i = 0; i < sp; i++) ts[stk[i]] = UNK; continue; }
    uint16_t vs = ts[g], vh = th[g];
    while (sp > 0) {
      uint32_t c = stk[--sp];
      vs = vs + 1; vh = vh + ((c & 1u) ? 2 : 1);
      ts[c] = vs; th[c] = vh;
    }
  }
  long insfail = 0;
  for (uint32_t f = 1; f < SZ; f++)
    if (ts[f] >= INP || (int)th[f] != 2 * (int)ts[f] - deg32(f)) insfail++;
  fails += degfail + runfail + cyclefail + insfail;
  printf("(2a) deg T(f) = deg f (f odd) / deg f - 1 (f even), all %u nonzero f with deg<=%d: %s\n",
         SZ - 1, D, degfail ? "FAIL" : "PASS");
  printf("(2b) initial odd-run length == v_x(f+1) for all f != 1 with deg<=%d: %s\n", D, runfail ? "FAIL" : "PASS");
  printf("(2c) every nonzero f with deg<=%d reaches 1 (no other cycle): %s\n", D, cyclefail ? "FAIL" : "PASS");
  printf("(2d) t_HMYZ(f) == 2 t_short(f) - deg f for all nonzero f with deg<=%d: %s\n", D, insfail ? "FAIL" : "PASS");
  printf("  d  #polys  maxT_short argmax_short                  maxT_HMYZ  d^2+2d  maxH/d^1.5 maxH/(d ln d) meanT_short/d meanT_HMYZ/d\n");
  for (int d = 1; d <= D; d++) {
    uint32_t lo = 1u << d, hi = 1u << (d + 1);
    unsigned ms = 0, mh = 0; uint32_t as = 0; double sums = 0, sumh = 0;
    for (uint32_t f = lo; f < hi; f++) {
      if (ts[f] > ms) { ms = ts[f]; as = f; }
      if (th[f] > mh) mh = th[f];
      sums += ts[f]; sumh += th[f];
    }
    char buf[512]; polystr(as, buf);
    double n = (double)(hi - lo);
    printf("%3d %7u %10u  %-30s %9u %7d %10.3f %13.3f %13.4f %12.4f\n", d, hi - lo, ms, buf, mh, d * d + 2 * d,
           mh / pow(d, 1.5), d > 1 ? mh / (d * log((double)d)) : 0.0, sums / n / d, sumh / n / d);
  }
  if (D >= 20) {
    unsigned mh20 = 0, ms20 = 0;
    for (uint32_t f = 2; f < (1u << 21); f++) { if (th[f] > mh20) mh20 = th[f]; if (ts[f] > ms20) ms20 = ts[f]; }
    printf("  all 2^21-1 nonzero f with deg<=20: max shortcut steps %u, max HMYZ steps %u\n", ms20, mh20);
  }
  free(ts); free(th); free(stk);

  // (3) parity-vector bijection and the exceptional ladder
  unsigned char *seen = malloc((size_t)1 << K);
  long bijfail = 0;
  printf("(3) k : bijection(f mod x^k -> first k parities) ; #classes mod x^k with k odd steps (exceptional count)\n");
  for (int k = 1; k <= K; k++) {
    uint32_t N = 1u << k; memset(seen, 0, N);
    long allodd = 0; uint32_t who = 0; int ok = 1;
    for (uint32_t rr = 0; rr < N; rr++) {
      uint32_t g = rr, w = 0;
      for (int i = 0; i < k; i++) { w |= (g & 1u) << i; g = T(g); }
      if (seen[w]) ok = 0;
      seen[w] = 1;
      if (w == N - 1) { allodd++; who = rr; }
    }
    if (!ok || allodd != 1 || who != 1) bijfail++;
    if (k <= 4 || k % 4 == 0 || k == K)
      printf("    k=%2d : %s ; exceptional classes = %ld (class of %u)\n", k, ok ? "bijection" : "NOT bijective", allodd, who);
  }
  fails += bijfail;
  printf("(3) all k<=%d: %s\n", K, bijfail ? "FAIL" : "PASS");
  free(seen);
  printf("F2[x] TOTAL: %s\n", fails ? "SOME CHECK FAILED" : "ALL CHECKS PASS");
  return fails ? 1 : 0;
}
