// collatz_procgen_20260922_ladder_eq_dfs.c -- exact per-class branch-and-bound for the choice games E_q^S.
//
// Independent code path for lane one's DP (exceptional_general.c).  Moves on Z_2 (all arithmetic is
// exact modulo 2^prec, prec = bits of precision still known):
//   odd x  -> q x + 1                       (forced, multiplier q, no bit consumed)
//   even x -> x/2                           (multiplier 1/2, one bit consumed)
//   even x -> q^2 x + q + 1  (excursion)    (multiplier q^2, no bit; mode 1: always, mode 2: iff x mod 2^J in S)
// A class mod 2^m is exceptional iff no legal path of precision <= m reaches a state with q^a < 2^b
// (a multiplications, b halvings).  Exact integer tests: bmin[a] = least b with 2^b > q^a.
//
// usage:  eq_dfs q mode exhaustive m [J r1 r2 ...]      -- all classes mod 2^m, exact count
//         eq_dfs q mode mc m nsamples seed [J r1 ...]    -- uniform random classes mod 2^m (m <= 64)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static uint64_t Q, Q2, ADD;
static int M, MODE, J = 1;
static unsigned char allowS[1 << 12];
static int bmin[256];
static uint64_t nodes;

static int dfs(uint64_t r, int prec, int a) {
  nodes++;
  if (a >= 255 || bmin[a] > M) return 0;  // even halving every remaining bit cannot descend
  if (prec == 0) return 0;                // parity unknown: no certified move
  uint64_t mask = (prec >= 64) ? ~0ULL : ((1ULL << prec) - 1);
  if (r & 1ULL) return dfs((Q * r + 1ULL) & mask, prec, a + 1);
  int b = M - prec;
  if (bmin[a] <= b + 1) return 1;         // halving now descends: q^a < 2^(b+1)
  if (dfs(r >> 1, prec - 1, a)) return 1;
  if (MODE == 1 || (MODE == 2 && prec >= J && allowS[r & ((1ULL << J) - 1)]))
    if (dfs((Q2 * r + ADD) & mask, prec, a + 2)) return 1;
  return 0;
}

static uint64_t sm_state;
static uint64_t splitmix64(void) {
  uint64_t z = (sm_state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "usage: see header\n"); return 1; }
  Q = strtoull(argv[1], 0, 10); MODE = atoi(argv[2]);
  int mc = strcmp(argv[3], "mc") == 0;
  M = atoi(argv[4]);
  if (M < 1 || M > 64) { fprintf(stderr, "m must be in 1..64\n"); return 1; }
  long nsamp = 0; uint64_t seed = 0; int argi = 5;
  if (mc) { nsamp = atol(argv[5]); seed = strtoull(argv[6], 0, 10); argi = 7; }
  if (MODE == 2) { J = atoi(argv[argi++]); for (; argi < argc; argi++) allowS[atoi(argv[argi]) & ((1 << J) - 1)] = 1; }
  Q2 = Q * Q; ADD = Q + 1;
  // bmin[a] = least b with 2^b > q^a, exact: long double is exact enough here (a*log2 q never within 1e-12 of an integer
  // for q <= 99, a < 255), and verified against integer arithmetic while q^a fits in 64 bits.
  for (int a = 0; a < 256; a++) {
    long double v = a * log2l((long double)Q);
    bmin[a] = (int)floorl(v) + 1;
  }
  { uint64_t p = 1; int a = 0;
    while (1) { int b = 0; while (b < 64 && (1ULL << b) <= p) b++; if (b != bmin[a]) { fprintf(stderr, "bmin check failed a=%d\n", a); return 1; }
      if (p > UINT64_MAX / Q / 2) break; p *= Q; a++; } }
  long bad = 0, n = 0;
  if (!mc) {
    if (M > 30) { fprintf(stderr, "exhaustive needs m<=30\n"); return 1; }
    for (uint64_t x = 0; x < (1ULL << M); x++) { n++; if (!dfs(x, M, 0)) bad++; }
    printf("q=%llu mode=%d exhaustive m=%d exceptional=%ld of %ld frac=%.6f nodes=%llu\n",
           (unsigned long long)Q, MODE, M, bad, n, (double)bad / n, (unsigned long long)nodes);
  } else {
    sm_state = seed;
    uint64_t maskM = (M >= 64) ? ~0ULL : ((1ULL << M) - 1);
    for (long i = 0; i < nsamp; i++) { uint64_t x = splitmix64() & maskM; n++; if (!dfs(x, M, 0)) bad++; }
    double f = (double)bad / n, se = sqrt(f * (1 - f) / n);
    printf("q=%llu mode=%d mc m=%d samples=%ld exceptional=%ld frac=%.5f +- %.5f (1 s.e.) nodes=%llu\n",
           (unsigned long long)Q, MODE, M, n, bad, f, se, (unsigned long long)nodes);
  }
  return 0;
}
