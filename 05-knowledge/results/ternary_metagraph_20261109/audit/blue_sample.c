/* blue_sample.c - auditor's independent sampled check of the sign of the blue-line (1,3) excess at odd n >= 13
 * (HYP-9540 (D); the session's later PREREGISTRATION-4, U1 REFUTED at n = 13).  Written from scratch.
 *
 * For uniform random y (bits of the N two-element sigma-orbits) it evaluates f = (H - 1)/2 mod 2 = alpha_1 mod 2 (odd
 * cycle parity; OCF) at (y, 0), (y, e_i), (ybar, 0), (ybar, e_i) and (ybar, 1...1).  y is in D iff B(y) = B(ybar),
 * where B_i(y) = f(y, e_i) + f(y, 0) (affinity in x is THM-4922 (4) step 3; verified exhaustively for n <= 11 by
 * blue_parity.c); for y in D the line type e(y) = f(y, 0) + f(ybar, 1...1) is read directly (no use of (P2)).
 * share of (1,3) lines = 1/2 + (|D_1| - |D_0|) / 2^(N + 1); estimator (d1 - d0) / (2 S), standard error
 * sqrt((d1 + d0) / S - ((d1 - d0) / S)^2) / (2 sqrt(S)).
 * usage: blue_sample n SAMPLES SEED
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
typedef unsigned long long u64;
#define MAXN 17

static int n, m, c, N;
static int TI[160], TJ[160], IDX[MAXN][MAXN];
static int mir[10], orbA[80], orbB[80];

static int odd_cycle_parity(const uint32_t *in, uint32_t *dp) {
  int par = 0;
  for (int s = 0; s + 2 < n; s++) {
    int r = n - 1 - s, F = 1 << r;
    memset(dp, 0, sizeof(uint32_t) * F);
    dp[0] = 1u << s;
    for (int T = 0; T < F; T++) {
      uint32_t d = dp[T];
      if (!d) continue;
      int size = 1 + __builtin_popcount(T);
      if (size >= 3 && (size & 1)) par ^= __builtin_popcount(d & in[s]) & 1;
      for (int k = 0; k < r; k++) {
        if ((T >> k) & 1) continue;
        int u = s + 1 + k;
        if (__builtin_popcount(d & in[u]) & 1) dp[T | (1 << k)] ^= 1u << u;
      }
    }
  }
  return par;
}
static int f_eval(u64 y, int x, uint32_t *dp) {
  uint32_t out[MAXN] = {0}, in[MAXN] = {0};
  for (int v = 0; v + 1 < n; v++) out[v] |= 1u << (v + 1);
  for (int o = 0; o < N; o++) {
    int r = (int)((y >> o) & 1ULL);
    for (int h = 0; h < 2; h++) { int b = h ? orbB[o] : orbA[o]; if (r) out[TJ[b]] |= 1u << TI[b]; else out[TI[b]] |= 1u << TJ[b]; }
  }
  for (int i = 0; i < c; i++) { int b = mir[i]; if ((x >> i) & 1) out[TJ[b]] |= 1u << TI[b]; else out[TI[b]] |= 1u << TJ[b]; }
  for (int v = 0; v < n; v++) { uint32_t o = out[v]; while (o) { int u = __builtin_ctz(o); o &= o - 1; in[u] |= 1u << v; } }
  return odd_cycle_parity(in, dp);
}
static u64 splitmix(u64 *s) { u64 z = (*s += 0x9E3779B97F4A7C15ULL); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }

int main(int argc, char **argv) {
  n = atoi(argv[1]); long long S = atoll(argv[2]); u64 seed = strtoull(argv[3], 0, 10);
  c = (n - 1) / 2; m = 0;
  for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { TI[m] = i; TJ[m] = j; IDX[i][j] = m; m++; }
  int seen[160] = {0}; N = 0;
  for (int b = 0; b < m; b++) {
    int sb = IDX[n - 1 - TJ[b]][n - 1 - TI[b]];
    if (sb == b || seen[b]) continue;
    seen[b] = seen[sb] = 1; orbA[N] = b; orbB[N] = sb; N++;
  }
  for (int i = 0; i < c; i++) mir[i] = IDX[i][n - 1 - i];
  u64 yall = (N == 64) ? ~0ULL : ((1ULL << N) - 1);
  int xall = (1 << c) - 1;
  long long d0 = 0, d1 = 0, dtot = 0;
#pragma omp parallel reduction(+:d0,d1,dtot)
  {
    uint32_t *dp = malloc(sizeof(uint32_t) << (n - 1));
    u64 st = seed * 1000003ULL + (u64)omp_get_thread_num() * 0x1234567ULL + 17;
#pragma omp for schedule(dynamic, 64)
    for (long long s = 0; s < S; s++) {
      u64 y = splitmix(&st) & yall, yb = y ^ yall;
      int A = f_eval(y, 0, dp), Ab = f_eval(yb, 0, dp), B = 0, Bb = 0;
      for (int i = 0; i < c; i++) { B |= (f_eval(y, 1 << i, dp) ^ A) << i; Bb |= (f_eval(yb, 1 << i, dp) ^ Ab) << i; }
      if (B == Bb) {
        dtot++;
        int e = A ^ f_eval(yb, xall, dp);       /* the line {(y, 0), (ybar, 1...1)}: 1 iff (1,3) */
        if (e) d1++; else d0++;
      }
    }
    free(dp);
  }
  double est = (double)(d1 - d0) / (2.0 * S);
  double p1 = (double)d1 / S, p0 = (double)d0 / S;
  double se = sqrt((p1 + p0 - (p1 - p0) * (p1 - p0)) / S) / 2.0;
  printf("n=%d c=%d N=%d: %lld samples (seed %llu): |D|/2^N ~ %.5f, D_1 %lld, D_0 %lld; (1,3)-share = 1/2 %+.6f (se %.6f, z = %.1f)\n",
         n, c, N, S, seed, (double)dtot / S, d1, d0, est, se, est / se);
  return 0;
}
