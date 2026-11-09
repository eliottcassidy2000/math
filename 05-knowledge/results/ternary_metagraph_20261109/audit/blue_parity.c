/* blue_parity.c - auditor's independent computation of THM-4951 (5c), (5d) at odd n (written from scratch).
 *
 * For every blue (sigma-fixed) tiling t at odd n = 2c + 1, with x = the c mirror tiles (i, n-1-i) and y = the bits of
 * the N two-element sigma-orbits, compute f(t) = (H(t) - 1)/2 mod 2 = alpha_1(t) mod 2 (OCF: H = 1 + 2 alpha_1 + 4(..)),
 * where alpha_1 mod 2 is the parity of the number of directed odd cycles (parity path DP, no Hamiltonian paths and no
 * mirror-path decomposition).  Then, from the full table f(y, x):
 *   - the blue lines {t, t + 1} by (H mod 4, H mod 4);
 *   - the affinity f(y, x) = A(y) + x.B(y) of THM-4922 (4) step 3, tested for EVERY (y, x);
 *   - (P2) sum_i B_i(y) = 1 for every y;  D = {y : B(y) = B(y + 1)} and D_e (THM-4951 (5d)), and the formula
 *     #(1,3) = 2^(N + c - 2) + 2^(c - 2) (|D_1| - |D_0|).
 * usage: blue_parity n
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
typedef unsigned long long u64;
#define MAXN 13

static int n, m, c, N;
static int TI[80], TJ[80], IDX[MAXN][MAXN];
static int mir[8];            /* tile index of mirror tile i */
static int orbA[40], orbB[40];/* the two tiles of orbit o */

static int odd_cycle_parity(const uint16_t *in) {
  /* parity of the number of directed cycles of odd length >= 3 */
  static uint16_t dp[1 << MAXN];
#pragma omp threadprivate(dp)
  int par = 0;
  for (int s = 0; s + 2 < n; s++) {
    int r = n - 1 - s;                       /* vertices s+1..n-1 */
    int F = 1 << r;
    for (int T = 0; T < F; T++) dp[T] = 0;
    /* dp[T]: bit v (absolute vertex) = parity of paths s -> ... -> v covering {s} + (T << (s+1)) */
    dp[0] = (uint16_t)(1u << s);
    for (int T = 0; T < F; T++) {
      uint16_t d = dp[T];
      if (!d) continue;
      int size = 1 + __builtin_popcount(T);
      if (size >= 3 && (size & 1)) par ^= __builtin_popcount(d & in[s]) & 1;
      for (int k = 0; k < r; k++) {
        if ((T >> k) & 1) continue;
        int u = s + 1 + k;
        if (__builtin_popcount(d & in[u]) & 1) dp[T | (1 << k)] ^= (uint16_t)(1u << u);
      }
    }
  }
  return par;
}

int main(int argc, char **argv) {
  n = atoi(argv[1]);
  c = (n - 1) / 2;
  m = 0;
  for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { TI[m] = i; TJ[m] = j; IDX[i][j] = m; m++; }
  int seen[80] = {0}; N = 0; int nm = 0;
  for (int b = 0; b < m; b++) {
    int sb = IDX[n - 1 - TJ[b]][n - 1 - TI[b]];
    if (sb == b) { mir[nm++] = b; continue; }
    if (seen[b]) continue;
    seen[b] = seen[sb] = 1; orbA[N] = b; orbB[N] = sb; N++;
  }
  if (nm != c) { fprintf(stderr, "mirror count %d != c %d\n", nm, c); return 1; }
  /* order mirror tiles by i */
  for (int i = 0; i < c; i++) mir[i] = IDX[i][n - 1 - i];
  u64 NY = 1ULL << N; int NX = 1 << c;
  uint8_t *f = calloc(NY * NX / 8 + 1, 1);
#pragma omp parallel for schedule(dynamic, 256)
  for (long long y = 0; y < (long long)NY; y++) {
    for (int x = 0; x < NX; x++) {
      uint16_t out[MAXN] = {0}, in[MAXN] = {0};
      for (int v = 0; v + 1 < n; v++) out[v] |= 1u << (v + 1);
      for (int o = 0; o < N; o++) {
        int r = (y >> o) & 1;
        for (int h = 0; h < 2; h++) { int b = h ? orbB[o] : orbA[o]; if (r) out[TJ[b]] |= 1u << TI[b]; else out[TI[b]] |= 1u << TJ[b]; }
      }
      for (int i = 0; i < c; i++) { int b = mir[i]; if ((x >> i) & 1) out[TJ[b]] |= 1u << TI[b]; else out[TI[b]] |= 1u << TJ[b]; }
      for (int v = 0; v < n; v++) for (int u = 0; u < n; u++) if ((out[v] >> u) & 1u) in[u] |= 1u << v;
      int p = odd_cycle_parity(in);
      if (p) {
        u64 idx = (u64)y * NX + x;
#pragma omp atomic
        f[idx >> 3] |= (uint8_t)(1u << (idx & 7));
      }
    }
  }
#define F_(y, x) ((f[(((u64)(y)) * NX + (x)) >> 3] >> ((((u64)(y)) * NX + (x)) & 7)) & 1)
  u64 yall = NY - 1; int xall = NX - 1;
  long long c11 = 0, c13 = 0, c33 = 0, affine_bad = 0, p2_bad = 0, D = 0, D1 = 0, D0 = 0, ones = 0;
  for (u64 y = 0; y < NY; y++) {
    int A = F_(y, 0), B = 0;
    for (int i = 0; i < c; i++) B |= (F_(y, 1 << i) ^ A) << i;
    for (int x = 0; x < NX; x++) {
      int pred = A ^ (__builtin_popcount(x & B) & 1);
      if (pred != F_(y, x)) affine_bad++;
      ones += F_(y, x);
      int g = F_(y ^ yall, x ^ xall);
      int a = F_(y, x);
      if (a == 0 && g == 0) c11++; else if (a == 1 && g == 1) c33++; else c13++;
    }
    if (!(__builtin_popcount(B) & 1)) p2_bad++;
    u64 yb = y ^ yall;
    int Ab = F_(yb, 0), Bb = 0;
    for (int i = 0; i < c; i++) Bb |= (F_(yb, 1 << i) ^ Ab) << i;
    if (B == Bb) { D++; if ((A ^ Ab ^ 1) == 1) D1++; else D0++; }
  }
  printf("n=%d c=%d N=%d: blue tilings %llu; H = 3 (mod 4) on %lld of them\n", n, c, N, NY * NX, ones);
  printf("  blue lines by (H mod 4) pair: (1,1) %lld (1,3) %lld (3,3) %lld (each line counted once)\n", c11 / 2, c13 / 2, c33 / 2);
  printf("  affinity f = A(y) + x.B(y) failures: %lld of %llu; (P2) failures: %lld of %llu\n", affine_bad, NY * NX, p2_bad, NY);
  printf("  |D| = %lld, |D_1| = %lld, |D_0| = %lld\n", D, D1, D0);
  long long pred13 = (1LL << (N + c - 2)) + ((c >= 2) ? (1LL << (c - 2)) * (D1 - D0) : 0);
  printf("  formula 2^(N+c-2) + 2^(c-2)(|D_1| - |D_0|) = %lld (direct %lld)\n", pred13, c13 / 2);
  return 0;
}
