/* Auditor's independent MITM for Z_K = #(zeroless K-digit x, 2^K | x), different split from the reader (m=15).
   Low part: m digits enumerated directly as integers s = sum a_i 10^i with pruning 2^(i+1) | s_{i+1}.
   Then x = s + 10^m y, y zeroless L digits; need 2^(m+L) | x <=> y == -5^(-m) (s/2^m) mod 2^L.
   H[r] = #low parts with -5^(-m) (s/2^m) == r mod 2^LMAX (uint64), folded to coarser L.
   B_L[t] = #(zeroless L-digit y == t mod 2^L), built by y = 10 y' + a.
   Z_{m+L} = sum_t H_L[t] B_L[t].   usage: a3_zk_mitm2 m LMAX                                   */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
static uint64_t *H; static uint64_t MASK, NEGINV; static int M;
static uint64_t P10[20], leaves = 0;
static void pr(u128 x){char b[64];int i=63;b[i]=0;if(!x){printf("0");return;}while(x){b[--i]='0'+(int)(x%10);x/=10;}printf("%s",b+i);}
static void dfs(int i, uint64_t s) {
  if (i == M) { uint64_t q = s >> M; H[(NEGINV * q) & MASK]++; leaves++; return; }
  for (uint64_t a = 1; a <= 9; a++) {
    uint64_t t = s + a * P10[i];
    if ((t & ((2ULL << i) - 1)) == 0) dfs(i + 1, t);   /* need 2^(i+1) | t */
  }
}
int main(int argc, char **argv) {
  M = atoi(argv[1]); int LMAX = atoi(argv[2]);
  P10[0] = 1; for (int i = 1; i < 20; i++) P10[i] = P10[i-1] * 10;
  uint64_t p5m = 1; for (int i = 0; i < M; i++) p5m *= 5;
  uint64_t inv = 1; for (int it = 0; it < 7; it++) inv *= 2 - p5m * inv;   /* inverse mod 2^64 */
  if (p5m * inv != 1) { fprintf(stderr, "inverse failed\n"); return 1; }
  NEGINV = 0 - inv; MASK = (1ULL << LMAX) - 1;
  H = calloc((size_t)1 << LMAX, sizeof(uint64_t));
  dfs(0, 0);
  printf("m=%d leaves=Z_m=%llu\n", M, (unsigned long long)leaves); fflush(stdout);
  uint64_t *B = calloc((size_t)1 << LMAX, sizeof(uint64_t)), *Bn = calloc((size_t)1 << LMAX, sizeof(uint64_t));
  uint64_t *Hc = malloc(((size_t)1 << LMAX) * sizeof(uint64_t));
  B[0] = 1;  /* L = 0: empty y */
  for (int L = 1; L <= LMAX; L++) {
    size_t sz = (size_t)1 << L, half = sz >> 1;
    memset(Bn, 0, sz * sizeof(uint64_t));
    for (size_t y = 0; y < half; y++) if (B[y]) for (uint64_t a = 1; a <= 9; a++) Bn[(10 * y + a) & (sz - 1)] += B[y];
    uint64_t *t = B; B = Bn; Bn = t;
    /* H folded to modulus 2^L: Hc[r] = sum over r' == r mod 2^L */
    memset(Hc, 0, sz * sizeof(uint64_t));
    for (size_t r = 0; r < ((size_t)1 << LMAX); r++) Hc[r & (sz - 1)] += H[r];
    u128 Z = 0; for (size_t r = 0; r < sz; r++) Z += (u128)Hc[r] * B[r];
    printf("Z_%d = ", M + L); pr(Z); printf("\n"); fflush(stdout);
  }
  return 0;
}
