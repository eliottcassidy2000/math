/* zk_mitm.c -- exact Z_K = #{zeroless K-digit strings x (leading zeros forbidden) with 2^K | x}
   = A181610(K) = #{n mod 4*5^(K-1): last K digits of 2^n are all nonzero}.
   Meet in the middle: low m digits enumerated by the q-process (q_{i+1} = (q_i + a 5^i)/2, a == q_i mod 2),
   giving the multiset Q_m = {x_low / 2^m}; high L digits y (zeroless, L digits) must satisfy
   q + 5^m y == 0 mod 2^L, i.e. y == -5^{-m} q mod 2^L.  B_L[s] = #{zeroless L-digit y == s mod 2^L}.
   Z_{m+L} = sum_s H_L[s] B_L[s], H_L = histogram of r(q) = -5^{-m} q mod 2^L.
   Also prints min/max of B_L (rigorous bounds Z_{k+L}/Z_k in [min B_L, max B_L]) and the DFS level counts.
   usage: zk_mitm m Lmax
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
static void print_u128(u128 x) { char buf[64]; int i = 63; buf[i] = 0; if (x == 0) { printf("0"); return; } while (x) { buf[--i] = '0' + (int)(x % 10); x /= 10; } printf("%s", buf + i); }
int main(int argc, char **argv) {
  int m = atoi(argv[1]), Lmax = atoi(argv[2]);
  uint64_t p5[64]; p5[0] = 1; for (int i = 1; i < 30; i++) p5[i] = p5[i-1]*5;
  uint64_t MASK = (Lmax == 64) ? ~0ULL : ((1ULL << Lmax) - 1);
  /* inverse of 5^m mod 2^64 by Newton */
  uint64_t a5m = p5[m], inv = 1; for (int i = 0; i < 7; i++) inv *= 2 - a5m * inv;
  uint64_t negInv = (0 - inv);
  /* H at finest level, then folded levels stored contiguously: level L at offset (1<<L) */
  uint32_t *H = calloc((size_t)1 << (Lmax + 1), sizeof(uint32_t));
  if (!H) { fprintf(stderr, "alloc H\n"); return 1; }
  uint32_t *HL = H + ((size_t)1 << Lmax);
  /* DFS over q-process to depth m */
  uint64_t levelcount[64] = {0}, levelodd[64] = {0};
  typedef struct { uint64_t q; int i; } node;
  node *stk = malloc(sizeof(node) * 64 * 10);
  int sp = 0; stk[sp].q = 0; stk[sp].i = 0; sp++;
  while (sp) {
    node nd = stk[--sp];
    levelcount[nd.i]++; if (nd.q & 1) levelodd[nd.i]++;
    if (nd.i == m) { uint64_t r = (negInv * nd.q) & MASK; HL[r]++; continue; }
    uint64_t par = nd.q & 1, P = p5[nd.i];
    for (uint64_t a = (par ? 1 : 2); a <= 9; a += 2) { stk[sp].q = (nd.q + a * P) >> 1; stk[sp].i = nd.i + 1; sp++; }
  }
  for (int i = 0; i <= m; i++) printf("DFS Z_%d = %llu  odd=%llu\n", i, (unsigned long long)levelcount[i], (unsigned long long)levelodd[i]);
  fflush(stdout);
  /* fold: level L-1 from level L */
  for (int L = Lmax; L >= 1; L--) {
    uint32_t *src = H + ((size_t)1 << L), *dst = H + ((size_t)1 << (L - 1));
    size_t half = (size_t)1 << (L - 1);
    for (size_t s = 0; s < half; s++) dst[s] = src[s] + src[s + half];
  }
  /* B_L upward */
  uint64_t *Bp = calloc((size_t)1 << Lmax, sizeof(uint64_t)), *Bn = calloc((size_t)1 << Lmax, sizeof(uint64_t));
  if (!Bp || !Bn) { fprintf(stderr, "alloc B\n"); return 1; }
  Bp[0] = 1; /* L = 0 */
  for (int L = 1; L <= Lmax; L++) {
    size_t sz = (size_t)1 << L, szp = sz >> 1; uint64_t msk = sz - 1;
    memset(Bn, 0, sz * sizeof(uint64_t));
    for (size_t y = 0; y < szp; y++) { uint64_t v = Bp[y]; if (!v) continue; uint64_t base = (10 * (uint64_t)y);
      for (uint64_t a = 1; a <= 9; a++) Bn[(a + base) & msk] += v; }
    uint64_t *t = Bp; Bp = Bn; Bn = t;
    /* stats */
    uint64_t mx = 0, mn = ~0ULL; size_t amx = 0, amn = 0;
    for (size_t s = 0; s < sz; s++) { if (Bp[s] > mx) { mx = Bp[s]; amx = s; } if (Bp[s] < mn) { mn = Bp[s]; amn = s; } }
    uint32_t *Hl = H + sz;
    u128 Z = 0; for (size_t s = 0; s < sz; s++) Z += (u128)Hl[s] * Bp[s];
    printf("L=%d B_L[0]=Z_L=%llu maxB=%llu (at s=%zu) minB=%llu (at s=%zu)  Z_%d = ", L, (unsigned long long)Bp[0],
           (unsigned long long)mx, amx, (unsigned long long)mn, amn, m + L);
    print_u128(Z); printf("\n"); fflush(stdout);
  }
  return 0;
}
