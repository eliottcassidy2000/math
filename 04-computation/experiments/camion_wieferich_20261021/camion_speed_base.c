/* camion_speed_base.c (from camion_speed.c) -- p-adic speed of the Camion functions.
 *
 * For a prime p and K >= 1, computes s(n) mod p^K (Wright's recurrence), and for every base m in
 * [2K-1, 2K-1 + p-2] (one per residue class mod p-1) the Mahler coefficients
 *   a_i(m) = Delta^i A(0),  A(l) = eps^-l D(m + l(p-1)),  i = 0..2K-1.
 * Prints, for each i, the minimum over all classes of v_p(a_i(m)) (capped at K), and how many classes attain it.
 * Theory (THM-4750 draft): v_p(a_i) >= ceil(i/2) always (half speed); if p is a Wieferich prime
 * (2^(p-1) = 1 mod p^2) then v_p(a_i) >= i (full speed).  The minima test whether these bounds are attained.
 * Generalised to base a (T(n) = a^C(n,2); the circulant correction is used only for a = 2; s itself obeys the same laws):
 * the curvature is then the Fermat quotient q_a(p), so base-a Wieferich primes (a^(p-1) = 1 mod p^2) give full speed.
 * mac-mini-2026-10-21-fanopsp.  Usage: ./camion_speed_base p K a */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 Q;
static inline u64 mulm(u64 a, u64 b){ return (u64)((u128)a * b % Q); }
static u64 pw(u64 b, u64 e){ u64 r = 1 % Q; b %= Q; while(e){ if(e & 1) r = mulm(r, b); b = mulm(b, b); e >>= 1; } return r; }
static int vp(u64 x, u64 p, int K){ if(x == 0) return K; int c = 0; while(x % p == 0){ x /= p; c++; } return c < K ? c : K; }
int main(int argc, char **argv){ u64 BASE = argc > 3 ? atoll(argv[3]) : 2;
  u64 p = atoll(argv[1]); int K = atoi(argv[2]);
  Q = 1; for(int i = 0; i < K; i++) Q *= p;
  int eps = (pw(BASE % p, (p - 1) / 2) % p == 1) ? 1 : -1;  /* (BASE/p) by Euler (pw works mod p^K; reduce mod p) */ int R = 2 * K;           /* coefficients i = 0..2K-1 */
  long m0 = 2 * K - 1, m1 = m0 + (long)p - 2; long M = m1 + (long)(R - 1) * (long)(p - 1);
  u64 *T = malloc((M + 1) * 8), *s = calloc(M + 1, 8), *row = calloc(M + 2, 8);
  for(long j = 0; j <= M; j++) T[j] = pw(BASE, (u64)j * (j - 1) / 2);
  row[0] = 1;
  for(long n = 1; n <= M; n++){
    for(long i = n - 1; i >= 1; i--){ row[i] += row[i - 1]; if(row[i] >= Q) row[i] -= Q; }
    row[n] = 1;
    u128 acc = 0;
    for(long k = 1; k < n; k++){ acc += (u128)mulm(row[k], s[k]) * T[n - k]; if((k & 63) == 0) acc %= Q; }
    s[n] = (T[n] + Q - (u64)(acc % Q)) % Q;
  }
  u64 q2 = (pw(BASE, p - 1) + Q - 1) % Q; int w = vp(q2, p, K) - 1;   /* v_p(2^(p-1) - 1) - 1 = v_p(q_2) */
  int minv[64]; long att[64]; for(int i = 0; i < R; i++){ minv[i] = K; att[i] = 0; }
  for(long m = m0; m <= m1; m++){
    u64 A[64], t[64];
    for(int l = 0; l < R; l++){ long n = m + (long)l * (long)(p - 1);
      u64 c = (BASE == 2 && (n & 1)) ? pw(2, (u64)(n - 1) / 2) : 0; u64 D = (s[n] + Q - c) % Q;
      A[l] = (eps == 1 || l % 2 == 0) ? D : (Q - D) % Q; }
    for(int l = 0; l < R; l++) t[l] = A[l];
    for(int i = 0; i < R; i++){ int v = vp(t[0], p, K);
      if(v < minv[i]){ minv[i] = v; att[i] = 1; } else if(v == minv[i]) att[i]++;
      for(int l = 0; l < R - 1 - i; l++) t[l] = (t[l + 1] + Q - t[l]) % Q; }
  }
  printf("base=%llu p=%llu K=%d eps=%+d v_p(q_a)=%s%d classes=%ld (m=%ld..%ld, s to n=%ld)\n", (unsigned long long)BASE, (unsigned long long)p, K, eps,
         (w >= K - 1) ? ">=" : "", w, m1 - m0 + 1, m0, m1, M);
  printf("  i      : "); for(int i = 0; i < R; i++) printf("%3d", i); printf("\n");
  printf("  min v  : "); for(int i = 0; i < R; i++) printf("%3d", minv[i]); printf("\n");
  printf("  ceil i/2:"); for(int i = 0; i < R; i++){ int c = (i + 1) / 2; printf("%3d", c < K ? c : K); } printf("\n");
  printf("  #attain: "); for(int i = 0; i < R; i++) printf(" %ld", att[i]); printf("\n");
  return 0;
}
