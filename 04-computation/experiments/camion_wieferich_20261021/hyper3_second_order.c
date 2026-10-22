/* hyper3_second_order.c -- the "degree = uniformity" test (HYP-9361 (d)) on connected 3-uniform hypergraphs.
 *
 * b(n) = connected labelled 3-uniform hypergraphs on n vertices:  T3(n) = 2^C(n,3) = sum_k C(n-1,k-1) b(k) T3(n-k).
 * For primes p with gcd(ord_p(2), 6) = 1 the Teichmueller part w^C(n,3) has period dividing p-1, so (as for tournaments,
 * but with multiplier 1) the class sequences A(l) = b(m + l(p-1)) should be p-adically interpolable at speed 1/3:
 * mod p^2 a CUBIC polynomial in l (a_1, a_2, a_3 divisible by p, a_4 = 0 mod p^2).  In a zero class (p | b(m)) the
 * l mod p with p^2 | b(m + l(p-1)) are then the roots of a cubic over F_p; for random cubics the root counts 0, 1, 3
 * have frequencies 1/3, 1/2, 1/6 (2 only for repeated roots).  Tournaments (C(n,2)) give quadratics: 0 or 2.
 * mac-mini-2026-10-21-fanopsp.  Usage: ./hyper3_second_order P0 P1 [VMAX] */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 pw(u64 b, u64 e, u64 m){ u64 r = 1 % m; b %= m; while(e){ if(e & 1) r = (u128)r * b % m; b = (u128)b * b % m; e >>= 1; } return r; }
static int isprime(u64 n){ if(n < 2) return 0; for(u64 d = 2; d * d <= n; d++) if(n % d == 0) return 0; return 1; }
static u64 ord2(u64 p){ u64 o = 1, x = 2 % p; while(x != 1){ x = x * 2 % p; o++; } return o; }
static u64 gcdu(u64 a, u64 b){ while(b){ u64 t = a % b; a = b; b = t; } return a; }
int main(int argc, char **argv){
  u64 P0 = atoll(argv[1]), P1 = atoll(argv[2]); u64 VMAX = argc > 3 ? atoll(argv[3]) : 0;
  long hist[5] = {0}, nprimes = 0, nzc = 0, bad = 0, cubic_zero = 0, nclasses = 0, direct = 0, direct_bad = 0;
  for(u64 p = (P0 < 5 ? 5 : P0); p <= P1; p++){
    if(!isprime(p)) continue;
    u64 t = ord2(p); if(gcdu(t, 6) != 1) continue;
    u64 Q = p * p, PH = p * (p - 1);
    long m0 = 6, m1 = m0 + (long)p - 2; long M = m1 + 4 * (long)(p - 1);
    if(VMAX && p <= VMAX) M = m1 + (long)(p + 2) * (long)(p - 1);
    u64 *T = malloc((M + 1) * 8), *b = calloc(M + 1, 8), *row = calloc(M + 2, 8);
    for(long n = 0; n <= M; n++){ u128 c3 = (u128)n * (n - 1) * (n - 2) / 6; T[n] = pw(2, (u64)(c3 % PH) + (c3 >= 3 ? PH : 0), Q); }
    T[0] = 1; T[1] = 1; T[2] = 1;
    row[0] = 1;   /* row = C(n-1, .) */
    for(long n = 1; n <= M; n++){
      if(n >= 2){ for(long i = n - 1; i >= 1; i--){ row[i] += row[i - 1]; if(row[i] >= Q) row[i] -= Q; } }
      u128 acc = 0; for(long k = 1; k < n; k++) acc += (u128)((u128)row[k - 1] * b[k] % Q) * T[n - k];
      b[n] = (T[n] + Q - (u64)(acc % Q)) % Q;
    }
    nprimes++;
    int zc = 0; int rootsum = 0; char buf[2048]; int bl = 0; buf[0] = 0;
    u64 inv2 = (Q + 1) / 2 % p, inv6 = pw(6, p - 2, p);
    for(long m = m0; m <= m1; m++){
      u64 A[5]; for(int l = 0; l < 5; l++) A[l] = b[m + (long)l * (long)(p - 1)];
      u64 a1 = (A[1] + Q - A[0]) % Q, a2 = (A[2] + 2 * Q - 2 * A[1] % Q + A[0]) % Q;
      u64 a3 = (A[3] + 3 * Q - 3 * A[2] % Q + 3 * A[1] % Q + Q - A[0]) % Q;
      u64 a4 = (A[4] + 4 * Q - 4 * A[3] % Q + 6 * A[2] % Q + 4 * Q - 4 * A[1] % Q + A[0]) % Q;
      nclasses++;
      if(a1 % p || a2 % p || a3 % p || a4 % Q){ bad++; continue; }
      if(a3 % Q == 0) cubic_zero++;
      if(VMAX && p <= VMAX) for(long l = 0; l <= (long)p + 1; l++){ long n = m + l * (long)(p - 1); if(n > M) break;
        u64 L = (u64)l; u64 c2 = (u64)(((long)l * (l - 1) / 2) % (long)Q), c3v = (u64)(((__int128)l * (l - 1) * (l - 2) / 6) % (__int128)Q);
        u64 pr = (A[0] + (u128)a1 * L % Q + (u128)a2 * c2 % Q + (u128)a3 * c3v % Q) % Q; direct++; if(pr != b[n]) direct_bad++; }
      if(A[0] % p) continue;
      zc++; nzc++;
      u64 c0 = A[0] / p, c1 = a1 / p, c2 = a2 / p, c3 = a3 / p;
      int roots = 0;
      for(u64 l = 0; l < p; l++){ u64 C2 = l * ((l + p - 1) % p) % p * inv2 % p; u64 C3 = C2 * ((l + p - 2) % p) % p * pw(3, p - 2, p) % p;
        u64 v = (c0 + c1 * l + c2 * C2 + c3 * C3) % p; if(v == 0) roots++; }
      (void)inv6; hist[roots < 4 ? roots : 4]++; rootsum += roots;
      if(bl < 1900) bl += sprintf(buf + bl, " %ld:%d%s", m % (long)(p - 1), roots, c3 % p ? "" : "(deg<3)");
    }
    printf("p=%llu ord2=%llu zero-classes:%s\n", (unsigned long long)p, (unsigned long long)t, buf); fflush(stdout);
    free(T); free(b); free(row);
  }
  printf("SUMMARY primes=%ld classes=%ld  cubic-law failures (p|a1,a2,a3; p^2|a4)=%ld  classes with p^2|a3 (not cubic)=%ld  direct checks=%ld failures=%ld\n",
         nprimes, nclasses, bad, cubic_zero, direct, direct_bad);
  printf("zero classes=%ld  root-count histogram: 0:%ld 1:%ld 2:%ld 3:%ld >3:%ld   (random cubic: 1/3, 1/2, ~0, 1/6)\n", nzc, hist[0], hist[1], hist[2], hist[3], hist[4]);
  return 0;
}
