/* hyper_r_second_order.c -- "degree = uniformity" (THM-4750 (7), HYP-9361 (d)) for connected r-uniform hypergraphs, general r.
 * b(n): 2^C(n,r) = sum_k C(n-1,k-1) b(k) 2^C(n-k,r).  At primes with gcd(ord_p(2), r!) = 1 the class sequences
 * A(l) = b(m + l(p-1)) are, mod p^2, polynomials of degree <= r in l (a_1..a_r divisible by p, a_(r+1) = 0 mod p^2).
 * Zero classes: count roots in F_p of g(l) = sum_i (a_i/p) C(l,i); compare with a random degree-r polynomial
 * (root counts distributed as fixed points of a random permutation in S_r).
 * mac-mini-2026-10-21-fanopsp.  Usage: ./hyper_r_second_order r P0 P1 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 pw(u64 b, u64 e, u64 m){ u64 r = 1 % m; b %= m; while(e){ if(e & 1) r = (u128)r * b % m; b = (u128)b * b % m; e >>= 1; } return r; }
static int isprime(u64 n){ if(n < 2) return 0; for(u64 d = 2; d * d <= n; d++) if(n % d == 0) return 0; return 1; }
static u64 ord2(u64 p){ u64 o = 1, x = 2 % p; while(x != 1){ x = x * 2 % p; o++; } return o; }
static u64 gcdu(u64 a, u64 b){ while(b){ u64 t = a % b; a = b; b = t; } return a; }
static u128 binom_u128(u64 n, int r){ u128 c = 1; for(int i = 0; i < r; i++){ c = c * (n - i) / (i + 1); } return c; } /* exact for n >= r */
int main(int argc, char **argv){
  int r = atoi(argv[1]); u64 P0 = atoll(argv[2]), P1 = atoll(argv[3]);
  u64 rf = 1; for(int i = 2; i <= r; i++) rf *= i;
  long hist[16] = {0}, nprimes = 0, nzc = 0, bad = 0, nclasses = 0, notfull = 0;
  for(u64 p = (P0 < 5 ? 5 : P0); p <= P1; p++){
    if(!isprime(p) || p <= (u64)r + 1) continue;
    u64 t = ord2(p); if(gcdu(t, rf) != 1) continue;
    u64 Q = p * p, PH = p * (p - 1);
    long m0 = r + 3, m1 = m0 + (long)p - 2; long M = m1 + (long)(r + 1) * (long)(p - 1);
    u64 *T = malloc((M + 1) * 8), *b = calloc(M + 1, 8), *row = calloc(M + 2, 8);
    for(long n = 0; n <= M; n++){ if(n < r){ T[n] = 1; continue; } u128 c = binom_u128((u64)n, r); T[n] = pw(2, (u64)(c % PH) + PH, Q); }
    row[0] = 1;
    for(long n = 1; n <= M; n++){
      if(n >= 2){ for(long i = n - 1; i >= 1; i--){ row[i] += row[i - 1]; if(row[i] >= Q) row[i] -= Q; } }
      u128 acc = 0; for(long k = 1; k < n; k++) acc += (u128)((u128)row[k - 1] * b[k] % Q) * T[n - k];
      b[n] = (T[n] + Q - (u64)(acc % Q)) % Q;
    }
    nprimes++;
    char buf[4096]; int bl = 0; buf[0] = 0;
    for(long m = m0; m <= m1; m++){
      u64 A[20], d[20]; int R = r + 2;
      for(int l = 0; l < R; l++) A[l] = b[m + (long)l * (long)(p - 1)];
      for(int i = 0; i < R; i++){ d[i] = A[0]; for(int l = 0; l < R - 1 - i; l++) A[l] = (A[l + 1] + Q - A[l]) % Q; }
      nclasses++;
      int ok = 1; for(int i = 1; i <= r; i++) if(d[i] % p) ok = 0; if(d[r + 1] % Q) ok = 0;
      if(!ok){ bad++; continue; }
      if(d[r] % Q == 0) notfull++;
      if(d[0] % p) continue;
      nzc++;
      u64 c[20]; for(int i = 0; i <= r; i++) c[i] = d[i] / p % p;
      u64 invf[20]; invf[0] = 1; for(int i = 1; i <= r; i++) invf[i] = invf[i - 1] * pw(i, p - 2, p) % p;
      int roots = 0;
      for(u64 l = 0; l < p; l++){ u64 v = 0, fall = 1;   /* fall = l(l-1)...(l-i+1) */
        for(int i = 0; i <= r; i++){ v = (v + c[i] * (fall % p) % p * invf[i]) % p; fall = fall * ((l + p - (u64)i) % p) % p; }
        if(v == 0) roots++; }
      hist[roots < 15 ? roots : 15]++;
      if(bl < 3900) bl += sprintf(buf + bl, " %ld:%d", m % (long)(p - 1), roots);
    }
    printf("p=%llu ord2=%llu zero-classes:%s\n", (unsigned long long)p, (unsigned long long)t, buf); fflush(stdout);
    free(T); free(b); free(row);
  }
  printf("SUMMARY r=%d primes=%ld classes=%ld degree-r law failures=%ld classes with p^2|a_r=%ld zero classes=%ld\n", r, nprimes, nclasses, bad, notfull, nzc);
  printf("root-count histogram:"); for(int i = 0; i <= r; i++) printf(" %d:%ld", i, hist[i]); printf("\n");
  return 0;
}
