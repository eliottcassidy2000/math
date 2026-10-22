/* c7_hyper_r.c -- second-pass audit of THM-4750 section 6b (r-uniform hypergraphs, r = 4 table), independent code.
 *  * 2^C(n,r) mod p^2 WITHOUT big exponents: P_j(n) = 2^C(n,j) satisfies P_j(n+1) = P_j(n) * P_(j-1)(n), P_0 = 2
 *    (Pascal's rule in the exponent); the session used C(n,r) mod p(p-1) + PH.
 *  * connected counts by the rooted identity b' = T' * T^(-1): W = Hurwitz inverse of T, b(n+1) = sum_k C(n,k) T(k+1) W(n-k)
 *    (the session used 2^C(n,r) = sum C(n-1,k-1) b(k) 2^C(n-k,r)).
 *  * prime set: 7 <= p, gcd(ord_p(2), r!) = 1 (ord by brute force).
 * Output per prime: one line "P p ord lawfail notfull", then for every class m in [r+3, r+p+1] that is a zero class a line
 * "Z m (m mod p-1) c_0 .. c_r" with c_i = a_i/p mod p; then lines "B m law(0/1) g(-1)" for the low bases m = 2, 3, i.e.
 * whether l -> b(m + l(p-1)) already obeys the degree-r law mod p^2 at base m, and the value at l = -1 of the class
 * polynomial of base m + (p-1) (mod p, after division by p; meaningful when the class is a zero class).
 * Usage: ./c7 r P0 P1        Second-pass auditor, 2026-10-21. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 M;
static u64 mulm(u64 x, u64 y){ return (u64)((u128)x * y % M); }
static int isprime(u64 n){ if(n < 2) return 0; for(u64 d = 2; d * d <= n; d++) if(n % d == 0) return 0; return 1; }
static u64 gcdu(u64 a, u64 b){ while(b){ u64 t = a % b; a = b; b = t; } return a; }
static u64 inv_p(u64 a, u64 p){ u64 r = 1, b = a % p, e = p - 2; while(e){ if(e & 1) r = r * b % p; b = b * b % p; e >>= 1; } return r; }
int main(int argc, char **argv){
  int r = atoi(argv[1]); u64 P0 = strtoull(argv[2], 0, 10), P1 = strtoull(argv[3], 0, 10);
  u64 rf = 1; for(int i = 2; i <= r; i++) rf *= (u64)i;
  for(u64 p = (P0 < 7 ? 7 : P0); p <= P1; p++){
    if(!isprime(p) || p <= (u64)r + 1) continue;
    u64 t = 1, x = 2 % p; while(x != 1){ x = x * 2 % p; t++; }
    if(gcdu(t, rf) != 1) continue;
    M = p * p;
    long mlo = 2, mhi = r + (long)p + 1;
    long Nmax = mhi + (long)(r + 1) * (long)(p - 1) + 2;
    u64 *P = malloc((size_t)(r + 1) * 8), *T = malloc((Nmax + 2) * 8), *W = calloc(Nmax + 2, 8), *b = calloc(Nmax + 2, 8), *row = calloc(Nmax + 2, 8);
    P[0] = 2; for(int j = 1; j <= r; j++) P[j] = 1;           /* P_j(0) = 2^C(0,j) */
    for(long n = 0; n <= Nmax + 1; n++){
      T[n] = P[r];
      for(int j = r; j >= 1; j--) P[j] = mulm(P[j], P[j - 1]); /* P_j(n+1) = P_j(n) P_(j-1)(n); P_0 stays 2 */
    }
    W[0] = 1; row[0] = 1;
    for(long n = 0; n < Nmax; n++){
      if(n >= 1){ for(long i = n; i >= 1; i--){ row[i] += row[i - 1]; if(row[i] >= M) row[i] -= M; } }
      if(n >= 1){ u128 acc = 0; for(long k = 1; k <= n; k++) acc += (u128)mulm(row[k], T[k]) * W[n - k]; W[n] = (M - (u64)(acc % M)) % M; }
      u128 acc = 0; for(long k = 0; k <= n; k++) acc += (u128)mulm(row[k], T[k + 1]) * W[n - k];
      b[n + 1] = (u64)(acc % M);
    }
    long lawfail = 0, notfull = 0;
    char *buf = malloc(1 << 16); int bl = 0;
    for(long m = mlo; m <= mhi; m++){
      u64 A[24], a[24]; int R = r + 2;
      for(int l = 0; l < R; l++) A[l] = b[m + (long)l * (long)(p - 1)];
      for(int i = 0; i < R; i++){ a[i] = A[0]; for(int l = 0; l < R - 1 - i; l++) A[l] = (A[l + 1] + M - A[l]) % M; }
      int law = 1; for(int i = 1; i <= r; i++) if(a[i] % p) law = 0; if(a[r + 1] % M) law = 0;
      if(m <= 3){                                              /* low bases, below the proved threshold r + 2 */
        /* class polynomial of base m + (p-1) evaluated at l = -1 (i.e. at n = m), divided by p, mod p; needs a zero class */
        long m2 = m + (long)(p - 1); u64 B[24], c[24];
        for(int l = 0; l < R; l++) B[l] = b[m2 + (long)l * (long)(p - 1)];
        for(int i = 0; i < R; i++){ c[i] = B[0]; for(int l = 0; l < R - 1 - i; l++) B[l] = (B[l + 1] + M - B[l]) % M; }
        int zc = (c[0] % p) == 0;
        /* sum_i c_i C(-1, i) = sum_i (-1)^i c_i  (mod p^2) */
        u64 g = 0; for(int i = 0; i <= r; i++) g = (i % 2) ? (g + M - c[i] % M) % M : (g + c[i]) % M;
        bl += sprintf(buf + bl, "B %ld law=%d zeroclass=%d gminus1_over_p=%llu value_mod_p2=%llu b(m)=%llu\n", m, law, zc,
                      (unsigned long long)((g % p) == 0 ? (g / p) % p : 999999), (unsigned long long)g, (unsigned long long)b[m]);
        continue;
      }
      if(m < r + 3) continue;
      if(!law){ lawfail++; continue; }
      if(a[r] % M == 0) notfull++;
      if(a[0] % p) continue;
      bl += sprintf(buf + bl, "Z %ld %ld", m, m % (long)(p - 1));
      for(int i = 0; i <= r; i++) bl += sprintf(buf + bl, " %llu", (unsigned long long)(a[i] / p % p));
      bl += sprintf(buf + bl, "\n");
    }
    printf("P %llu %llu %ld %ld\n%s", (unsigned long long)p, (unsigned long long)t, lawfail, notfull, buf); fflush(stdout);
    free(buf); free(P); free(T); free(W); free(b); free(row);
  }
  (void)inv_p;
  return 0;
}
