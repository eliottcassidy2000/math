/* camion_second_order.c -- second-order (mod p^2) structure of the Camion function.
 *
 * Theory (THM-4750 draft): for m >= 3 the sequence A(l) = eps^-l D(m + l(p-1)) satisfies, mod p^2,
 *   A(l) = a0 + a1 l + a2 C(l,2)   with p | a1, p | a2   (higher Mahler coefficients vanish mod p^2),
 * and the half-speed coefficient is a Fermat quotient times a tournament statistic:
 *   a2/p = -q_2(p) X(m-2) (mod p),   X = Delta^2 T * T^(-2) (Hurwitz product), i.e. EGF T''/T^2.
 * For a zero class (p | a0) the reduction g(l) = a0/p + (a1/p) l + (a2/p) C(l,2) is a quadratic over F_p,
 * whose roots are the l mod p with p^2 | D(m + l(p-1)).  Prediction: 0 or 2 roots (about 1/2 each),
 * except at Wieferich primes, where g is linear (1 root).
 *
 * For each prime p in [P0, P1]: s(n) mod p^2 for n <= 3p (Wright's recurrence), all classes m = 3..p+1,
 * checks p | a1, p | a2 and the a2 formula, counts roots of g for every zero class, and (if p <= VMAX)
 * verifies the quadratic law directly for l = 0..p+1.
 * Output: one line per prime with zero classes and root counts; summary histogram.
 * mac-mini-2026-10-21-fanopsp.  Usage: ./camion_second_order P0 P1 [VMAX] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 Q, P;
static u64 pw(u64 b, u64 e, u64 m){ u64 r = 1 % m; b %= m; while(e){ if(e & 1) r = (u128)r * b % m; b = (u128)b * b % m; e >>= 1; } return r; }
static int isprime(u64 n){ if(n < 2) return 0; for(u64 d = 2; d * d <= n; d++) if(n % d == 0) return 0; return 1; }
/* s(n) mod q for n <= M by Wright's recurrence T(n) = sum_{k=1}^n C(n,k) s(k) T(n-k) */
static u64 *s_upto(long M, u64 q){
  u64 *T = malloc((M + 1) * 8), *s = calloc(M + 1, 8), *row = calloc(M + 2, 8);
  for(long j = 0; j <= M; j++) T[j] = pw(2, (u64)j * (j - 1) / 2, q);
  row[0] = 1;
  for(long n = 1; n <= M; n++){
    for(long i = n - 1; i >= 1; i--){ row[i] += row[i - 1]; if(row[i] >= q) row[i] -= q; }
    row[n] = 1;
    u128 acc = 0;
    for(long k = 1; k < n; k++) acc += (u128)((u128)row[k] * s[k] % q) * T[n - k];
    s[n] = (T[n] + q - (u64)(acc % q)) % q;
  }
  free(T); free(row); return s;
}
static u64 cfun(long n, u64 q){ return (n & 1) ? pw(2, (n - 1) / 2, q) : 0; }
int main(int argc, char **argv){
  u64 P0 = atoll(argv[1]), P1 = atoll(argv[2]); u64 VMAX = (argc > 3) ? atoll(argv[3]) : 0;
  long hist[8] = {0}, histW[8] = {0}; long nzc = 0, nprimes = 0, bad = 0, formula_bad = 0, direct_bad = 0, direct_checked = 0;
  long lam_even = 0, lam_odd = 0, nlin = 0, ndbl = 0;
  for(P = (P0 < 3 ? 3 : P0); P <= P1; P++){
    if(!isprime(P)) continue;
    u64 p = P; Q = p * p; int eps = (p % 8 == 1 || p % 8 == 7) ? 1 : -1;
    long M = 3 * (long)p + 2; if(VMAX && p <= VMAX) M = (long)(p + 1) + (long)(p + 2) * (long)(p - 1) + 2;
    u64 *s = s_upto(M, Q);
    /* X(n) mod p for n <= p: Y = (e - S)^2, X = Delta^2 T * Y (binomial convolutions) */
    long NX = p + 2; u64 *E = malloc((NX + 1) * 8), *Y = calloc(NX + 1, 8), *X = calloc(NX + 1, 8), *T2 = malloc((NX + 3) * 8), *rw = calloc(NX + 2, 8);
    E[0] = 1; for(long n = 1; n <= NX; n++) E[n] = (p - s[n] % p) % p;
    for(long n = 0; n <= NX + 2; n++) T2[n] = pw(2, (u64)n * (n - 1) / 2, p);
    rw[0] = 1;   /* rolling Pascal row mod p: Y = (e - S)^2, X = Delta^2 T * Y */
    for(long n = 0; n <= NX; n++){
      if(n){ for(long i = n - 1; i >= 1; i--) rw[i] = (rw[i] + rw[i - 1]) % p; rw[n] = 1; }
      u64 a = 0; for(long k = 0; k <= n; k++) a = (a + rw[k] * E[k] % p * E[n - k]) % p; Y[n] = a;
      a = 0; for(long k = 0; k <= n; k++) a = (a + rw[k] * T2[k + 2] % p * Y[n - k]) % p; X[n] = a; }
    u64 q2 = ((pw(2, p - 1, Q) + Q - 1) % Q) / p;   /* Fermat quotient mod p */
    int wief = (q2 % p == 0);
    printf("p=%llu eps=%+d q2=%llu%s zero-classes:", (unsigned long long)p, eps, (unsigned long long)q2, wief ? " WIEFERICH" : "");
    int lam = 0;
    for(long m = 3; m <= (long)p + 1; m++){
      u64 A[3];
      for(int l = 0; l < 3; l++){ long n = m + (long)l * (long)(p - 1); u64 D = (s[n] + Q - cfun(n, Q)) % Q; A[l] = (eps == 1 || l % 2 == 0) ? D : (Q - D) % Q; }
      u64 a0 = A[0], a1 = (A[1] + Q - A[0]) % Q, a2 = (A[2] + 2 * Q - 2 * A[1] % Q + A[0]) % Q;
      if(a1 % p || a2 % p){ bad++; printf(" [BAD m=%ld]", m); continue; }
      u64 c2 = a2 / p;
      u64 pred = (p - q2 % p) % p * X[m - 2] % p;
      if(c2 % p != pred){ formula_bad++; }
      if(VMAX && p <= VMAX){ /* direct check of the quadratic law mod p^2 for l = 0..p+1 */
        for(long l = 0; l <= (long)p + 1; l++){ long n = m + l * (long)(p - 1); if(n > M) break;
          u64 D = (s[n] + Q - cfun(n, Q)) % Q; u64 Al = (eps == 1 || l % 2 == 0) ? D : (Q - D) % Q;
          u64 Cl2 = (u64)((l * (l - 1) / 2) % (long)Q);
          u64 pr = (a0 + (u128)a1 * (l % Q) % Q + (u128)a2 * Cl2 % Q) % Q;
          direct_checked++; if(pr != Al) direct_bad++; } }
      if(a0 % p) continue;                      /* not a zero class */
      nzc++;
      u64 c0 = a0 / p, c1 = a1 / p; u64 inv2 = (p + 1) / 2;
      int roots = 0; for(u64 l = 0; l < p; l++){ u64 v = (c0 + c1 * l % p + c2 * (l * ((l + p - 1) % p) % p) % p * inv2) % p; if(v == 0) roots++; }
      /* g = (c2/2) l^2 + (c1 - c2/2) l + c0: tag L if linear (c2 = 0), D if a double root (discriminant 0) */
      u64 A2 = c2 * inv2 % p, B1 = (c1 + p - A2) % p; u64 disc = (B1 * B1 % p + p - 4 * A2 % p * c0 % p) % p;
      const char *tag = (c2 % p == 0) ? "L" : (disc == 0 ? "D" : "");
      int cls = (int)(m % (p - 1)); if(cls == 0) cls = (int)(p - 1);
      printf(" %d:%d%s%s", cls, roots, tag, (m == (long)p && c0 == 0) ? "CW" : "");
      if(c2 % p == 0) nlin++; if(c2 % p && disc == 0) ndbl++;
      int hb = roots < 7 ? roots : 7; if(wief) histW[hb]++; else hist[hb]++;
      lam += roots;
    }
    if(lam % 2) lam_odd++; else lam_even++;
    printf("  lambda_C=%d\n", lam); fflush(stdout);
    nprimes++;
    free(rw); free(E); free(Y); free(X); free(T2); free(s);
  }
  printf("SUMMARY primes=%ld zero-classes=%ld  checks: p|a1,p|a2 failures=%ld  a2-formula failures=%ld  direct-law checks=%ld failures=%ld\n",
         nprimes, nzc, bad, formula_bad, direct_checked, direct_bad);
  printf("root-count histogram (non-Wieferich):"); for(int i = 0; i < 8; i++) printf(" %d:%ld", i, hist[i]); printf("\n");
  printf("root-count histogram (Wieferich):"); for(int i = 0; i < 8; i++) printf(" %d:%ld", i, histW[i]); printf("\n");
  printf("lambda_C parity: even=%ld odd=%ld   linear classes (c2 = 0): %ld   double roots: %ld\n", lam_even, lam_odd, nlin, ndbl);
  return 0;
}
