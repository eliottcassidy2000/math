/* camion_residue.c -- the Camion residue rho_C(p) = (eps*s_omega(p) - 1)/p mod p, where s_omega is the strong count with the
 * Teichmueller base omega = omega(2) (T(n) = omega^C(n,2), omega = 2^p mod p^2), and the decomposition of the Camion quotient
 *   q_C(p) := eps*D(p)/p = rho_C(p) - q_2(p)/2   (mod p)      [THM-4750 continuation: Lambda_p = sigma_p - <2>^(-1/2)].
 * Checks the decomposition against a direct D(p) mod p^2 and prints rho_C, q_2, the Wilson quotient w_p, and
 * B_(p-3) mod p (via the classical sum), for a fishing test of simple relations.
 * mac-mini-2026-10-21-fanopsp.  Usage: ./camion_residue P0 P1 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 pw(u64 b, u64 e, u64 m){ u64 r = 1 % m; b %= m; while(e){ if(e & 1) r = (u128)r * b % m; b = (u128)b * b % m; e >>= 1; } return r; }
static int isprime(u64 n){ if(n < 2) return 0; for(u64 d = 2; d * d <= n; d++) if(n % d == 0) return 0; return 1; }
static u64 strong_at(u64 base, long N, u64 Q){   /* s(N) mod Q for T(n) = base^C(n,2) */
  u64 *T = malloc((N + 1) * 8), *s = calloc(N + 1, 8), *row = calloc(N + 2, 8);
  u64 PH = 0; (void)PH;
  for(long j = 0; j <= N; j++) T[j] = pw(base, (u64)j * (j - 1) / 2, Q);
  row[0] = 1;
  for(long n = 1; n <= N; n++){
    for(long i = n - 1; i >= 1; i--){ row[i] += row[i - 1]; if(row[i] >= Q) row[i] -= Q; }
    row[n] = 1; u128 acc = 0;
    for(long k = 1; k < n; k++) acc += (u128)((u128)row[k] * s[k] % Q) * T[n - k];
    s[n] = (T[n] + Q - (u64)(acc % Q)) % Q; }
  u64 r = s[N]; free(T); free(s); free(row); return r; }
int main(int argc, char **argv){
  u64 P0 = atoll(argv[1]), P1 = atoll(argv[2]); long bad = 0, np = 0;
  for(u64 p = (P0 < 5 ? 5 : P0); p <= P1; p++){ if(!isprime(p)) continue;
    u64 Q = p * p; int eps = (p % 8 == 1 || p % 8 == 7) ? 1 : -1;
    u64 om = pw(2, p, Q);                                   /* Teichmueller lift of 2 mod p^2 */
    u64 so = strong_at(om, (long)p, Q), s2 = strong_at(2, (long)p, Q);
    u64 sig = eps == 1 ? so : (Q - so) % Q;                 /* eps * s_omega(p) */
    u64 rhoC = ((sig + Q - 1) % Q) / p;                      /* (sigma - 1)/p mod p */
    u64 c = pw(2, (p - 1) / 2, Q); u64 D = (s2 + Q - c) % Q; u64 epsD = eps == 1 ? D : (Q - D) % Q;
    u64 qC = epsD / p;                                       /* eps D(p)/p mod p */
    u64 q2 = ((pw(2, p - 1, Q) + Q - 1) % Q) / p;
    u64 inv2 = (p + 1) / 2;
    u64 pred = (rhoC + p - q2 * inv2 % p) % p;
    if(pred != qC % p) bad++;
    u64 f = 1; for(u64 i = 2; i < p; i++) f = (u128)f * i % Q; u64 wp = ((f + 1) % Q) / p;   /* Wilson quotient */
    /* B_(p-3) mod p via sum_{a=1}^{p-1} a^(p-3) ... use the classical congruence sum a^k = p B_k mod p^2 for k = p-3 (even, p >= 5) */
    u64 S = 0; for(u64 a = 1; a < p; a++) S = (S + pw(a, p - 3, Q)) % Q; u64 Bp3 = S / p % p;  /* = B_(p-3) mod p (Bernoulli, up to the standard congruence) */
    printf("p=%llu eps=%+d rho_C=%llu q_2=%llu q_C=%llu w_p=%llu B_(p-3)=%llu %s\n", (unsigned long long)p, eps, (unsigned long long)rhoC,
           (unsigned long long)q2, (unsigned long long)(qC % p), (unsigned long long)wp, (unsigned long long)Bp3, pred == qC % p ? "" : "DECOMP-FAIL");
    np++; }
  printf("SUMMARY primes=%ld decomposition failures=%ld\n", np, bad); return 0; }
