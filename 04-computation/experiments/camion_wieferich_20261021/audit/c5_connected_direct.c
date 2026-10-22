/* c5_connected_direct.c -- second-pass audit: connected labelled graphs mod p^2 at n = p, by the rooted identity
 *   b' = T' * T^(-1)  (EGF; Session25 section 7 "Delta b = Delta T * T^(-1)"), i.e.
 *   W(0) = 1, W(n) = -sum_{k=1}^n C(n,k) T(k) W(n-k)        (Hurwitz inverse of T(n) = 2^C(n,2))
 *   b(n+1) = sum_{k=0}^n C(n,k) T(k+1) W(n-k)
 * This differs from cw_connected.c (which uses 2^C(n,2) = sum C(n-1,k-1) b(k) 2^C(n-k,2)).
 * Prints, for each prime p in [P0,P1]: b(p) mod p^2, the rotation congruence b(p) = 2^((p-1)/2) - 1 (mod p), and whether
 * p^2 | b(p) - 2^((p-1)/2) + 1.  Usage: ./c5c P0 P1      Second-pass auditor, 2026-10-21. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 M;
static u64 mulm(u64 x, u64 y){ return (u64)((u128)x * y % M); }
static u64 pw(u64 b, u64 e){ u64 r = 1 % M; b %= M; while(e){ if(e & 1) r = mulm(r, b); b = mulm(b, b); e >>= 1; } return r; }
static int isprime(u64 n){ if(n < 2) return 0; for(u64 d = 2; d * d <= n; d++) if(n % d == 0) return 0; return 1; }
int main(int argc, char **argv){
  u64 P0 = strtoull(argv[1], 0, 10), P1 = strtoull(argv[2], 0, 10); long np = 0, rotbad = 0, hits = 0;
  for(u64 p = (P0 < 3 ? 3 : P0); p <= P1; p++){
    if(!isprime(p)) continue;
    M = p * p; long N = (long)p;
    u64 *T = malloc((N + 2) * 8), *W = calloc(N + 2, 8), *row = calloc(N + 2, 8);
    u64 an = 1; T[0] = 1;
    for(long n = 0; n <= N; n++){ T[n + 1] = mulm(T[n], an); an = mulm(an, 2); }
    W[0] = 1; row[0] = 1;
    u64 bp = 0;
    for(long n = 0; n <= N - 1; n++){
      if(n >= 1){ for(long i = n; i >= 1; i--){ row[i] += row[i - 1]; if(row[i] >= M) row[i] -= M; } }   /* row = C(n,.) */
      if(n >= 1){ u128 acc = 0; for(long k = 1; k <= n; k++) acc += (u128)mulm(row[k], T[k]) * W[n - k]; W[n] = (M - (u64)(acc % M)) % M; }
      if(n == N - 1){ u128 acc = 0; for(long k = 0; k <= n; k++) acc += (u128)mulm(row[k], T[k + 1]) * W[n - k]; bp = (u64)(acc % M); }
    }
    u64 h = pw(2, (p - 1) / 2);
    u64 target = (h + M - 1) % M;
    int rot = ((bp + M - target) % M) % p == 0; if(!rot) rotbad++;
    int hit = (bp == target); if(hit) hits++;
    printf("p=%llu b(p) mod p^2=%llu  2^((p-1)/2)-1 mod p^2=%llu  rotation(mod p)=%s  %s\n", (unsigned long long)p, (unsigned long long)bp,
           (unsigned long long)target, rot ? "ok" : "FAIL", hit ? "CONNECTED-CW" : "");
    np++; free(T); free(W); free(row);
  }
  printf("SUMMARY primes=%ld rotation-congruence failures=%ld connected-CW=%ld\n", np, rotbad, hits);
  return 0;
}
