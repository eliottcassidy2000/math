/* cw_connected.c -- "Camion-Wieferich" primes for connected labelled graphs: the rotation argument gives
 * cg(p) = #(nonempty circulant graphs on Z_p) = 2^((p-1)/2) - 1 (mod p); test p^2 | cg(p) - 2^((p-1)/2) + 1.
 * cg from 2^C(n,2) = sum_k C(n-1,k-1) cg(k) 2^C(n-k,2).  Compare with the tournament list {3, 13, 22409}.
 * mac-mini-2026-10-21-fanopsp.  Usage: ./cw_connected P0 P1 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 pw(u64 b, u64 e, u64 m){ u64 r = 1 % m; b %= m; while(e){ if(e & 1) r = (u128)r * b % m; b = (u128)b * b % m; e >>= 1; } return r; }
static int isprime(u64 n){ if(n < 2) return 0; for(u64 d = 2; d * d <= n; d++) if(n % d == 0) return 0; return 1; }
int main(int argc, char **argv){
  u64 P0 = atoll(argv[1]), P1 = atoll(argv[2]); long hits = 0, np = 0, eulerbad = 0;
  for(u64 p = (P0 < 3 ? 3 : P0); p <= P1; p++){ if(!isprime(p)) continue;
    u64 Q = p * p; long N = (long)p;
    u64 *T = malloc((N + 1) * 8), *b = calloc(N + 1, 8), *row = calloc(N + 2, 8);
    for(long j = 0; j <= N; j++) T[j] = pw(2, (u64)j * (j - 1) / 2, Q);
    row[0] = 1;                                  /* row = C(n-1, .) */
    for(long n = 1; n <= N; n++){ if(n >= 2){ for(long i = n - 1; i >= 1; i--){ row[i] += row[i - 1]; if(row[i] >= Q) row[i] -= Q; } }
      u128 acc = 0; for(long k = 1; k < n; k++) acc += (u128)((u128)row[k - 1] * b[k] % Q) * T[n - k]; b[n] = (T[n] + Q - (u64)(acc % Q)) % Q; }
    u64 D = (b[N] + 1 + Q - pw(2, (p - 1) / 2, Q)) % Q;
    if(D % p) eulerbad++;
    if(D == 0){ printf("CONNECTED-CW p=%llu\n", (unsigned long long)p); hits++; }
    np++; free(T); free(b); free(row); }
  printf("connected graphs, primes in [%llu,%llu]: %ld tested, %ld with p^2 | cg(p) - 2^((p-1)/2) + 1, rotation-congruence failures %ld\n",
         (unsigned long long)P0, (unsigned long long)P1, np, hits, eulerbad);
  return 0; }
