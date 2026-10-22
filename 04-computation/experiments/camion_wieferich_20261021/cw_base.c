/* cw_base.c -- Camion-Wieferich primes to base a: p^2 | s_a(p) - a^((p-1)/2), where s_a = e - T_a^(-1), T_a(n) = a^C(n,2)
 * (a = 2 tournaments, a = 3 strong semicomplete digraphs).  Direct Wright recurrence mod p^2 (O(p^2) per prime).
 * Also reports whether p is base-a Wieferich.  mac-mini-2026-10-21-fanopsp.  Usage: ./cw_base a P0 P1 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 pw(u64 b, u64 e, u64 m){ u64 r = 1 % m; b %= m; while(e){ if(e & 1) r = (u128)r * b % m; b = (u128)b * b % m; e >>= 1; } return r; }
static int isprime(u64 n){ if(n < 2) return 0; for(u64 d = 2; d * d <= n; d++) if(n % d == 0) return 0; return 1; }
int main(int argc, char **argv){
  u64 a = atoll(argv[1]), P0 = atoll(argv[2]), P1 = atoll(argv[3]); long hits = 0, np = 0;
  for(u64 p = (P0 < 3 ? 3 : P0); p <= P1; p++){ if(!isprime(p) || a % p == 0) continue;
    u64 Q = p * p; long N = (long)p;
    u64 *T = malloc((N + 1) * 8), *s = calloc(N + 1, 8), *row = calloc(N + 2, 8);
    for(long j = 0; j <= N; j++) T[j] = pw(a, (u64)j * (j - 1) / 2, Q);
    row[0] = 1;
    for(long n = 1; n <= N; n++){ for(long i = n - 1; i >= 1; i--){ row[i] += row[i - 1]; if(row[i] >= Q) row[i] -= Q; } row[n] = 1;
      u128 acc = 0; for(long k = 1; k < n; k++) acc += (u128)((u128)row[k] * s[k] % Q) * T[n - k]; s[n] = (T[n] + Q - (u64)(acc % Q)) % Q; }
    u64 D = (s[N] + Q - pw(a, (p - 1) / 2, Q)) % Q; int wief = (pw(a, p - 1, Q) == 1);
    if(D % p) printf("  p=%llu: Euler criterion FAILS (bug)\n", (unsigned long long)p);
    if(D == 0){ printf("CAMION-WIEFERICH base %llu: p=%llu%s\n", (unsigned long long)a, (unsigned long long)p, wief ? " (also base-a Wieferich)" : ""); hits++; }
    else if(wief) printf("  base-%llu Wieferich prime p=%llu is not Camion-Wieferich\n", (unsigned long long)a, (unsigned long long)p);
    np++; free(T); free(s); free(row); }
  printf("base %llu, primes in [%llu,%llu]: %ld tested, %ld Camion-Wieferich\n", (unsigned long long)a, (unsigned long long)P0, (unsigned long long)P1, np, hits);
  return 0; }
