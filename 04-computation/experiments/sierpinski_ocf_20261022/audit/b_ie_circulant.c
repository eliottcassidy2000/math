/* Audit: exact Hamiltonian path count of a circulant tournament Cay(Z_m, D) by inclusion-exclusion
 *   H = sum_{S subset Z_m} (-1)^(m-|S|) 1^T A_S^(m-1) 1,
 * summed over rotation-orbit representatives (N(S) is rotation invariant), in unsigned __int128
 * arithmetic (exact modulo 2^128; the true H is far below 2^127 for m <= 29).
 * Usage: ./b_ie_circulant m d1 d2 ...   (connection set D)
 * Independent of the session's ie_ham code (written from scratch for the audit). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static int m;
static uint32_t outm[32];
static uint32_t FULL;
static inline uint32_t rot(uint32_t S, int t){ return ((S << t) | (S >> (m - t))) & FULL; }
static void print_u128(u128 x){ char buf[64]; int i = 63; buf[i] = 0; if(x == 0){ printf("0"); return; } while(x){ buf[--i] = '0' + (int)(x % 10); x /= 10; } printf("%s", buf + i); }
int main(int argc, char **argv){
  m = atoi(argv[1]); FULL = (m == 32) ? 0xffffffffu : ((1u << m) - 1);
  int inD[64] = {0};
  for(int a = 2; a < argc; a++) inD[atoi(argv[a]) % m] = 1;
  for(int r = 1; r < m; r++) if(inD[r] + inD[m - r] != 1){ fprintf(stderr, "not a tournament connection set (r=%d)\n", r); return 1; }
  for(int u = 0; u < m; u++){ outm[u] = 0; for(int r = 1; r < m; r++) if(inD[r]) outm[u] |= 1u << ((u + r) % m); }
  u128 total = 0; long long reps = 0;
  u128 v[32], w[32];
  for(uint64_t S64 = 1; S64 <= FULL; S64++){
    uint32_t S = (uint32_t)S64;
    /* canonical representative: S minimal among its rotations; orbit size = least period */
    int canon = 1, period = m;
    for(int t = 1; t < m; t++){ uint32_t R = rot(S, t); if(R < S){ canon = 0; break; } if(R == S && period == m) period = t; }
    if(!canon) continue;
    reps++;
    /* N(S) = number of walks with m vertices inside S */
    int k = __builtin_popcount(S);
    for(int u = 0; u < m; u++) v[u] = (S >> u & 1) ? 1 : 0;
    for(int it = 0; it < m - 1; it++){
      for(int u = 0; u < m; u++){ if(!(S >> u & 1)){ w[u] = 0; continue; } uint32_t nb = outm[u] & S; u128 s = 0; while(nb){ int x = __builtin_ctz(nb); nb &= nb - 1; s += v[x]; } w[u] = s; }
      for(int u = 0; u < m; u++) v[u] = w[u];
    }
    u128 N = 0; for(int u = 0; u < m; u++) N += v[u];
    u128 contrib = N * (u128)period;
    if((m - k) & 1) total -= contrib; else total += contrib;
  }
  printf("m=%d D=", m); for(int r = 1; r < m; r++) if(inD[r]) printf("%d,", r);
  printf(" orbit_reps=%lld H=", reps); print_u128(total); printf("\n");
  return 0;
}
