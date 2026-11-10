// x5_altrank.c -- AUDIT 2 (auditor's own code).  Brute force over ALL alternating m x m matrices over F_2 (all graphs
// on m vertices): number of matrices by (rank, parity of the number of edges, [all degrees even]).
// usage: x5_altrank m [r k]      (graphs with code = r mod k; parts can be added)       m <= 8
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64;
int main(int argc, char **argv) { int m = atoi(argv[1]); int r = argc > 3 ? atoi(argv[2]) : 0, k = argc > 3 ? atoi(argv[3]) : 1;
  int np = m * (m - 1) / 2; static u64 cnt[12][2][2]; int pi[40], pj[40], b = 0;
  for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { pi[b] = i; pj[b] = j; b++; }
  for (u64 code = (u64)r; code < (1ULL << np); code += (u64)k) { unsigned rows[10] = {0};
    for (int t = 0; t < np; t++) if (code >> t & 1) { rows[pi[t]] |= 1u << pj[t]; rows[pj[t]] |= 1u << pi[t]; }
    int ev = 1; for (int i = 0; i < m; i++) if (__builtin_popcount(rows[i]) & 1) ev = 0;
    int rk = 0; for (int c = 0; c < m; c++) { int p = -1; for (int i = rk; i < m; i++) if (rows[i] >> c & 1) { p = i; break; } if (p < 0) continue;
      unsigned t = rows[p]; rows[p] = rows[rk]; rows[rk] = t; for (int i = rk + 1; i < m; i++) if (rows[i] >> c & 1) rows[i] ^= rows[rk]; rk++; }
    cnt[rk][__builtin_popcountll(code) & 1][ev]++; }
  printf("m %d", m); for (int rk = 0; rk <= m; rk++) for (int e = 0; e < 2; e++) for (int f = 0; f < 2; f++) if (cnt[rk][e][f]) printf(" r%d_e%d_f%d %llu", rk, e, f, cnt[rk][e][f]);
  printf("\n"); return 0; }
