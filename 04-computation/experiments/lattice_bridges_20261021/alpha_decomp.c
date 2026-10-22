/* alpha_decomp.c -- the OCF decomposition H = sum_k 2^k alpha_k of a circulant tournament Cay(Z_n, S), n <= 21:
 *   cyc[U] = # directed Hamiltonian cycles of T[U] (|U| odd >= 3), by the start-at-min-vertex path DP;
 *   F[W](x) = sum over sets of vertex-disjoint odd cycles with union exactly W of x^(#cycles), by the recursion on the
 *   cycle through min(W) (O(3^n)); alpha_k = sum_W [x^k] F[W]; also c_k = sum_{|U| = k} cyc[U].
 * Output: c_k, alpha_k, the contributions 2^k alpha_k and their sum (= H, checked against the class tables).
 * mac-mini-2026-10-21-fanopsp.  usage: ./alpha_decomp n s1 s2 ... s_m */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static void p128(u128 x){ char b[64]; int i=63; b[i]=0; if(!x){printf("0");return;} while(x){ b[--i]='0'+(int)(x%10); x/=10; } printf("%s",b+i); }
int main(int argc, char **argv){
  int n = atoi(argv[1]); if (n > 21) { fprintf(stderr, "n <= 21\n"); return 1; }
  int inS[64] = {0}; for (int i = 2; i < argc; i++) inS[atoi(argv[i]) % n] = 1;
  uint32_t out[32];
  for (int i = 0; i < n; i++){ out[i] = 0; for (int j = 0; j < n; j++) if (j != i && inS[((j - i) % n + n) % n]) out[i] |= 1u << j; }
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (((out[i] >> j) & 1) == ((out[j] >> i) & 1)) { fprintf(stderr, "not a tournament\n"); return 1; }
  uint32_t N = 1u << n; int K = n / 3;
  u64 *dp = calloc((size_t)N * n, sizeof(u64)), *cyc = calloc(N, sizeof(u64));
  if (!dp || !cyc) { fprintf(stderr, "alloc\n"); return 1; }
  u128 ck[32] = {0};
  for (uint32_t mask = 1; mask < N; mask++){
    int v = __builtin_ctz(mask);
    if (mask == (1u << v)) { dp[(size_t)mask * n + v] = 1; continue; }
    uint32_t rest = mask & ~(1u << v);
    for (uint32_t r = rest; r; r &= r - 1){
      int u = __builtin_ctz(r); uint32_t pm = mask & ~(1u << u); u64 s = 0;
      for (uint32_t q = pm; q; q &= q - 1){ int w = __builtin_ctz(q); if ((out[w] >> u) & 1) s += dp[(size_t)pm * n + w]; }
      dp[(size_t)mask * n + u] = s;
    }
    int sz = __builtin_popcount(mask);
    if (sz >= 3){ u64 c = 0; for (uint32_t r = rest; r; r &= r - 1){ int u = __builtin_ctz(r); if ((out[u] >> v) & 1) c += dp[(size_t)mask * n + u]; }
      ck[sz] += c; if (sz & 1) cyc[mask] = c; }
  }
  free(dp);
  u64 *F = calloc((size_t)N * (K + 1), sizeof(u64)); if (!F) { fprintf(stderr, "alloc F\n"); return 1; }
  F[0] = 1;
  u128 alpha[16] = {0}; alpha[0] = 1;
  for (uint32_t W = 1; W < N; W++){
    int v = __builtin_ctz(W); uint32_t Wp = W & ~(1u << v);
    u64 *FW = F + (size_t)W * (K + 1);
    for (uint32_t sub = Wp; ; sub = (sub - 1) & Wp){
      uint32_t U = sub | (1u << v);
      u64 c = cyc[U];
      if (c){ u64 *FR = F + (size_t)(W & ~U) * (K + 1); for (int k = 1; k <= K; k++) FW[k] += c * FR[k - 1]; }
      if (!sub) break;
    }
    for (int k = 1; k <= K; k++) alpha[k] += FW[k];
  }
  printf("n=%d S={", n); int first = 1; for (int i = 1; i < n; i++) if (inS[i]) { printf(first ? "%d" : ",%d", i); first = 0; } printf("}\n");
  printf("  c_k:"); for (int k = 3; k <= n; k++) { printf(" c%d=", k); p128(ck[k]); } printf("\n");
  u128 H = 0;
  for (int k = 0; k <= K; k++){ u128 t = alpha[k] << k; H += t; printf("  alpha_%d=", k); p128(alpha[k]); printf("  2^%d*alpha_%d=", k, k); p128(t); printf("\n"); }
  printf("  H=sum="); p128(H); printf("\n");
  return 0;
}
