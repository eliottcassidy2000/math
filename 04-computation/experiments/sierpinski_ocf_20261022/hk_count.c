/* hk_count.c -- Held-Karp Hamiltonian path count for each tournament of a gentourng stream (n <= 16); prints H per line.
 * mac-mini-2026-10-22-sierp. */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64;
static u64 dp[1 << 16][16];
int main(void){
  char line[256];
  while (fgets(line, sizeof line, stdin)){
    if (line[0] == '>') continue;
    int L = (int)strlen(line); while (L && line[L-1] == '\n') line[--L] = 0;
    int n = 0; while (n * (n - 1) / 2 < L) n++;
    uint32_t out[16] = {0}; int p = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, p++){ if (line[p] == '1') out[i] |= 1u << j; else out[j] |= 1u << i; }
    uint32_t N = 1u << n;
    for (uint32_t m = 0; m < N; m++) for (int v = 0; v < n; v++) dp[m][v] = 0;
    for (int v = 0; v < n; v++) dp[1u << v][v] = 1;
    for (uint32_t m = 1; m < N; m++) for (int v = 0; v < n; v++){ u64 c = dp[m][v]; if (!c) continue;
      for (uint32_t q = out[v] & ~m; q; q &= q - 1){ int w = __builtin_ctz(q); dp[m | (1u << w)][w] += c; } }
    u64 H = 0; for (int v = 0; v < n; v++) H += dp[N - 1][v];
    printf("%llu\n", (unsigned long long)H);
  }
  return 0;
}
