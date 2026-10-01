/* Orchestrator's independent arc-HP engine for the audit of the petersen lane (2026-10-01).
 * General digraphs (not only tournaments). Written from the definition; the lane's code was not read.
 * usage: arcs N m_0 m_1 ... m_{N-1}    (m_i = out-neighbourhood bitmask of vertex i, hex)
 * prints: H <total directed Hamiltonian paths>, then one line per arc "u v c(u->v)".
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
int main(int argc, char **argv) {
  int N = atoi(argv[1]);
  if (N < 1 || N > 22 || argc < 2 + N) { fprintf(stderr, "bad args\n"); return 1; }
  uint32_t OUT[32] = {0}, IN[32] = {0};
  for (int i = 0; i < N; i++) { OUT[i] = (uint32_t)strtoul(argv[2 + i], 0, 16); }
  for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (OUT[i] >> j & 1) IN[j] |= 1u << i;
  size_t M = (size_t)1 << N;
  u64 *F = calloc(M * N, 8), *B = calloc(M * N, 8);
  if (!F || !B) { fprintf(stderr, "oom\n"); return 1; }
  for (int v = 0; v < N; v++) { F[((size_t)1 << v) * N + v] = 1; B[((size_t)1 << v) * N + v] = 1; }
  for (size_t m = 1; m < M; m++) for (int v = 0; v < N; v++) {
    if (!(m >> v & 1)) continue;
    u64 f = F[m * N + v], b = B[m * N + v];
    if (f) { uint32_t t = OUT[v] & ~(uint32_t)m; while (t) { int w = __builtin_ctz(t); t &= t - 1; F[(m | (size_t)1 << w) * N + w] += f; } }
    if (b) { uint32_t t = IN[v] & ~(uint32_t)m; while (t) { int w = __builtin_ctz(t); t &= t - 1; B[(m | (size_t)1 << w) * N + w] += b; } }
  }
  size_t full = M - 1; u64 H = 0;
  for (int v = 0; v < N; v++) H += F[full * N + v];
  printf("H %llu\n", H);
  static u64 C[32][32]; memset(C, 0, sizeof C);
  for (size_t A = 1; A < full; A++) {
    size_t R = full ^ A;
    for (int u = 0; u < N; u++) if (A >> u & 1) {
      u64 f = F[A * N + u]; if (!f) continue;
      uint32_t t = OUT[u] & (uint32_t)R;
      while (t) { int v = __builtin_ctz(t); t &= t - 1; C[u][v] += f * B[R * N + v]; }
    }
  }
  for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if (OUT[u] >> v & 1) printf("%d %d %llu\n", u, v, C[u][v]);
  return 0;
}
