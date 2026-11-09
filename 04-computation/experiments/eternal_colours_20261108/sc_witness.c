/* sc_witness.c - mac-mini-2026-11-08-firstrealizers. Constructive self-converse witnesses for odd Hamiltonian-path counts.
 * The transitive tournament TT_n (i -> j for i < j) is self-converse under i -> n-1-i. Reversing a set of arcs that is
 * closed under the mirror (i,j) -> (n-1-j, n-1-i) keeps it self-converse. We enumerate every mirror-closed set made of at
 * most K mirror-orbits of arcs, for n = NMIN..NMAX, compute H by subset DP, and print, for each H <= HMAX that occurs,
 * the least order n, a strongness flag, and the reversed arcs of one witness.
 * usage: sc_witness NMIN NMAX K HMAX > sc_witness.out */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
#define MAXN 14
static int n;
static uint16_t adj[MAXN];
static u64 f[1 << MAXN][MAXN];
static u64 hpaths(void){
  int F = 1 << n;
  for (int S = 0; S < F; S++) for (int v = 0; v < n; v++) f[S][v] = 0;
  for (int v = 0; v < n; v++) f[1 << v][v] = 1;
  for (int S = 1; S < F; S++) for (int v = 0; v < n; v++) if (f[S][v]) {
    u64 x = f[S][v]; uint16_t out = adj[v] & ~S;
    while (out) { int u = __builtin_ctz(out); out &= out - 1; f[S | 1 << u][u] += x; } }
  u64 h = 0; for (int v = 0; v < n; v++) h += f[F - 1][v]; return h;
}
static int strongq(void){
  uint16_t all = (uint16_t)((1u << n) - 1), seen = 1, fr = 1;
  while (fr) { uint16_t nx = 0; for (int v = 0; v < n; v++) if (fr >> v & 1) nx |= adj[v]; nx &= ~seen; seen |= nx; fr = nx; }
  if (seen != all) return 0;
  seen = 1; fr = 1;
  while (fr) { uint16_t nx = 0; for (int v = 0; v < n; v++) if (!(seen >> v & 1)) for (int u = 0; u < n; u++) if ((fr >> u & 1) && (adj[v] >> u & 1)) { nx |= 1u << v; break; } seen |= nx; fr = nx; }
  return seen == all;
}
int main(int argc, char **argv){
  int NMIN = atoi(argv[1]), NMAX = atoi(argv[2]), K = atoi(argv[3]); long HMAX = atol(argv[4]);
  int *best_n = calloc(HMAX + 1, sizeof(int)); int *best_s = calloc(HMAX + 1, sizeof(int));
  char (*wit)[160] = calloc(HMAX + 1, 160);
  for (n = NMIN; n <= NMAX; n++) {
    int oi[200], oj[200], no = 0;                     /* mirror orbits of arcs (i<j) */
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
      int mi = n - 1 - j, mj = n - 1 - i;             /* mirror of (i,j) */
      if (mi < i || (mi == i && mj < j)) continue;    /* keep one representative per orbit */
      oi[no] = i; oj[no] = j; no++; }
    long cnt = 0;
    int sel[8];
    for (int k = 0; k <= K; k++) {
      for (int t = 0; t < k; t++) sel[t] = t;
      while (1) {
        if (k > no) break;
        for (int v = 0; v < n; v++) adj[v] = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) adj[i] |= 1u << j;
        for (int t = 0; t < k; t++) {
          int i = oi[sel[t]], j = oj[sel[t]], mi = n - 1 - j, mj = n - 1 - i;
          adj[i] &= ~(1u << j); adj[j] |= 1u << i;
          if (!(mi == i && mj == j)) { adj[mi] &= ~(1u << mj); adj[mj] |= 1u << mi; } }
        u64 h = hpaths(); cnt++;
        if (h <= (u64)HMAX && (best_n[h] == 0 || (best_n[h] == n && !best_s[h] && strongq()))) {
          /* first witness of the least order, or a strong one replacing a non-strong one of the same order: the printed
             arcs always describe a tournament whose strongness is the printed flag (audit K1 correction 18) */
          best_n[h] = n; best_s[h] = strongq(); int p = 0; p += snprintf(wit[h] + p, 160 - p, "rev");
          for (int t = 0; t < k && p < 150; t++) p += snprintf(wit[h] + p, 160 - p, " %d%c%d", oi[sel[t]], '-', oj[sel[t]]);
        }
        if (k == 0) break;
        int t = k - 1; while (t >= 0 && sel[t] == no - k + t) t--;
        if (t < 0) break;
        sel[t]++; for (int u = t + 1; u < k; u++) sel[u] = sel[u - 1] + 1;
      }
    }
    fprintf(stderr, "n=%d: %d mirror orbits, %ld self-converse tournaments enumerated\n", n, no, cnt);
  }
  for (long h = 1; h <= HMAX; h += 2) if (best_n[h]) printf("W %ld %d %d %s\n", h, best_n[h], best_s[h], wit[h]);
  return 0;
}
