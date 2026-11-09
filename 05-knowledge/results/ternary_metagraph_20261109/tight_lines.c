/* tight_lines.c - opus-2026-11-09-S87.
 * Tight lines of the tiling line metagraph (THM-643 frame): lines {t, t + 1} with nu(t) + nu(t + 1) = floor(n/3),
 * the equality case of the bound nu(t) + nu(t + 1) >= floor(n/3) (Session26's consecutive-triples argument;
 * THM-4951).  nu = maximum number of
 * vertex-disjoint cycles = maximum number of vertex-disjoint 3-cycles.
 * Tiling t: bit (i, j), j >= i + 2, set iff the arc is reversed (j -> i); path arcs i -> i + 1 are fixed.
 * Search: assign the tiles (i, v), i <= v - 2, vertex by vertex.  nu of the prefix tournament on {0..v} never exceeds
 * nu of the whole tournament, so a prefix with nu_t + nu_u > floor(n/3) is pruned (u = complement).  Each line is
 * counted once: its representative t has tile (0, 2) reversed.
 * Output: counts by unordered nu-pair and colour (blue = sigma-fixed, sigma(i, j) = (n-1-j, n-1-i)), and the
 * histogram of the triple word x_a = bit (a, a + 2), a = 0..n-3, of the representative.
 * usage: tight_lines n SLICE NSLICES SPLIT_V   (prefixes at vertex SPLIT_V are dealt round-robin to slices)
 *        tight_lines n 0 1 -1 WORD            (first-hit mode: triple word fixed to WORD, a string of n - 2 bits with
 *                                              WORD[0] = 1; prints the first tight line found and stops)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXN 16
static int n, K, slice, ns, splitv;
static uint32_t ot[MAXN], ou[MAXN];        /* out-neighbourhoods of t and of u = t + 1 (prefix) */
static int bit[MAXN][MAXN];                /* bit[i][j] for j >= i + 2 */
static long long cnt[8][8][2], nodes, splitcount;
static long long wordhist[1 << 14][2];
static uint32_t tri[MAXN][128]; static int ntri[MAXN];
static int best_found;
static int forced, fw[MAXN];
static void packdfs(uint32_t avail, int cur, int cap){
  if (cur > best_found) best_found = cur;
  if (best_found >= cap) return;
  if (cur + __builtin_popcount(avail) / 3 <= best_found) return;
  if (!avail) return;
  int a = __builtin_ctz(avail);
  for (int k = 0; k < ntri[a]; k++) if ((tri[a][k] & avail) == tri[a][k]) packdfs(avail & ~tri[a][k], cur + 1, cap);
  packdfs(avail & ~(1u << a), cur, cap);
}
static int nu(const uint32_t *out, int V, int cap){   /* max disjoint 3-cycles on {0..V-1}, capped at cap */
  uint32_t in[MAXN]; memset(in, 0, sizeof in);
  for (int a = 0; a < V; a++) for (int b = 0; b < V; b++) if (out[a] >> b & 1) in[b] |= 1u << a;
  for (int a = 0; a < V; a++){ ntri[a] = 0; uint32_t above = (V >= 32 ? 0xffffffffu : ((1u << V) - 1)) & ~((2u << a) - 1);
    uint32_t bs = out[a] & above;
    while (bs){ int b = __builtin_ctz(bs); bs &= bs - 1; uint32_t cs = out[b] & in[a] & above;
      while (cs){ int c = __builtin_ctz(cs); cs &= cs - 1; tri[a][ntri[a]++] = (1u << a) | (1u << b) | (1u << c); } } }
  best_found = 0; packdfs((1u << V) - 1, 0, cap); return best_found;
}
static void record(void){
  int a = nu(ot, n, K + 1), b = nu(ou, n, K + 1);
  if (a + b != K) return;                  /* cannot be < K by the bound; > K excluded */
  int blue = 1;
  for (int i = 0; i < n && blue; i++) for (int j = i + 2; j < n; j++){ int si = n - 1 - j, sj = n - 1 - i; if (bit[i][j] != bit[si][sj]){ blue = 0; break; } }
  int lo = a < b ? a : b, hi = a < b ? b : a;
  cnt[lo][hi][blue]++;
  int w = 0; for (int x = 0; x + 2 < n; x++) w |= bit[x][x + 2] << x;
  if (n - 2 <= 14) wordhist[w][blue]++;
  if (forced){
    printf("FOUND n=%d word ", n); for (int x = 0; x + 2 < n; x++) putchar('0' + bit[x][x + 2]);
    printf(" nu %d %d %s tiles(i,j reversed):", a, b, blue ? "blue" : "black");
    for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) if (bit[i][j]) printf(" %d,%d", i, j);
    printf("\n"); fflush(stdout); exit(0);
  }
}
static void dfs(int v){
  if (v == n){ record(); return; }
  if (v == splitv){ if ((splitcount++ % ns) != slice) return; }
  nodes++;
  int nt = v - 1;                          /* tiles (i, v), i = 0..v-2 */
  for (uint32_t mask = 0; mask < (1u << nt); mask++){
    if (v == 2 && mask != 1) continue;     /* representative: tile (0, 2) reversed */
    if (forced && (int)(mask >> (v - 2) & 1) != fw[v - 2]) continue;   /* triple tile (v - 2, v) fixed */
    ot[v] = 0; ou[v] = 0;
    for (int i = 0; i < v; i++){ ot[i] &= ~(1u << v); ou[i] &= ~(1u << v); }
    ot[v - 1] |= 1u << v; ou[v - 1] |= 1u << v;   /* path arc */
    for (int i = 0; i < nt; i++){ int r = mask >> i & 1; bit[i][v] = r;
      if (r){ ot[v] |= 1u << i; ou[i] |= 1u << v; } else { ot[i] |= 1u << v; ou[v] |= 1u << i; } }
    int a = nu(ot, v + 1, K + 1); if (a > K) continue;
    int b = nu(ou, v + 1, K + 1 - a); if (a + b > K) continue;
    dfs(v + 1);
  }
  ot[v] = ou[v] = 0;
  for (int i = 0; i < v; i++){ ot[i] &= ~(1u << v); ou[i] &= ~(1u << v); }
}
int main(int argc, char **argv){
  n = atoi(argv[1]); slice = argc > 2 ? atoi(argv[2]) : 0; ns = argc > 3 ? atoi(argv[3]) : 1; splitv = argc > 4 ? atoi(argv[4]) : -1;
  K = n / 3;
  if (argc > 5){ forced = 1; for (int x = 0; x + 2 < n; x++) fw[x] = argv[5][x] - '0'; }
  memset(ot, 0, sizeof ot); memset(ou, 0, sizeof ou);
  if (n >= 2){ ot[0] |= 2; ou[0] |= 2; }
  dfs(2);
  long long tot = 0, totb = 0;
  for (int i = 0; i < 8; i++) for (int j = i; j < 8; j++) for (int c = 0; c < 2; c++) if (cnt[i][j][c]){
    printf("pair %d %d %s %lld\n", i, j, c ? "blue" : "black", cnt[i][j][c]); tot += cnt[i][j][c]; if (c) totb += cnt[i][j][c]; }
  if (n - 2 <= 14) for (int w = 0; w < (1 << (n - 2)); w++) for (int c = 0; c < 2; c++) if (wordhist[w][c]){
    printf("word "); for (int x = 0; x + 2 < n; x++) putchar('0' + (w >> x & 1)); printf(" %s %lld\n", c ? "blue" : "black", wordhist[w][c]); }
  fprintf(stderr, "n=%d slice %d/%d: tight lines %lld (blue %lld), search nodes %lld\n", n, slice, ns, tot, totb, nodes);
  return 0;
}
