/* alpha_census.c - opus-2026-11-08-S86.
 * Reads the class files written by kind-pasteur's classgen.c (05-knowledge/results/first_creators_genealogy_20261107)
 * line format: code H c3 strong SC aut ccode | del_code:del_H ...
 * (code = canonical upper-triangle bit code: bit b for the b-th pair (i,j), i<j in row order, set iff i -> j).
 * For every class it computes the odd-cycle packing vector alpha_k (k = 0..4) of OCF, H = sum alpha_k 2^k (checked
 * against the H column), whether alpha_2 = 0 ("intersecting": no two vertex-disjoint odd cycles), whether some vertex
 * meets every cycle (T - a acyclic, i.e. has no cyclic triangle), and the strong-component profile.
 * usage: alpha_census n < cls_n.txt > alpha_n.txt
 * output line: code H strong SC a1 a2 a3 a4 inter apexvertex ncomp compsizes(sorted, '-' separated)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MAXN 10
static int n;
static unsigned out[MAXN];
static long long cyc[1 << MAXN];      /* directed odd cycles with vertex set exactly S */
static long long dp[1 << MAXN][MAXN]; /* paths from min(S) to v covering S */
static long long F[1 << MAXN][6];     /* packings covering exactly U, by number of cycles */

static void decode(u64 code){
  memset(out, 0, sizeof out); int b = 0;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++){ if (code >> b & 1) out[i] |= 1u << j; else out[j] |= 1u << i; b++; }
}
static int popc(unsigned x){ return __builtin_popcount(x); }
static void cycles(void){
  int N = 1 << n;
  for (int S = 0; S < N; S++){ cyc[S] = 0; for (int v = 0; v < n; v++) dp[S][v] = 0; }
  for (int s = 0; s < n; s++) dp[1u << s][s] = 1;
  for (int S = 1; S < N; S++){
    int s = __builtin_ctz(S);
    for (int v = 0; v < n; v++){ long long c = dp[S][v]; if (!c) continue;
      unsigned nb = out[v] & ~(unsigned)S & ~((1u << (s + 1)) - 1);   /* only vertices > s */
      while (nb){ int u = __builtin_ctz(nb); nb &= nb - 1; dp[S | (1u << u)][u] += c; } }
    int k = popc(S);
    if (k >= 3 && (k & 1)){ long long t = 0; for (int v = 0; v < n; v++) if (dp[S][v] && (out[v] >> s & 1)) t += dp[S][v]; cyc[S] = t; }
  }
}
static void packings(long long *alpha){
  int N = 1 << n;
  for (int U = 0; U < N; U++) for (int k = 0; k < 6; k++) F[U][k] = 0;
  F[0][0] = 1;
  for (int U = 1; U < N; U++){
    int m = __builtin_ctz(U); unsigned rest = U & ~(1u << m);
    /* S = {m} u T, T subset of rest */
    for (unsigned T = rest;; T = (T - 1) & rest){
      unsigned S = T | (1u << m);
      if (cyc[S]) for (int k = 0; k < 5; k++) if (F[U & ~S][k]) F[U][k + 1] += cyc[S] * F[U & ~S][k];
      if (!T) break;
    }
  }
  for (int k = 0; k < 6; k++){ alpha[k] = 0; for (int U = 0; U < N; U++) alpha[k] += F[U][k]; }
}
static int has_cyclic_triangle(unsigned mask){
  for (int a = 0; a < n; a++) if (mask >> a & 1) for (int b = 0; b < n; b++) if ((mask >> b & 1) && (out[a] >> b & 1))
    for (int c = 0; c < n; c++) if ((mask >> c & 1) && (out[b] >> c & 1) && (out[c] >> a & 1)) return 1;
  return 0;
}
/* strong components via reachability closure */
static int comps(int *sizes){
  unsigned reach[MAXN];
  for (int v = 0; v < n; v++){ unsigned r = 1u << v, f = r; while (f){ unsigned g = 0; for (int u = 0; u < n; u++) if (f >> u & 1) g |= out[u]; g &= ~r; r |= g; f = g; } reach[v] = r; }
  int seen = 0, c = 0;
  for (int v = 0; v < n; v++){ if (seen >> v & 1) continue; unsigned cl = 0; for (int u = 0; u < n; u++) if ((reach[v] >> u & 1) && (reach[u] >> v & 1)) cl |= 1u << u; seen |= cl; sizes[c++] = popc(cl); }
  for (int i = 1; i < c; i++){ int x = sizes[i], j = i - 1; while (j >= 0 && sizes[j] < x){ sizes[j + 1] = sizes[j]; j--; } sizes[j + 1] = x; }
  return c;
}
int main(int argc, char **argv){
  if (argc < 2){ fprintf(stderr, "usage: alpha_census n < cls_n.txt\n"); return 1; }
  n = atoi(argv[1]);
  char line[1 << 16]; long rows = 0, bad = 0;
  while (fgets(line, sizeof line, stdin)){
    u64 code, H, cc; int c3, strong, sc; long aut;
    if (sscanf(line, "%llu %llu %d %d %d %ld %llu", &code, &H, &c3, &strong, &sc, &aut, &cc) != 7) continue;
    decode(code); cycles(); long long al[6]; packings(al);
    long long Hc = 0; for (int k = 0; k < 6; k++) Hc += al[k] << k;
    if ((u64)Hc != H) bad++;
    int apex = 0; unsigned full = (1u << n) - 1;
    if (al[1] == 0) apex = 1; else for (int a = 0; a < n; a++) if (!has_cyclic_triangle(full & ~(1u << a))){ apex = 1; break; }
    int sz[MAXN]; int nc = comps(sz);
    printf("%llu %llu %d %d %lld %lld %lld %lld %d %d %d ", code, H, strong, sc, al[1], al[2], al[3], al[4], al[2] == 0, apex, nc);
    for (int i = 0; i < nc; i++) printf(i ? "-%d" : "%d", sz[i]);
    printf("\n"); rows++;
  }
  fprintf(stderr, "order %d: %ld classes, %ld H mismatches\n", n, rows, bad);
  return bad ? 2 : 0;
}
