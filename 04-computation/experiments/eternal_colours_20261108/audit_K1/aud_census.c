/* aud_census.c - audit K1 of THM-4930 (2026-11-08). Written from scratch.
 * Reads THM-4922's census file cls_n.txt (line: code H c3 strong SC aut ccode | ...; code = upper-triangle bits, bit b
 * for the b-th pair i<j in lexicographic order, set iff i -> j), and recomputes for every class, independently:
 *   H (pull subset DP), strong (Landau's score criterion), c3 (from scores), and self-converseness by a backtracking
 *   search for an anti-automorphism (a bijection p with u -> v iff p(v) -> p(u)); also |Aut| by backtracking.
 * Prints mismatches against the census fields, and a summary line per order; with argv[3] = "vals" it prints
 * "n h SCclasses strongSCclasses allclasses" for every attained h.
 * usage: aud_census n cls_n.txt [vals]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef uint32_t u32; typedef unsigned long long u64;
#define MAXN 12
static int n; static uint16_t outm[MAXN], inm[MAXN]; static u32 F[(1 << MAXN) * MAXN];
static u32 hpaths(void){
  int FN = 1 << n;
  for (int v = 0; v < n; v++) F[(1u << v) * n + v] = 1;
  for (int S = 3; S < FN; S++){ if (!(S & (S - 1))) continue;
    for (int v = 0; v < n; v++) if (S >> v & 1){ int R = S & ~(1 << v); uint16_t m = inm[v] & (uint16_t)R; u32 s = 0;
      while (m){ int u = __builtin_ctz(m); m &= m - 1; s += F[(u32)R * n + u]; } F[(u32)S * n + v] = s; } }
  u32 h = 0; for (int v = 0; v < n; v++) h += F[(u32)(FN - 1) * n + v]; return h;
}
static int landau_strong(void){
  if (n == 1) return 1;
  int sc[MAXN]; for (int v = 0; v < n; v++) sc[v] = __builtin_popcount(outm[v]);
  for (int i = 1; i < n; i++){ int x = sc[i], j = i - 1; while (j >= 0 && sc[j] > x){ sc[j + 1] = sc[j]; j--; } sc[j + 1] = x; }
  int s = 0; for (int k = 1; k < n; k++){ s += sc[k - 1]; if (s == k * (k - 1) / 2) return 0; }
  return 1;
}
static int c3(void){ int s = n * (n - 1) * (n - 2) / 6; for (int v = 0; v < n; v++){ int d = __builtin_popcount(outm[v]); s -= d * (d - 1) / 2; } return s; }
/* backtracking: anti = 1 searches anti-automorphisms (count up to 1), anti = 0 counts automorphisms */
static int p[MAXN], used[MAXN]; static long found; static int stop_at_one;
static int arc(int a, int b){ return outm[a] >> b & 1; }
static void bt(int k, int anti){
  if (stop_at_one && found) return;
  if (k == n){ found++; return; }
  int dk = __builtin_popcount(outm[k]);
  for (int w = 0; w < n; w++) if (!used[w]){
    int dw = __builtin_popcount(outm[w]);
    if (anti ? (dw != n - 1 - dk) : (dw != dk)) continue;
    int ok = 1;
    for (int j = 0; j < k && ok; j++){
      int a = arc(j, k);                       /* j -> k in T */
      int b = anti ? arc(w, p[j]) : arc(p[j], w);   /* image condition */
      if (a != b) ok = 0;
    }
    if (!ok) continue;
    p[k] = w; used[w] = 1; bt(k + 1, anti); used[w] = 0;
    if (stop_at_one && found) return;
  }
}
int main(int argc, char **argv){
  n = atoi(argv[1]); FILE *fi = fopen(argv[2], "r"); int vals = argc > 3 && !strcmp(argv[3], "vals");
  static char line[1 << 16];
  static long allc[1 << 16], scc[1 << 16], sscc[1 << 16];
  long ncls = 0, mism = 0, nsc = 0; u64 orbsum_num = 0; /* sum of n!/|Aut| must be 2^C(n,2) */
  u64 nfact = 1; for (int i = 2; i <= n; i++) nfact *= i;
  while (fgets(line, sizeof line, fi)){
    u64 code, H, cc; int c3c, strong, SC; long aut;
    if (sscanf(line, "%llu %llu %d %d %d %ld %llu", &code, &H, &c3c, &strong, &SC, &aut, &cc) != 7) continue;
    memset(outm, 0, sizeof outm); memset(inm, 0, sizeof inm);
    int b = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++){ if (code >> b & 1){ outm[i] |= 1u << j; inm[j] |= 1u << i; } else { outm[j] |= 1u << i; inm[i] |= 1u << j; } b++; }
    u32 h = hpaths(); int st = landau_strong(); int t3 = c3();
    memset(used, 0, sizeof used); found = 0; stop_at_one = 1; bt(0, 1); int sc = found > 0;
    memset(used, 0, sizeof used); found = 0; stop_at_one = 0; bt(0, 0); long au = found;
    if (h != H || st != strong || t3 != c3c || sc != SC || au != aut){ mism++; if (mism < 10) fprintf(stderr, "mismatch: %s", line); }
    ncls++; nsc += sc; orbsum_num += nfact / (u64)au;
    allc[h]++; if (sc){ scc[h]++; if (st) sscc[h]++; }
  }
  int nb = n * (n - 1) / 2;
  fprintf(stderr, "order %d: classes %ld, self-converse %ld, mismatches %ld, sum n!/|Aut| = %llu vs 2^%d = %llu\n", n, ncls, nsc, mism,
          orbsum_num, nb, 1ULL << nb);
  if (vals) for (int h = 0; h < (1 << 16); h++) if (allc[h]) printf("%d %d %ld %ld %ld\n", n, h, scc[h], sscc[h], allc[h]);
  return 0;
}
