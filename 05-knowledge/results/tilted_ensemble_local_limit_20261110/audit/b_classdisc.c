// b_classdisc.c -- AUDIT: exact LABELLED distribution of H mod 8 at order n from the isomorphism classes (gentourng),
// each class weighted by n!/|Aut|.  Independent cross-check of the exact mod-4 discrepancies that b_thmC.c derives at
// order m+1 from the vertex cubes of the order-m classes (the identity in the proof of Theorem C (4)).
// usage: gentourng -q n [res/mod] | ./b_classdisc n      (prints raw weights; parts can be added)
#include <stdio.h>
#include <stdlib.h>
typedef long long i64; typedef unsigned long long u64;
#define MAXN 10
static int n, full, adj[MAXN][MAXN], inm[MAXN], sc[MAXN], perm_[MAXN], used_;
static i64 E[1 << MAXN][MAXN];
static i64 autcount(int k) {
  if (k == n) return 1; i64 c = 0;
  for (int w = 0; w < n; w++) { if (used_ >> w & 1) continue; if (sc[w] != sc[k]) continue; int ok = 1;
    for (int j = 0; j < k && ok; j++) if (adj[j][k] != adj[perm_[j]][w]) ok = 0;
    if (ok) { perm_[k] = w; used_ |= 1 << w; c += autcount(k + 1); used_ &= ~(1 << w); } }
  return c;
}
int main(int argc, char **argv) {
  n = atoi(argv[1]); full = (1 << n) - 1; i64 nf = 1; for (int i = 2; i <= n; i++) nf *= i;
  char line[256]; u64 nt = 0; i64 w8[8] = {0}, tot = 0;
  while (fgets(line, sizeof line, stdin)) {
    if (line[0] != '0' && line[0] != '1') continue; int k = 0;
    for (int i = 0; i < n; i++) { inm[i] = 0; sc[i] = 0; adj[i][i] = 0; }
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { int b = line[k++] == '1'; adj[i][j] = b; adj[j][i] = !b; if (b) { inm[j] |= 1 << i; sc[i]++; } else { inm[i] |= 1 << j; sc[j]++; } }
    for (int S = 1; S <= full; S++) { int vv = S; while (vv) { int v = __builtin_ctz(vv); vv &= vv - 1; int R = S ^ (1 << v);
        if (!R) { E[S][v] = 1; continue; } i64 e = 0; int pp = R & inm[v]; while (pp) { int u = __builtin_ctz(pp); pp &= pp - 1; e += E[R][u]; } E[S][v] = e; } }
    i64 H = 0; for (int v = 0; v < n; v++) H += E[full][v];
    used_ = 0; i64 au = autcount(0); if (nf % au) { fprintf(stderr, "bad aut\n"); return 1; }
    w8[H & 7] += nf / au; tot += nf / au; nt++;
  }
  printf("order %d: classes %llu, total labelled weight %lld; weight by H mod 8 (1,3,5,7): %lld %lld %lld %lld; even residues: %lld; W(1 mod 4) - W(3 mod 4) = %lld\n",
         n, nt, tot, w8[1], w8[3], w8[5], w8[7], w8[0] + w8[2] + w8[4] + w8[6], w8[1] + w8[5] - w8[3] - w8[7]);
  return 0;
}
