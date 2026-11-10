/* a_classes.c (S20 audit): certify a list of tournament class representatives of order n.
   Reads codes (first token per line; bit k over pairs (i<j) in order i then j, bit set = i->j),
   computes for each the canonical code = min over all n! relabellings, and reports the number of codes,
   the number of distinct canonical codes, and whether every input is a valid tournament code.
   If #distinct == A000568(n), the list is a complete set of isomorphism-class representatives.
   usage: a_classes n file */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n, perm[10], used[10], adj[10][10], pairidx[10][10];
static u64 best;
static int cmpu(const void *a, const void *b){ u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }
static void rec(int k){
  if (k == n) {
    /* relabel: new vertex k has old vertex perm[k]; code bit for pair (a<b) = adj[perm[a]][perm[b]] */
    u64 c = 0;
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (adj[perm[a]][perm[b]]) c |= 1ULL << pairidx[a][b];
    if (c < best) best = c;
    return;
  }
  for (int x = 0; x < n; x++) if (!used[x]) { used[x] = 1; perm[k] = x; rec(k + 1); used[x] = 0; }
}
int main(int argc, char **argv){
  n = atoi(argv[1]); FILE *f = fopen(argv[2], "r"); if (!f) { perror("open"); return 1; }
  int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) pairidx[i][j] = k++;
  int npairs = k;
  u64 *canon = malloc(sizeof(u64) * 300000); long cnt = 0; char line[4096];
  while (fgets(line, sizeof line, f)) {
    u64 code; if (sscanf(line, "%llu", &code) != 1) continue;
    if (npairs < 64 && (code >> npairs)) { printf("invalid code %llu\n", code); return 1; }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) adj[i][j] = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if ((code >> pairidx[i][j]) & 1) adj[i][j] = 1; else adj[j][i] = 1; }
    best = ~0ULL; memset(used, 0, sizeof used); rec(0);
    canon[cnt++] = best;
  }
  qsort(canon, cnt, sizeof(u64), cmpu);
  long d = 0; for (long i = 0; i < cnt; i++) if (i == 0 || canon[i] != canon[i - 1]) d++;
  printf("n=%d codes=%ld distinct_iso_classes=%ld\n", n, cnt, d);
  return 0;
}
