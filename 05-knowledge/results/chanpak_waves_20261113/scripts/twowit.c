// Witnesses for sampled targets in the marked two-block family; prints L code/mask and B code/mask.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
typedef struct { i64 x0, x1, x2, x3; int m; char code[64]; int V; } p4;
static p4 *load(int maxm, int right, int *cnt) {
  p4 *v = malloc(sizeof(p4) * 2000000); int c = 0;
  for (int m = 1; m <= maxm; m++) { char fn[64]; sprintf(fn, "pal_%d.txt", m); FILE *fp = fopen(fn, "r"); char line[512];
    while (fgets(line, sizeof line, fp)) { i64 a, q, s, e, t, f, n; char code[64] = ""; int V;
      int k = sscanf(line, "%lld %lld %lld %lld %lld %lld %lld %63s %d", &a, &q, &s, &e, &t, &f, &n, code, &V);
      if (k < 9) { V = atoi(code); strcpy(code, "-"); }   // order 1 has an empty code field
      p4 p; p.x0 = a; p.x1 = n; p.x2 = right ? s : t; p.x3 = right ? e : f; p.m = m; strcpy(p.code, code); p.V = V; v[c++] = p; }
    fclose(fp); }
  *cnt = c; return v; }
int main(int argc, char **argv) {
  int ML = atoi(argv[1]), MB = atoi(argv[2]); int nt = argc - 3; i64 *tg = malloc(sizeof(i64) * nt); int *done = calloc(nt, sizeof(int));
  for (int i = 0; i < nt; i++) tg[i] = atoll(argv[3 + i]);
  int nL, nB; p4 *L = load(ML, 1, &nL), *B = load(MB, 0, &nB);
  int left = nt;
  for (int i = 0; i < nL && left; i++) for (int j = 0; j < nB && left; j++) {
    i64 N = B[j].x0 * (L[i].x1 - L[i].x0) + L[i].x0 * B[j].x1 + L[i].x2 * B[j].x2 + L[i].x3 * B[j].x3;
    for (int k = 0; k < nt; k++) if (!done[k] && tg[k] == N) { done[k] = 1; left--;
      printf("%lld %d %s %d %d %s %d\n", N, L[i].m, L[i].code, L[i].V, B[j].m, B[j].code, B[j].V); } }
  for (int k = 0; k < nt; k++) if (!done[k]) printf("%lld NONE\n", tg[k]);
  return 0; }
