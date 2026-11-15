// Coverage of the marked two-block one-apex family: N = a_B (n_L - g) + g n_B + s_L t_B + d_L f_B,
// L ranges over right-profiles (g,n_L,s_L,d_L) of marked tournaments of order <= ML, B over left-profiles
// (a_B,n_B,t_B,f_B) of order <= MB. Reports the coverage of odd N up to X.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
typedef struct { i64 x0, x1, x2, x3; } p4;
static int cmpB(const void *A, const void *B) { const p4 *a = A, *b = B; return (a->x1 > b->x1) - (a->x1 < b->x1); }
static p4 *load(int maxm, int right, int *cnt) {
  size_t cap = 4000000; p4 *v = malloc(sizeof(p4) * cap); int c = 0;
  for (int m = 1; m <= maxm; m++) {
    char fn[64]; sprintf(fn, "pal_%d.txt", m); FILE *fp = fopen(fn, "r"); char line[512];
    while (fgets(line, sizeof line, fp)) { i64 a, q, s, e, t, f, n; sscanf(line, "%lld %lld %lld %lld %lld %lld %lld", &a, &q, &s, &e, &t, &f, &n);
      p4 p; if (right) { p.x0 = a; p.x1 = n; p.x2 = s; p.x3 = e; } else { p.x0 = a; p.x1 = n; p.x2 = t; p.x3 = f; }
      v[c++] = p; }
    fclose(fp);
  }
  // dedupe
  qsort(v, c, sizeof(p4), cmpB);
  int k = 0; for (int i = 0; i < c; i++) { int dup = 0; for (int j = k - 1; j >= 0 && v[j].x1 == v[i].x1; j--) if (!memcmp(&v[j], &v[i], sizeof(p4))) { dup = 1; break; } if (!dup) v[k++] = v[i]; }
  *cnt = k; return v;
}
int main(int argc, char **argv) {
  int ML = atoi(argv[1]), MB = atoi(argv[2]); i64 X = atoll(argv[3]);
  int nL, nB; p4 *L = load(ML, 1, &nL), *B = load(MB, 0, &nB);
  fprintf(stderr, "L profiles %d, B profiles %d\n", nL, nB);
  unsigned char *hit = calloc(X / 16 + 2, 1);   // odd numbers only: bit (N-1)/2
  long long pairs = 0;
  for (int i = 0; i < nL; i++) {
    i64 g = L[i].x0, nl = L[i].x1, sl = L[i].x2, dl = L[i].x3;
    for (int j = 0; j < nB; j++) {
      i64 nb = B[j].x1; if (g * nb > X) break;     // B sorted by n_B, all terms nonnegative
      i64 N = B[j].x0 * (nl - g) + g * nb + sl * B[j].x2 + dl * B[j].x3;
      if (N <= X) { i64 k = (N - 1) / 2; hit[k >> 3] |= 1 << (k & 7); }
      pairs++;
    }
  }
  i64 miss = 0, firstbig = -1; i64 lastmiss = 0; int shown = 0;
  printf("misses: ");
  for (i64 N = 1; N <= X; N += 2) { i64 k = (N - 1) / 2; if (!(hit[k >> 3] & (1 << (k & 7)))) { miss++; lastmiss = N; if (shown < 40) { printf("%lld ", N); shown++; } } }
  printf("\nML=%d MB=%d X=%lld pairs=%lld odd misses=%lld largest miss=%lld\n", ML, MB, X, pairs, miss, lastmiss);
  // first gap structure: report the longest initial run of coverage beyond 21
  i64 run_end = 0; for (i64 N = 23; N <= X; N += 2) { i64 k = (N - 1) / 2; if (!(hit[k >> 3] & (1 << (k & 7)))) { run_end = N; break; } }
  printf("all odd N in [23, %lld) covered\n", run_end ? run_end : X);
  return 0;
}
