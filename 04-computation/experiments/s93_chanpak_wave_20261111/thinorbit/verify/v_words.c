/* Verifier: what does the lane's "number of words with readout <= X" count?
 * (a) lane-style count: from each orbit state with h = 1, words whose first letter leaves h = 1, readout <= X,
 *     no letter acting as a fixed point (child == parent).
 * (b) words from (1,0,1) itself with readout <= X, no fixed-point letter (prefixes may wander among h = 1 states).
 * usage: v_words MATFILE X
 */
#include <stdio.h>
#include <stdlib.h>
typedef long long ll;
static ll M[4096][9]; static int nm; static ll LIM;
static unsigned long long dfs(ll h, ll e, ll o, int depth, int *maxdepth) {
  unsigned long long v = 1;
  if (depth > *maxdepth) *maxdepth = depth;
  if (depth > 2000) { fprintf(stderr, "depth > 2000\n"); exit(1); }
  for (int k = 0; k < nm; k++) {
    const ll *m = M[k];
    ll h2 = m[0] * h + m[1] * e + m[2] * o;
    if (h2 > LIM) continue;
    ll e2 = m[3] * h + m[4] * e + m[5] * o, o2 = m[6] * h + m[7] * e + m[8] * o;
    if (h2 == h && e2 == e && o2 == o) continue;
    v += dfs(h2, e2, o2, depth + 1, maxdepth);
  }
  return v;
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); LIM = atoll(argv[2]);
  while (fscanf(f, "%lld %lld %lld %lld %lld %lld %lld %lld %lld", &M[nm][0], &M[nm][1], &M[nm][2], &M[nm][3], &M[nm][4], &M[nm][5], &M[nm][6], &M[nm][7], &M[nm][8]) == 9) nm++;
  fclose(f);
  /* states with h = 1: closure */
  ll S[64][3]; int ns = 1; S[0][0] = 1; S[0][1] = 0; S[0][2] = 1;
  for (int i = 0; i < ns; i++) for (int k = 0; k < nm; k++) {
    const ll *m = M[k]; ll h = S[i][0], e = S[i][1], o = S[i][2];
    ll h2 = m[0] * h + m[1] * e + m[2] * o; if (h2 != 1) continue;
    ll e2 = m[3] * h + m[4] * e + m[5] * o, o2 = m[6] * h + m[7] * e + m[8] * o;
    int dup = 0; for (int j = 0; j < ns; j++) if (S[j][1] == e2 && S[j][2] == o2) dup = 1;
    if (!dup && ns < 64) { S[ns][0] = 1; S[ns][1] = e2; S[ns][2] = o2; ns++; }
  }
  printf("%d generators, X = %lld; orbit states with h = 1: %d:", nm, LIM, ns);
  for (int i = 0; i < ns; i++) printf(" (1,%lld,%lld)", S[i][1], S[i][2]);
  printf("\n");
  unsigned long long a = 0; int md = 0;
  for (int i = 0; i < ns; i++) for (int k = 0; k < nm; k++) {
    const ll *m = M[k]; ll h = S[i][0], e = S[i][1], o = S[i][2];
    ll h2 = m[0] * h + m[1] * e + m[2] * o; if (h2 <= 1 || h2 > LIM) continue;
    ll e2 = m[3] * h + m[4] * e + m[5] * o, o2 = m[6] * h + m[7] * e + m[8] * o;
    a += dfs(h2, e2, o2, 1, &md);
  }
  printf("(a) lane-style visits: %llu (max depth %d)\n", a, md);
  md = 0;
  unsigned long long b = dfs(1, 0, 1, 0, &md) - 1;
  printf("(b) nonempty words from (1,0,1) with readout <= X, fixed-point letters excluded: %llu (max word length %d)\n", b, md);
  return 0;
}
