/* THM-5071 checks.
 * (1) Weighted cycle formula: for a loopless digraph D,
 *       ham(D) = sum over permutations sigma whose nontrivial cycles are each a directed cycle of D or of
 *                its complement Dbar, of prod over cycles: +1 (cycle of D) or (-1)^(len-1) (cycle of Dbar).
 *     Checked exhaustively on every loopless digraph of order <= 4, every semicomplete digraph and every
 *     oriented graph of order <= 5, and on random digraphs of orders 6 and 7.
 * (2) One-defect spectra: H of a tournament with one arc deleted (one 00 pair) and with one reverse arc
 *     added (one 11 pair), orders 3..7.
 * Usage: digraph_weights SAMPLES SEED
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n;
static int A[8][8];
static long long dp[1 << 8][8];

static long long ham(void) {
  int N = 1 << n;
  for (int s = 0; s < N; s++) for (int v = 0; v < n; v++) dp[s][v] = 0;
  for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
  for (int s = 1; s < N; s++) for (int v = 0; v < n; v++) if (dp[s][v])
    for (int w = 0; w < n; w++) if (!(s >> w & 1) && A[v][w]) dp[s | 1 << w][w] += dp[s][v];
  long long t = 0;
  for (int v = 0; v < n; v++) t += dp[N - 1][v];
  return t;
}

static int perm[8], used[8];
static long long gs_sum;
static void gs_rec(int i) {
  if (i == n) {
    int seen[8] = {0};
    long long w = 1;
    for (int s = 0; s < n && w; s++) if (!seen[s]) {
      int len = 0, inD = 1, inDbar = 1, v = s;
      do { seen[v] = 1; int u = perm[v]; len++;
           if (u != v) { if (!A[v][u]) inD = 0; if (A[v][u]) inDbar = 0; }
           v = u; } while (v != s);
      if (len == 1) continue;
      long long c = (inD ? 1 : 0) + (inDbar ? ((len - 1) % 2 ? -1 : 1) : 0);
      w *= c;
    }
    gs_sum += w;
    return;
  }
  for (int v = 0; v < n; v++) if (!used[v]) { used[v] = 1; perm[i] = v; gs_rec(i + 1); used[v] = 0; }
}
static long long gs(void) { gs_sum = 0; memset(used, 0, sizeof used); gs_rec(0); return gs_sum; }

static unsigned long long rng = 88172645463325252ULL;
static unsigned long long xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

static long long fails, checks;
static void check(void) { long long a = ham(), b = gs(); checks++; if (a != b) { fails++; if (fails < 5) {
  printf("FAIL n=%d ham=%lld gs=%lld\n", n, a, b); } } }

int main(int argc, char **argv) {
  long long samples = argc > 1 ? atoll(argv[1]) : 100000;
  if (argc > 2) rng = strtoull(argv[2], 0, 10) | 1;
  /* (1a) all loopless digraphs, n <= 4 */
  for (n = 1; n <= 4; n++) {
    int m = n * (n - 1) / 2;
    long long total = 1; for (int i = 0; i < m; i++) total *= 4;
    for (long long code = 0; code < total; code++) {
      long long c = code; int k = 0;
      memset(A, 0, sizeof A);
      for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, k++) {
        int st = c % 4; c /= 4;
        A[i][j] = st & 1; A[j][i] = st >> 1 & 1;
      }
      check();
    }
    printf("all loopless digraphs n=%d: checks so far %lld fails %lld\n", n, checks, fails);
  }
  /* (1b) semicomplete (states 1,2,3) and oriented (states 0,1,2) n = 5 */
  for (int kind = 0; kind < 2; kind++) {
    n = 5; int m = 10; long long total = 1; for (int i = 0; i < m; i++) total *= 3;
    for (long long code = 0; code < total; code++) {
      long long c = code; memset(A, 0, sizeof A);
      for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        int st = c % 3; c /= 3;
        if (kind == 0) st += 1;           /* 1: i->j, 2: j->i, 3: both */
        A[i][j] = st & 1; A[j][i] = st >> 1 & 1;
      }
      check();
    }
    printf("%s n=5: checks so far %lld fails %lld\n", kind ? "oriented" : "semicomplete", checks, fails);
  }
  /* (1c) random general digraphs, n = 6, 7 */
  for (n = 6; n <= 7; n++) {
    long long reps = n == 6 ? samples : samples / 10;
    for (long long r = 0; r < reps; r++) {
      memset(A, 0, sizeof A);
      for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j) A[i][j] = xr() & 1;
      check();
    }
    printf("random digraphs n=%d: checks so far %lld fails %lld\n", n, checks, fails);
  }
  printf("GS weighted formula: %lld checks, %lld failures\n", checks, fails);
  /* (2) one-defect spectra */
  for (n = 3; n <= 7; n++) {
    int m = n * (n - 1) / 2;
    static char del[6000], add[6000];
    memset(del, 0, sizeof del); memset(add, 0, sizeof add);
    for (long long code = 0; code < (1LL << m); code++) {
      int pi[32], pj[32], k = 0;
      memset(A, 0, sizeof A);
      for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, k++) {
        pi[k] = i; pj[k] = j;
        if (code >> k & 1) A[i][j] = 1; else A[j][i] = 1;
      }
      for (k = 0; k < m; k++) {
        int i = pi[k], j = pj[k];
        int a = A[i][j], b = A[j][i];
        A[i][j] = 0; A[j][i] = 0;
        long long h = ham(); if (h < 6000) del[h] = 1;
        A[i][j] = 1; A[j][i] = 1;
        h = ham(); if (h < 6000) add[h] = 1;
        A[i][j] = a; A[j][i] = b;
      }
    }
    printf("n=%d one arc deleted: H=7 %s, H=21 %s; first absent values:", n, del[7] ? "yes" : "no",
           del[21] ? "yes" : "no");
    int shown = 0; for (int h = 0; h < 200 && shown < 12; h++) if (!del[h]) { printf(" %d", h); shown++; }
    printf("\nn=%d one reverse arc added: H=7 %s, H=21 %s; first absent values:", n, add[7] ? "yes" : "no",
           add[21] ? "yes" : "no");
    shown = 0; for (int h = 1; h < 200 && shown < 12; h++) if (!add[h]) { printf(" %d", h); shown++; }
    printf("\n");
  }
  return 0;
}
