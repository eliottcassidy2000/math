/* Independent packing counts for small tournaments (order <= 16), used to check THM-5070.
 * Input (stdin): lines "n row0row1...row(n-1)" with the n*n adjacency bits in row order.
 * Output: "ham I2 I3" per line: direct Hamiltonian path count, I(Omega,2) and I(Omega,3), where
 * c[S] = number of directed Hamiltonian cycles on vertex set S (|S| odd >= 3), and
 * I(x) = sum over set partitions into packed odd-cycle blocks: F[S] = F[S - min] + sum_{T ∋ min} x c[T] F[S-T].
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n, A[16][16];
static long long *cyc, *pathdp, *F;

int main(void) {
  char buf[400];
  cyc = malloc(sizeof(long long) << 16);
  F = malloc(sizeof(long long) << 16);
  pathdp = malloc(sizeof(long long) * (1 << 16) * 16);
  while (scanf("%d %s", &n, buf) == 2) {
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] = buf[i * n + j] - '0';
    int N = 1 << n;
    /* Hamiltonian cycles per subset: fix the minimum vertex s as start */
    for (int S = 0; S < N; S++) cyc[S] = 0;
    for (int s = 0; s < n; s++) {
      int base = 1 << s;
      /* dp over subsets of {s..n-1} containing s, paths from s */
      for (int S = 0; S < N; S++) if ((S & base) && !(S & (base - 1))) for (int v = 0; v < n; v++) pathdp[(long long)S * n + v] = 0;
      pathdp[(long long)base * n + s] = 1;
      for (int S = 0; S < N; S++) {
        if (!(S & base) || (S & (base - 1))) continue;
        for (int v = 0; v < n; v++) {
          long long c = pathdp[(long long)S * n + v];
          if (!c) continue;
          for (int w = s + 1; w < n; w++) if (!(S >> w & 1) && A[v][w]) pathdp[(long long)(S | 1 << w) * n + w] += c;
          if (A[v][s] && __builtin_popcount(S) >= 3) cyc[S] += c;
        }
      }
    }
    long long res[2];
    for (int xi = 0; xi < 2; xi++) {
      long long x = xi + 2;
      F[0] = 1;
      for (int S = 1; S < N; S++) {
        int m = __builtin_ctz(S);
        long long val = F[S & ~(1 << m)];
        int rest = S & ~(1 << m);
        for (int T = rest; ; T = (T - 1) & rest) {
          int U = T | (1 << m);
          if (__builtin_popcount(U) % 2 == 1 && __builtin_popcount(U) >= 3 && cyc[U]) val += x * cyc[U] * F[S & ~U];
          if (T == 0) break;
        }
        F[S] = val;
      }
      res[xi] = F[N - 1];
    }
    /* direct Hamiltonian path count */
    long long *dp = pathdp;
    for (long long i = 0; i < (long long)N * n; i++) dp[i] = 0;
    for (int v = 0; v < n; v++) dp[(long long)(1 << v) * n + v] = 1;
    for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) {
      long long c = dp[(long long)S * n + v]; if (!c) continue;
      for (int w = 0; w < n; w++) if (!(S >> w & 1) && A[v][w]) dp[(long long)(S | 1 << w) * n + w] += c;
    }
    long long h = 0; for (int v = 0; v < n; v++) h += dp[(long long)(N - 1) * n + v];
    printf("%lld %lld %lld\n", h, res[0], res[1]);
  }
  return 0;
}
