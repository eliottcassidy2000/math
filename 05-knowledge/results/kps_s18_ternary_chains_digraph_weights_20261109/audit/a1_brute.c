/* AUDIT (independent): brute-force port state of a ported tournament, all bases x = 2..8.
 * Input lines: "n port bits" (n*n adjacency bits, row-major), n <= 16.
 * Output lines: "H  h2 E2 O2  h3 E3 O3 ... h8 E8 O8" where H is the direct Hamiltonian path count and
 *   h_x = I(Omega(T), x), E_x / O_x = sum over paths P from the port with |P| even / odd of I(Omega(T - V(P)), x).
 * Method: directed odd Hamiltonian cycles per vertex set by a forward DP from the set minimum; packing counts
 * F[S] = F[S - m] + x * sum_{m in U subset S} cyc[U] F[S - U] (m = min S); paths from the port by a subset DP.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;
static int n, A[16][16];
static ll *cyc, *F, *dp;

int main(void) {
  char bits[300];
  int port;
  cyc = malloc(sizeof(ll) << 16);
  F = malloc(sizeof(ll) << 16);
  dp = malloc(sizeof(ll) * 16 * (1 << 16));
  while (scanf("%d %d %s", &n, &port, bits) == 3) {
    int N = 1 << n, full = N - 1;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] = bits[i * n + j] - '0';
    memset(cyc, 0, sizeof(ll) * N);
    for (int s = 0; s < n; s++) {
      /* subsets whose minimum is s; dp[S*n+v] = paths s..v covering S */
      for (int S = 0; S < N; S++) if ((S >> s & 1) && !(S & ((1 << s) - 1))) memset(dp + (ll)S * n, 0, sizeof(ll) * n);
      dp[(ll)(1 << s) * n + s] = 1;
      for (int S = 1 << s; S < N; S++) {
        if (!(S >> s & 1) || (S & ((1 << s) - 1))) continue;
        int pc = __builtin_popcount(S);
        for (int v = 0; v < n; v++) {
          ll c = dp[(ll)S * n + v];
          if (!c) continue;
          if (pc >= 3 && (pc & 1) && A[v][s]) cyc[S] += c;
          for (int t = s + 1; t < n; t++) if (!(S >> t & 1) && A[v][t]) dp[(ll)(S | 1 << t) * n + t] += c;
        }
      }
    }
    /* Hamiltonian paths */
    memset(dp, 0, sizeof(ll) * N * n);
    for (int v = 0; v < n; v++) dp[(ll)(1 << v) * n + v] = 1;
    for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) {
      ll c = dp[(ll)S * n + v];
      if (!c) continue;
      for (int t = 0; t < n; t++) if (!(S >> t & 1) && A[v][t]) dp[(ll)(S | 1 << t) * n + t] += c;
    }
    ll H = 0;
    for (int v = 0; v < n; v++) H += dp[(ll)full * n + v];
    printf("%lld", H);
    /* paths from the port: tot[S] stored in cyc-independent array (reuse dp) */
    memset(dp, 0, sizeof(ll) * N * n);
    dp[(ll)(1 << port) * n + port] = 1;
    for (int S = 1; S < N; S++) {
      if (!(S >> port & 1)) continue;
      for (int v = 0; v < n; v++) {
        ll c = dp[(ll)S * n + v];
        if (!c) continue;
        for (int t = 0; t < n; t++) if (!(S >> t & 1) && A[v][t]) dp[(ll)(S | 1 << t) * n + t] += c;
      }
    }
    for (ll x = 2; x <= 8; x++) {
      F[0] = 1;
      for (int S = 1; S < N; S++) {
        int m = __builtin_ctz(S), rest = S & ~(1 << m);
        ll val = F[rest];
        for (int T = rest;; T = (T - 1) & rest) {
          int U = T | (1 << m);
          if (cyc[U]) val += x * cyc[U] * F[S & ~U];
          if (!T) break;
        }
        F[S] = val;
      }
      ll E = 0, O = 0;
      for (int S = 1; S < N; S++) {
        if (!(S >> port & 1)) continue;
        ll tot = 0;
        for (int v = 0; v < n; v++) tot += dp[(ll)S * n + v];
        if (!tot) continue;
        if (__builtin_popcount(S) & 1) O += tot * F[full & ~S]; else E += tot * F[full & ~S];
      }
      printf("  %lld %lld %lld", F[full], E, O);
    }
    printf("\n");
    fflush(stdout);
  }
  return 0;
}
