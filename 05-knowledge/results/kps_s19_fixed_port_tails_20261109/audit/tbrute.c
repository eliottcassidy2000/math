/* S19 audit: brute-force statistics of an explicit tournament with a port.
 * stdin, one tournament per line:  n u out_0 out_1 ... out_{n-1}   (out_v = out-neighbour bitmask of v)
 * stdout, one line per input:
 *   H  hend_u  [I_x  E_x  O_x  hu_x  for x = 2, 3, 4]
 * H = Hamiltonian paths (direct DP), hend_u = Hamiltonian paths ending at u,
 * I_x = I(Omega(T), x), E_x/O_x = port sums, hu_x = I(Omega(T - u), x).
 */
#include <stdio.h>
#include "tcount.h"

int main(void) {
  int n, u;
  static ll Ix[1 << MAXN];
  while (scanf("%d %d", &n, &u) == 2) {
    for (int v = 0; v < n; v++) { unsigned long long m; if (scanf("%llu", &m) != 1) return 1; TOUT[v] = (uint32_t)m; }
    TN = n; tc_alloc(n);
    /* sanity: tournament */
    for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if (a != b) {
      int ab = TOUT[a] >> b & 1, ba = TOUT[b] >> a & 1;
      if (ab + ba != 1) { fprintf(stderr, "not a tournament\n"); return 2; }
    }
    ll hend[MAXN];
    ll H = tc_ham(hend);
    tc_cycles();
    printf("%lld %lld", H, hend[u]);
    size_t FULL = ((size_t)1 << n) - 1;
    for (ll x = 2; x <= 4; x++) {
      tc_packing(x, Ix);
      ll E, O; tc_port(u, Ix, &E, &O);
      printf(" %lld %lld %lld %lld", Ix[FULL], E, O, Ix[FULL & ~((size_t)1 << u)]);
    }
    printf("\n");
  }
  return 0;
}
