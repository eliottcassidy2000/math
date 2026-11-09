/* S19 audit: find a labelled tournament of order N with port 0 whose base-2 port state is (h, R, D, hu); print its out-masks. */
#include <stdio.h>
#include "tcount.h"
int main(int argc, char **argv) {
  int n = atoi(argv[1]); ll th = atoll(argv[2]), tR = atoll(argv[3]), tD = atoll(argv[4]), tu = atoll(argv[5]);
  static ll I2[1 << MAXN]; int np = n * (n - 1) / 2; size_t FULL = ((size_t)1 << n) - 1;
  for (unsigned long long c = 0; c < (1ULL << np); c++) {
    uint32_t out[16] = {0}; int b = 0;
    for (int a = 0; a < n; a++) for (int d = a + 1; d < n; d++, b++) { if (c >> b & 1) out[a] |= 1u << d; else out[d] |= 1u << a; }
    TN = n; tc_alloc(n); for (int v = 0; v < n; v++) TOUT[v] = out[v];
    tc_cycles(); tc_packing(2, I2);
    if (I2[FULL] != th || I2[FULL & ~1] != tu) continue;
    ll E, O; tc_port(0, I2, &E, &O);
    if (E + O == tR && E - O == tD) { printf("%d 0", n); for (int v = 0; v < n; v++) printf(" %u", out[v]); printf("\n"); return 0; }
  }
  printf("none\n"); return 1;
}
