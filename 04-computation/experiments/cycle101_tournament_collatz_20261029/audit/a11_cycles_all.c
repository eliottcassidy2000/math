/* Audit: positive integer cycles of n -> (3n+b)/2^v (n odd) found from every odd start n <= N.
   Independent implementation (Brent cycle detection per start; a bitmap marks odd values <= N already
   known to lead to a recorded cycle, so later starts stop there).
   Prints, per cycle: min, k (odd steps), H (halvings), and the points in orbit order from the min.
   usage: a7_cycles b N */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned long long u64;
static u64 b;
static inline u64 U(u64 n, int *v) { u64 t = 3 * n + b; int c = __builtin_ctzll(t); *v = c; return t >> c; }
#define MAXC 4096
static u64 cyc_min[MAXC]; static int ncyc = 0;
int main(int argc, char **argv) {
  b = strtoull(argv[1], 0, 10); u64 N = strtoull(argv[2], 0, 10);
  unsigned char *known = calloc(N / 2 + 2, 1);
  long escapes = 0;
  for (u64 n0 = 1; n0 <= N; n0 += 2) {
    if (known[n0 / 2]) continue;
    /* Brent: find a point on the eventual cycle (or hit a known value / escape) */
    u64 power = 1, lam = 1, tort = n0, hare; int v; int stop = 0;
    hare = U(n0, &v);
    long steps = 0;
    while (tort != hare) {
      if (hare <= N && known[hare / 2]) { stop = 1; break; }
      if (hare > (1ULL << 60) || steps > 50000000L) { stop = 2; break; }
      if (power == lam) { tort = hare; power *= 2; lam = 0; }
      hare = U(hare, &v); lam++; steps++;
    }
    u64 m = 0;
    if (stop == 0) {
      /* hare is on a cycle: get its min, record if new */
      m = hare; u64 x = hare; do { x = U(x, &v); if (x < m) m = x; } while (x != hare);
      int seen = 0; for (int i = 0; i < ncyc; i++) if (cyc_min[i] == m) seen = 1;
      if (!seen && ncyc < MAXC) {
        cyc_min[ncyc++] = m;
        int k = 0, H = 0; x = m;
        do { x = U(x, &v); k++; H += v; } while (x != m);
        printf("b=%llu cycle min=%llu k=%d H=%d points:", b, m, k, H);
        x = m; do { printf(" %llu", x); x = U(x, &v); } while (x != m);
        printf("\n");
      }
    } else if (stop == 2) escapes++;
    /* mark the whole trajectory up to the first known value or once around the cycle */
    u64 x = n0; long guard = 0;
    while (1) {
      if (x <= N) { if (known[x / 2]) break; known[x / 2] = 1; }
      if (stop == 0 && x == m) { /* on the cycle: mark all its points <= N, then stop */
        u64 y = U(m, &v); while (y != m) { if (y <= N) known[y / 2] = 1; y = U(y, &v); } break; }
      if (x > (1ULL << 60) || ++guard > 50000000L) break;
      x = U(x, &v);
    }
  }
  printf("b=%llu N=%llu cycles=%d escapes=%ld\n", b, N, ncyc, escapes);
  return 0;
}
