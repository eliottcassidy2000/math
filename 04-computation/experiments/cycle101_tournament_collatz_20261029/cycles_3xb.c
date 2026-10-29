/* Positive integer cycles of U_b(n) = (3n+b)/2^v (n odd, v = v_2(3n+b)), b odd > 0.
   Iterates every odd n <= N; a trajectory stops at a number already classified, at a repeat (new cycle), or
   above 2^61 (reported as escape).  Prints each cycle: min point, k (odd steps), H (halvings), max point,
   gcd(min, b), and the points.  usage: cycles_3xb b N */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned long long u64;
static u64 gcdu(u64 a, u64 b){ while(b){ u64 t=a%b; a=b; b=t;} return a; }
int main(int argc, char **argv){
  u64 b = strtoull(argv[1],0,10), N = strtoull(argv[2],0,10);
  /* status: 0 unknown, 1 known (reaches a known cycle or escapes) for odd n <= N, index n/2 */
  unsigned char *st = calloc(N/2+2, 1);
  u64 *path = malloc(sizeof(u64)*100000);
  long ncyc = 0, nesc = 0;
  for (u64 n0 = 1; n0 <= N; n0 += 2){
    if (st[n0/2]) continue;
    int len = 0; u64 n = n0; int stop = 0;
    while (1){
      if (n <= N && st[n/2]) { stop = 1; break; }
      for (int i = 0; i < len; i++) if (path[i] == n) { /* cycle found: path[i..len-1] */
          u64 mn = n, mx = n; int k = 0; int H = 0;
          for (int j = i; j < len; j++){ if (path[j] < mn) mn = path[j]; if (path[j] > mx) mx = path[j]; }
          /* recompute k, H around the cycle from its min */
          u64 m = mn; do { u64 t = 3*m + b; int v = 0; while (!(t & 1)){ t >>= 1; v++; } H += v; k++; m = t; } while (m != mn);
          printf("cycle min=%llu k=%d H=%d max=%llu gcd(min,b)=%llu points:", mn, k, H, mx, gcdu(mn,b));
          m = mn; int cnt = 0; do { printf(" %llu", m); u64 t = 3*m + b; while (!(t & 1)) t >>= 1; m = t; cnt++; } while (m != mn && cnt < 64);
          if (cnt >= 64) printf(" ...");
          printf("\n"); ncyc++; stop = 2; break; }
      if (stop) break;
      if (len >= 99999 || n > (1ULL<<61)) { nesc++; stop = 3; break; }
      path[len++] = n;
      u64 t = 3*n + b; while (!(t & 1)) t >>= 1; n = t;
    }
    for (int i = 0; i < len; i++) if (path[i] <= N) st[path[i]/2] = 1;
  }
  fprintf(stderr, "b=%llu N=%llu cycles=%ld escapes=%ld\n", b, N, ncyc, nesc);
  return 0;
}
