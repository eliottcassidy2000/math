/* Inversion symmetry of positive rational 3x+1 periodic points.  For every odd b <= B, find all positive cycles
   of U_b(n) = (3n+b)/2^v with all points n <= N (trajectory search from odd n <= N); each point gives the
   rational n/b, reduced p/q.  A pair p/q, q/p is "mutual" if both are periodic.  Output: all periodic
   reduced rationals p/q with p, q <= L (as "P p q"), from which mutual pairs are read off. */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64;
static u64 g(u64 a, u64 b){ while(b){u64 t=a%b;a=b;b=t;} return a; }
int main(int argc, char **argv){
  u64 B = strtoull(argv[1],0,10), N = strtoull(argv[2],0,10), L = strtoull(argv[3],0,10);
  unsigned char *st = malloc(N/2+2);
  u64 *path = malloc(sizeof(u64)*200000);
  for (u64 b = 1; b <= B; b += 2){
    for (u64 i = 0; i <= N/2; i++) st[i] = 0;
    for (u64 n0 = 1; n0 <= N; n0 += 2){
      if (st[n0/2]) continue;
      int len = 0; u64 n = n0; int stop = 0;
      while (1){
        if (n <= N && st[n/2]) break;
        for (int i = 0; i < len; i++) if (path[i] == n){
            u64 m = n; do { u64 gg = g(m, b); u64 p = m/gg, q = b/gg; if (p <= L && q <= L) printf("P %llu %llu b=%llu\n", p, q, b);
                            u64 t = 3*m + b; while(!(t&1)) t >>= 1; m = t; } while (m != n);
            stop = 1; break; }
        if (stop) break;
        if (len >= 199999 || n > (1ULL<<61)) break;
        path[len++] = n;
        u64 t = 3*n + b; while (!(t & 1)) t >>= 1; n = t;
      }
      for (int i = 0; i < len; i++) if (path[i] <= N) st[path[i]/2] = 1;
    }
  }
  return 0;
}
