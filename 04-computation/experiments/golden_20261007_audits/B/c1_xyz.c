/* THM-4566 addendum: (1) brute force which n <= N are NOT xy+yz+zx with x,y,z >= 1;
 * (2) independent idoneal test for every n = 2 mod 4 (n <= N): n idoneal iff no primitive reduced
 *     form (a, 2b, c), ac - b^2 = n, is non-ambiguous (0 < 2b < a < c, gcd(a,2b,c) = 1).
 * usage: c1_xyz N */
#include <stdio.h>
#include <stdlib.h>
static long g(long a, long b) { while (b) { long t = a % b; a = b; b = t; } return a; }
int main(int argc, char **argv) {
    long N = atol(argv[1]);
    char *rep = calloc(N + 1, 1), *bad = calloc(N + 1, 1);
    for (long x = 1; 3 * x * x <= N; x++)
        for (long y = x; ; y++) {
            long base = x * y, s = x + y;
            if (base + y * s > N) break;
            for (long n = base + y * s; n <= N; n += s) rep[n] = 1;
        }
    printf("non-representable n <= %ld:", N);
    long cnt = 0;
    for (long n = 1; n <= N; n++) if (!rep[n]) { printf(" %ld", n); cnt++; }
    printf("\ncount = %ld\n", cnt);
    /* non-ambiguous primitive reduced forms with ac - b^2 = n, n = 2 mod 4 */
    for (long a = 1; 3 * a * a <= 4 * N; a++)
        for (long b = 1; 2 * b < a; b++)
            for (long c = a + 1; a * c - b * b <= N; c++) {
                long n = a * c - b * b;
                if ((n & 3) != 2 || bad[n]) continue;
                if (g(g(a, 2 * b), c) == 1) bad[n] = 1;
            }
    printf("idoneal n = 2 mod 4, n <= %ld:", N);
    long ci = 0;
    for (long n = 2; n <= N; n += 4) if (!bad[n]) { printf(" %ld", n); ci++; }
    printf("\ncount = %ld\n", ci);
    /* comparison */
    long mism = 0;
    for (long n = 1; n <= N; n++) {
        int exc = !rep[n];
        int pred = (n == 1 || n == 4) || ((n & 3) == 2 && !bad[n]);
        if (exc != pred) { mism++; if (mism < 20) printf("MISMATCH n=%ld exc=%d pred=%d\n", n, exc, pred); }
    }
    printf("characterization mismatches for n <= %ld: %ld\n", N, mism);
    return 0;
}
