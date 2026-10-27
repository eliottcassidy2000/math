/* independent: max over primitive odd s < N of sum of native prices lambda = 2^a/(6X + 2^a) along the orbit to 1,
   and of the forcing sum 2 artanh(lambda) = log(1 + 2^a/(3X)); the root is computed by valuations (3z+1 = 4^j (3y+1)). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef unsigned __int128 u128;
static int v2(u128 x) { int v = 0; while (!(x & 1)) { x >>= 1; v++; } return v; }
int main(int argc, char **argv) {
    unsigned long long N = strtoull(argv[1], 0, 10);
    long double best = 0, bestF = 0; unsigned long long arg = 0, argF = 0;
    for (unsigned long long s0 = 3; s0 < N; s0 += 2) {
        if (s0 % 8 == 5) continue;
        u128 s = s0; long double P = 0, F = 0; int steps = 0;
        while (s != 1) {
            u128 X = 3 * s + 1;
            int a = (s % 8 == 1) ? 2 : 1;
            u128 z = X >> a;               /* odd */
            u128 w3 = 3 * z + 1; int w = v2(w3); int j = (w - 1) / 2;
            u128 Y = w3 >> (2 * j);        /* 3y + 1 */
            u128 y = (Y - 1) / 3;
            long double lam = (long double)(1 << a) / (6.0L * (long double)X + (long double)(1 << a));
            P += lam;
            F += log1pl((long double)(1 << a) / (3.0L * (long double)X));
            s = y; steps++;
            if (steps > 100000) { printf("long orbit %llu\n", s0); return 1; }
        }
        if (P > best) { best = P; arg = s0; }
        if (F > bestF) { bestF = F; argF = s0; }
    }
    printf("N=%llu: max sum lambda = %.12Lf at s=%llu; max sum 2artanh(lambda) = %.12Lf at s=%llu\n", N, best, arg, bestF, argF);
    return 0;
}
