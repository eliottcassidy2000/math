// R(X) = S(X) - (# root-ray triples U^k, D^k with c <= X), where S(X) = sum_{c <= X} (-1)^#A over primitive triples.
// Prints R(X), R/sqrt(X), R/X^(1/4) at X = 10^k * {1, 2, 5}.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static long gcdl(long a, long b) { while (b) { long t = a % b; a = b; b = t; } return a; }
int main(int argc, char **argv) {
    long X = atol(argv[1]);
    long M = (long)sqrt((double)X) + 2;
    int *sgn = calloc(X + 1, sizeof(int));
    for (long m = 2; m <= M; m++) for (long n = (m % 2 ? 2 : 1); n < m; n += 2) {
        long c = m * m + n * n; if (c > X) break;
        if (gcdl(m, n) != 1) continue;
        long a = m, b = n; int nA = 0;
        while (!(a == 2 && b == 1)) {
            if (a > 3 * b) a -= 2 * b;
            else if (a > 2 * b) { long t = a - 2 * b; a = b; b = t; nA++; }
            else { long t = 2 * b - a; a = b; b = t; }
        }
        // root rays: D^k: (m, 1) with m even;  U^k: (m, m-1)
        int ray = (n == 1) || (n == m - 1);
        sgn[c] += ((nA & 1) ? -1 : 1) - (ray ? 1 : 0);
    }
    long R = 0; double base[3] = {1, 2, 5};
    printf("X R(X) R/sqrt(X) R/X^0.25 R/X^0.4\n");
    for (long p = 10; p <= X; p *= 10) for (int i = 0; i < 3; i++) {}
    long next[64]; int nn = 0;
    for (long p = 100; p <= X; p *= 10) for (int i = 0; i < 3; i++) { long v = (long)(p * base[i]); if (v <= X) next[nn++] = v; }
    int k = 0;
    for (long c = 1; c <= X && k < nn; c++) {
        R += sgn[c];
        if (c == next[k]) { printf("%ld %ld %.5f %.3f %.4f\n", c, R, R / sqrt((double)c), R / pow((double)c, 0.25), R / pow((double)c, 0.4)); k++; }
    }
    return 0;
}
