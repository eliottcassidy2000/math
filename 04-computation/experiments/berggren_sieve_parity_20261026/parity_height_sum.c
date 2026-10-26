// Signed count of the branch character by height: S(X) = sum over primitive triples with hypotenuse c <= X of
// (-1)^#A, where #A = number of middle Berggren steps; also the same sum over prime hypotenuses.
// Enumerates Euclid pairs m > n >= 1, coprime, opposite parity, m^2 + n^2 <= X; #A by the descent.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static long gcdl(long a, long b) { while (b) { long t = a % b; a = b; b = t; } return a; }
int main(int argc, char **argv) {
    long X = atol(argv[1]);
    long M = (long)sqrt((double)X) + 2;
    // prime sieve up to X for prime-hypotenuse sums
    char *comp = calloc(X + 1, 1);
    for (long i = 2; i * i <= X; i++) if (!comp[i]) for (long j = i * i; j <= X; j += i) comp[j] = 1;
    long checkpoints[64]; int nc = 0;
    for (long c = 10; c <= X; c *= 10) checkpoints[nc++] = c;
    // accumulate per hypotenuse value in an array of signed counts (multiple triples can share c)
    int *sgn = calloc(X + 1, sizeof(int)); int *cnt = calloc(X + 1, sizeof(int));
    for (long m = 2; m <= M; m++) for (long n = (m % 2 ? 2 : 1); n < m; n += 2) {
        long c = m * m + n * n; if (c > X) break;
        if (gcdl(m, n) != 1) continue;
        long a = m, b = n; int nA = 0;
        while (!(a == 2 && b == 1)) {
            if (a > 3 * b) a -= 2 * b;
            else if (a > 2 * b) { long t = a - 2 * b; a = b; b = t; nA++; }
            else { long t = 2 * b - a; a = b; b = t; }
        }
        sgn[c] += (nA & 1) ? -1 : 1; cnt[c] += 1;
    }
    long S = 0, N = 0, Sp = 0, Np = 0; int k = 0;
    printf("X, N(X) triples, S(X)=sum (-1)^#A, S/sqrt(N), prime-hypotenuse: Np, Sp, Sp/sqrt(Np)\n");
    for (long c = 1; c <= X; c++) {
        S += sgn[c]; N += cnt[c];
        if (c > 1 && !comp[c]) { Sp += sgn[c]; Np += cnt[c]; }
        if (k < nc && c == checkpoints[k]) {
            printf("%ld %ld %ld %.4f %ld %ld %.4f\n", c, N, S, N ? S / sqrt((double)N) : 0.0, Np, Sp, Np ? Sp / sqrt((double)Np) : 0.0);
            k++;
        }
    }
    return 0;
}
