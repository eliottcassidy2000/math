/* Multiplicative-window and large-modulus equidistribution of the three negative basins at scale [2^27, 2^28). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(int argc, char **argv) {
    long N = atol(argv[1]);
    unsigned char *b = calloc(N + 1, 1);
    long c5[] = {5, 7, 10}; long c17[] = {17,25,37,55,82,41,61,91,136,68,34};
    b[1] = 1; b[2] = 1;
    for (int i = 0; i < 3; i++) b[c5[i]] = 5;
    for (int i = 0; i < 11; i++) b[c17[i]] = 17;
    for (long m = 1; m <= N; m++) {
        if (b[m]) continue;
        long x = m; while (1) { x = (x % 2 == 0) ? x / 2 : (3 * x - 1) / 2; if (x <= N && b[x]) break; }
        unsigned char id = b[x]; long y = m;
        while (!(y <= N && b[y])) { if (y <= N) b[y] = id; y = (y % 2 == 0) ? y / 2 : (3 * y - 1) / 2; }
    }
    long lo = (N + 1) / 2, hi = N + 1;  /* [2^27, 2^28) when N = 2^28 - 1 */
    /* 1. multiplicative windows: 64 log-spaced windows across [lo, hi) */
    printf("window j: y = lo*2^(j/64), eta = 2^(1/64)-1; densities (-1, -5, -17)\n");
    double mn[3] = {1,1,1}, mx[3] = {0,0,0};
    for (int j = 0; j < 64; j++) {
        long a = (long)(lo * pow(2.0, j / 64.0)), z = (long)(lo * pow(2.0, (j + 1) / 64.0));
        long c1 = 0, c5c = 0, c17c = 0, t = 0;
        for (long m = a; m < z && m < hi; m++) { t++; if (b[m] == 1) c1++; else if (b[m] == 5) c5c++; else c17c++; }
        double d[3] = {(double)c1 / t, (double)c5c / t, (double)c17c / t};
        for (int i = 0; i < 3; i++) { if (d[i] < mn[i]) mn[i] = d[i]; if (d[i] > mx[i]) mx[i] = d[i]; }
        if (j % 8 == 0) printf("  j=%2d  %.4f %.4f %.4f  (n=%ld)\n", j, d[0], d[1], d[2], t);
    }
    printf("range over 64 windows: -1 [%.4f, %.4f]  -5 [%.4f, %.4f]  -17 [%.4f, %.4f]\n", mn[0], mx[0], mn[1], mx[1], mn[2], mx[2]);
    /* 2. large moduli: max deviation of basin(-1) density over residue classes mod M */
    long Ms[] = {16, 64, 256, 1024, 27, 81, 243, 729, 11, 121, 7, 49};
    for (int q = 0; q < 12; q++) {
        long M = Ms[q]; long *tot = calloc(M, sizeof(long)), *one = calloc(M, sizeof(long)); long T1 = 0, TT = 0;
        for (long m = lo; m < hi; m++) { tot[m % M]++; TT++; if (b[m] == 1) { one[m % M]++; T1++; } }
        double d1 = (double)T1 / TT, dev = 0;
        for (long r = 0; r < M; r++) { double d = (double)one[r] / tot[r]; if (fabs(d - d1) > dev) dev = fabs(d - d1); }
        printf("mod %4ld: max |density(-1 basin | r) - overall| = %.4f  (expected sampling sd ~ %.4f)\n", M, dev, sqrt(d1 * (1 - d1) * M / (double)TT));
        free(tot); free(one);
    }
    return 0;
}
