/* Basins of the three negative cycles of T(x) = x/2, (3x+1)/2 on Z_<0.
   For m = -n in [1, N]: basin id 1 (cycle -1), 5 (cycle -5,-7,-10), 17 (cycle -17,...).
   Work with m = |n|: T on negatives: n even -> n/2 ; n odd -> (3n+1)/2.  With n = -m:
   m even -> m/2 ; m odd -> (3m-1)/2  (the 3x-1 map on positives). Cycles of 3x-1: {1}, {5,7,10}, {17,25,37,55,82,41,61,91,136,68,34}. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main(int argc, char **argv) {
    long N = atol(argv[1]);
    unsigned char *b = calloc(N + 1, 1);
    /* seed cycles */
    long c5[] = {5, 7, 10}; long c17[] = {17,25,37,55,82,41,61,91,136,68,34};
    b[1] = 1; b[2] = 1; /* 2 -> 1 */
    for (int i = 0; i < 3; i++) b[c5[i]] = 5;
    for (int i = 0; i < 11; i++) b[c17[i]] = 17;
    for (long m = 1; m <= N; m++) {
        if (b[m]) continue;
        /* follow orbit until hitting a value < m (already known) or known */
        long x = m; int steps = 0;
        while (1) {
            x = (x % 2 == 0) ? x / 2 : (3 * x - 1) / 2;
            steps++;
            if (x <= N && b[x]) break;
            if (steps > 100000) { fprintf(stderr, "long orbit at %ld\n", m); return 1; }
        }
        unsigned char id = b[x];
        /* second pass: fill along orbit while within range */
        long y = m;
        while (!(y <= N && b[y])) { if (y <= N) b[y] = id; y = (y % 2 == 0) ? y / 2 : (3 * y - 1) / 2; }
    }
    /* stats by dyadic scale */
    printf("scale k: [2^(k-1), 2^k)  densities of basins (-1, -5, -17)   cut fraction (basin(m) != basin(m+1))   residue spread mod 3, 4, 9\n");
    for (int k = 2; (1L << k) <= N + 1; k++) {
        long lo = 1L << (k - 1), hi = (1L << k);
        long cnt[256] = {0}, cuts = 0, tot = 0;
        long r3[3][256] = {{0}}, r4[4][256] = {{0}}, r9[9][256] = {{0}};
        for (long m = lo; m < hi && m <= N; m++) {
            unsigned char id = b[m]; cnt[id]++; tot++;
            if (m + 1 <= N && b[m + 1] != id) cuts++;
            r3[m % 3][id]++; r4[m % 4][id]++; r9[m % 9][id]++;
        }
        double d1 = (double)cnt[1] / tot, d5 = (double)cnt[5] / tot, d17 = (double)cnt[17] / tot;
        /* residue spread: max over residues of |density of basin -1 in residue class - overall| */
        double s3 = 0, s4 = 0, s9 = 0;
        for (int r = 0; r < 3; r++) { long t = r3[r][1] + r3[r][5] + r3[r][17]; double d = (double)r3[r][1] / t; if (d - d1 > s3) s3 = d - d1; if (d1 - d > s3) s3 = d1 - d; }
        for (int r = 0; r < 4; r++) { long t = r4[r][1] + r4[r][5] + r4[r][17]; double d = (double)r4[r][1] / t; if (d - d1 > s4) s4 = d - d1; if (d1 - d > s4) s4 = d1 - d; }
        for (int r = 0; r < 9; r++) { long t = r9[r][1] + r9[r][5] + r9[r][17]; double d = (double)r9[r][1] / t; if (d - d1 > s9) s9 = d - d1; if (d1 - d > s9) s9 = d1 - d; }
        printf("k=%2d  %.4f %.4f %.4f   cuts %.4f   spread(-1 basin) mod3 %.4f mod4 %.4f mod9 %.4f\n", k, d1, d5, d17, (double)cuts / tot, s3, s4, s9);
    }
    return 0;
}
