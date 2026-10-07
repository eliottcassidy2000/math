/* Log-periodic structure of the -1 basin density: for each dyadic scale k (16..28), split [2^(k-1), 2^k) into W=512
   log-equal windows, compute the -1-basin density per window, and print Fourier amplitudes |c_f| for f = 1..40
   (frequency f = f cycles per octave). Also the -5 and -17 basins' dominant modes. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(int argc, char **argv) {
    long N = atol(argv[1]); int W = 2048;
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
    double *dens = malloc(sizeof(double) * W);
    for (int k = 16; (1L << k) <= N + 1; k += 2) {
        long lo = 1L << (k - 1);
        double mean = 0;
        for (int j = 0; j < W; j++) {
            long a = (long)ceil(lo * pow(2.0, (double)j / W)), z = (long)ceil(lo * pow(2.0, (double)(j + 1) / W));
            long c = 0, t = 0;
            for (long m = a; m < z; m++) { t++; if (b[m] == 1) c++; }
            dens[j] = t ? (double)c / t : 0; mean += dens[j] / W;
        }
        printf("k=%d mean(-1 basin)=%.4f  top Fourier modes (freq per octave: amplitude):", k, mean);
        double amp[71];
        for (int f = 1; f <= 70; f++) {
            double re = 0, im = 0;
            for (int j = 0; j < W; j++) { double th = 2 * M_PI * f * (j + 0.5) / W; re += (dens[j] - mean) * cos(th); im += (dens[j] - mean) * sin(th); }
            amp[f] = 2 * sqrt(re * re + im * im) / W;
        }
        /* print top 6 */
        for (int r = 0; r < 9; r++) { int best = 1; for (int f = 1; f <= 70; f++) if (amp[f] > amp[best]) best = f; printf("  %d:%.4f", best, amp[best]); amp[best] = -1; }
        printf("\n");
    }
    return 0;
}
