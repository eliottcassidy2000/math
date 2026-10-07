/* A2 audit: independent negative-integer basin labels for T on Z_<0, written as the 3m-1 map on m=-n>0.
   Labelling by first descent below m (stopping-time recursion), NOT by path filling.
   Output: per-octave densities, cut fraction, residue deviation of -1 basin mod 3,4,9;
   log-Fourier spectrum (W windows/octave, direct DFT) at chosen k; window range at top octave. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
static uint8_t *lab;
int main(int argc, char **argv) {
    int KMAX = atoi(argv[1]);
    int W = argc > 2 ? atoi(argv[2]) : 2048;
    long N = (1L << KMAX);              /* labels for m in [1, N] */
    lab = calloc(N + 2, 1);
    if (!lab) { fprintf(stderr, "alloc\n"); return 1; }
    long unknown_cycles = 0;
    for (long m = 1; m <= N; m++) {
        if (m == 1) { lab[1] = 1; continue; }
        if (m == 5) { lab[5] = 2; continue; }
        if (m == 17) { lab[17] = 3; continue; }
        unsigned __int128 x = m; long steps = 0;
        while (1) {
            if (x & 1) x = (3 * x - 1) / 2; else x = x / 2;
            steps++;
            if (x < (unsigned __int128)m) { lab[m] = lab[(long)x]; break; }
            if (x == (unsigned __int128)m) { unknown_cycles++; lab[m] = 9; fprintf(stderr, "new cycle min %ld\n", m); break; }
            if (steps > 10000000) { fprintf(stderr, "runaway %ld\n", m); return 1; }
        }
    }
    printf("unknown cycles: %ld\n", unknown_cycles);
    printf(" k   d(-1)   d(-5)   d(-17)  cuts    dev3    dev4    dev9\n");
    for (int k = 2; k <= KMAX; k++) {
        long lo = 1L << (k - 1), hi = 1L << k;
        long c[4] = {0}, cuts = 0, tot = 0;
        long r3[3] = {0}, n3[3] = {0}, r4[4] = {0}, n4[4] = {0}, r9[9] = {0}, n9[9] = {0};
        for (long m = lo; m < hi; m++) {
            int L = lab[m]; c[L]++; tot++;
            if (lab[m + 1] != L) cuts++;          /* m+1 <= N always since hi-1+1 = hi <= N */
            n3[m % 3]++; n4[m % 4]++; n9[m % 9]++;
            if (L == 1) { r3[m % 3]++; r4[m % 4]++; r9[m % 9]++; }
        }
        double d1 = (double)c[1] / tot, dv3 = 0, dv4 = 0, dv9 = 0;
        for (int r = 0; r < 3; r++) { double d = fabs((double)r3[r] / n3[r] - d1); if (d > dv3) dv3 = d; }
        for (int r = 0; r < 4; r++) { double d = fabs((double)r4[r] / n4[r] - d1); if (d > dv4) dv4 = d; }
        for (int r = 0; r < 9; r++) { double d = fabs((double)r9[r] / n9[r] - d1); if (d > dv9) dv9 = d; }
        printf("%2d  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n", k, d1, (double)c[2] / tot, (double)c[3] / tot,
               (double)cuts / tot, dv3, dv4, dv9);
    }
    /* Fourier: windows j: [lo*2^(j/W), lo*2^((j+1)/W)), density of label 1; also labels 2,3 */
    double *dens = malloc(sizeof(double) * W);
    for (int k = 12; k <= KMAX; k += 2) {
        long lo = 1L << (k - 1);
        for (int which = 1; which <= 1; which++) {
            double mean = 0;
            for (int j = 0; j < W; j++) {
                long a = (long)ceil(lo * exp2((double)j / W)), z = (long)ceil(lo * exp2((double)(j + 1) / W));
                long cc = 0, t = 0;
                for (long m = a; m < z; m++) { t++; if (lab[m] == which) cc++; }
                dens[j] = t ? (double)cc / t : 0; mean += dens[j] / W;
            }
            printf("FOURIER k=%d W=%d label=%d mean=%.4f :", k, W, which, mean);
            int FM = 140; double amp[141], ph[141];
            for (int f = 1; f <= FM; f++) {
                double re = 0, im = 0;
                for (int j = 0; j < W; j++) { double th = 2 * M_PI * f * (j + 0.5) / W; re += (dens[j] - mean) * cos(th); im -= (dens[j] - mean) * sin(th); }
                amp[f] = 2 * sqrt(re * re + im * im) / W; ph[f] = atan2(im, re);
            }
            /* print all amps for f in list, and top 12 */
            int list[] = {5, 7, 12, 17, 24, 29, 36, 41, 46, 53, 58, 65, 70, 94, 106, 118};
            for (int q = 0; q < 16; q++) printf(" %d:%.4f", list[q], amp[list[q]]);
            printf("\n  top:");
            double tmp[141]; memcpy(tmp, amp, sizeof(amp));
            for (int r = 0; r < 14; r++) { int b = 1; for (int f = 1; f <= FM; f++) if (tmp[f] > tmp[b]) b = f; printf(" %d:%.4f", b, tmp[b]); tmp[b] = -1; }
            printf("\n  phases(12,53,41,65): %.3f %.3f %.3f %.3f\n", ph[12], ph[53], ph[41], ph[65]);
        }
    }
    /* multiplicative windows of 1/64 octave at top octave */
    {
        long lo = 1L << (KMAX - 1); double mn = 1, mx = 0;
        for (int j = 0; j < 64; j++) {
            long a = (long)(lo * exp2(j / 64.0)), z = (long)(lo * exp2((j + 1) / 64.0));
            long cc = 0, t = 0;
            for (long m = a; m < z; m++) { t++; if (lab[m] == 1) cc++; }
            double d = (double)cc / t; if (d < mn) mn = d; if (d > mx) mx = d;
        }
        printf("top octave k=%d, 64 windows: -1 basin range [%.4f, %.4f]\n", KMAX, mn, mx);
    }
    /* write labels for octave KMAX-? to file for python FFT check */
    FILE *fo = fopen(argv[3] ? argv[3] : "/dev/null", "wb");
    if (fo) { fwrite(lab, 1, N + 1, fo); fclose(fo); }
    return 0;
}
