/* A2 audit: windowless log-Fourier spectrum of class densities.
   mode 0: -1 basin on Z_<0 (3m-1 map on m>0); mode 1: positive entry class via 85.
   S(f) = sum_{m in [2^(k-1),2^k)} (1_B(m) - mean) * w_m * exp(-2 pi i f log2 m),  w_m = log2(1+1/m),
   amplitude = 2|S(f)| (comparable to the windowed 2|sum_j (d_j-mean) e^{..}|/W). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <string.h>
static uint8_t *lab;
static int tlab(int mode) { return mode == 0 ? 1 : 4; }
int main(int argc, char **argv) {
    int mode = atoi(argv[1]), KMAX = atoi(argv[2]), FMAX = atoi(argv[3]);
    long N = 1L << KMAX;
    lab = calloc(N + 2, 1);
    for (long m = 1; m <= N; m++) {
        if (mode == 0) {
            if (m == 1) { lab[1] = 1; continue; }
            if (m == 5) { lab[5] = 2; continue; }
            if (m == 17) { lab[17] = 3; continue; }
            unsigned __int128 x = m;
            while (1) {
                x = (x & 1) ? (3 * x - 1) / 2 : x / 2;
                if (x < (unsigned __int128)m) { lab[m] = lab[(long)x]; break; }
                if (x == (unsigned __int128)m) { fprintf(stderr, "cycle %ld\n", m); return 1; }
            }
        } else {
            if ((m & (m - 1)) == 0) { lab[m] = 0; continue; }       /* powers of 2 incl. 1 */
            unsigned __int128 v = m, u = 0;
            while (1) {
                unsigned __int128 prev = v;
                if (v & 1) v = (3 * v + 1) / 2; else v = v / 2;
                if ((v & (v - 1)) == 0) {               /* reached a power of 2: prev must be odd */
                    if (!(prev & 1)) { fprintf(stderr, "even prev?? %ld\n", m); return 1; }
                    unsigned __int128 t = 3 * prev + 1; int j = 0; while (t > 1) { t >>= 2; j++; }
                    lab[m] = (uint8_t)j; break;
                }
                if (v < (unsigned __int128)m) { lab[m] = lab[(long)v]; break; }
                (void)u;
            }
        }
    }
    int T = tlab(mode);
    for (int k = 16; k <= KMAX; k += 2) {
        long lo = 1L << (k - 1), hi = 1L << k;
        double tot = 0, cnt = 0, wsum = 0, wcnt = 0;
        for (long m = lo; m < hi; m++) { double w = log2(1.0 + 1.0 / m); wsum += w; tot++; if (lab[m] == T) { cnt++; wcnt += w; } }
        double mean = wcnt / wsum;   /* log-weighted mean */
        double complex *S = calloc(FMAX + 1, sizeof(double complex));
        for (long m = lo; m < hi; m++) {
            double w = log2(1.0 + 1.0 / m);
            double a = ((lab[m] == T) ? 1.0 : 0.0) - mean;
            if (a == 0) continue;
            double ph = log2((double)m) - (k - 1);       /* in [0,1) */
            double complex z = cexp(-2 * M_PI * I * ph), zz = 1;
            double aw = a * w;
            for (int f = 1; f <= FMAX; f++) { zz *= z; S[f] += aw * zz; }
        }
        printf("k=%d mode=%d natural_density=%.5f logmean=%.5f\n", k, mode, cnt / tot, mean);
        double *amp = malloc(sizeof(double) * (FMAX + 1));
        for (int f = 1; f <= FMAX; f++) amp[f] = 2 * cabs(S[f]);
        printf("  list:");
        int list[] = {5, 12, 17, 24, 29, 36, 41, 53, 65, 70, 94, 106, 118, 147, 159, 200, 212, 253, 265, 306, 318};
        for (int q = 0; q < 21; q++) if (list[q] <= FMAX) printf(" %d:%.4f", list[q], amp[list[q]]);
        printf("\n  top20:");
        for (int r = 0; r < 20; r++) { int b = 1; for (int f = 1; f <= FMAX; f++) if (amp[f] > amp[b]) b = f; printf(" %d:%.4f", b, amp[b]); amp[b] = -1; }
        printf("\n");
        free(S); free(amp);
    }
    return 0;
}
