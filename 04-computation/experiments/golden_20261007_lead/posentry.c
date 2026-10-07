/* Positive integers: "entry class" e(n) = j where the last odd number > 1 in n's orbit is (4^j-1)/3 (j = 2: 5, j = 3: 21, ...),
   e(n) = 1 for powers of 2. Classes are closed under T away from the end. Density of e = 3 (entry via 21) and e >= 3,
   per dyadic scale, with Fourier analysis in log2 position (512 windows), and cut fraction (e(n) != e(n+1)). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(int argc, char **argv) {
    long N = atol(argv[1]); int W = 512;
    unsigned char *e = calloc(N + 1, 1);
    e[1] = 1;
    for (long m = 2; m <= N; m++) {
        if (e[m]) continue;
        /* follow orbit (standard 3x+1 with halving) until reaching a known value; remember last odd > 1 if reaching 1 directly */
        long x = m; unsigned char id = 0;
        while (1) {
            if (x <= N && e[x]) { id = e[x]; break; }
            if (x % 2) {
                long y = 3 * x + 1;
                if ((y & (y - 1)) == 0) { /* 3x+1 is a power of 2: x is the last odd */
                    int jj = 0; long t = y; while (t > 1) { t >>= 2; jj++; } id = (unsigned char)jj; break; }
                x = y;
            } else x /= 2;
        }
        long y = m;
        while (!(y <= N && e[y])) {
            if (y <= N) e[y] = id;
            if (y % 2) { long z = 3 * y + 1; if ((z & (z - 1)) == 0) break; y = z; } else y /= 2;
        }
    }
    double *dens = malloc(sizeof(double) * W);
    for (int k = 16; (1L << k) <= N + 1; k += 2) {
        long lo = 1L << (k - 1), cnt3 = 0, cnt2 = 0, tot = 0, cuts = 0;
        for (long m = lo; m < 2 * lo; m++) { tot++; if (e[m] == 3) cnt3++; if (e[m] == 2) cnt2++; if (m + 1 <= N && e[m] != e[m + 1]) cuts++; }
        double mean = 0;
        for (int j = 0; j < W; j++) {
            long a = (long)ceil(lo * pow(2.0, (double)j / W)), z = (long)ceil(lo * pow(2.0, (double)(j + 1) / W));
            long c = 0, t = 0;
            for (long m = a; m < z; m++) { t++; if (e[m] == 3) c++; }
            dens[j] = t ? (double)c / t : 0; mean += dens[j] / W;
        }
        printf("k=%d  P(entry via 5)=%.4f P(via 21)=%.4f  cuts %.4f  top modes of P(via 21) (freq/octave:amp):", k, (double)cnt2 / tot, (double)cnt3 / tot, (double)cuts / tot);
        double amp[41];
        for (int f = 1; f <= 40; f++) {
            double re = 0, im = 0;
            for (int j = 0; j < W; j++) { double th = 2 * M_PI * f * (j + 0.5) / W; re += (dens[j] - mean) * cos(th); im += (dens[j] - mean) * sin(th); }
            amp[f] = 2 * sqrt(re * re + im * im) / W;
        }
        for (int r = 0; r < 5; r++) { int best = 1; for (int f = 1; f <= 40; f++) if (amp[f] > amp[best]) best = f; printf("  %d:%.4f", best, amp[best]); amp[best] = -1; }
        printf("\n");
    }
    return 0;
}
