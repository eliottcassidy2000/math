// pi2(X) = #{p <= X : p and p + 2 prime} at X = 10^3..10^K by a segmented odd-only sieve, and the Mertens ratio
//   rho(X) = (X/2) * prod_{3 <= q <= sqrt X} (1 - 2/q) / pi2(X)
// (the naive independent-residue density of the exact level-sqrt(X) twin sieve, over the true count),
// compared with its limit 4 e^(-2 gamma) and with the Hardy-Littlewood model 4 e^(-2 gamma) (X / log^2 X) / Li2(X).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
int main(int argc, char **argv) {
    int K = atoi(argv[1]);
    long X = 1; for (int i = 0; i < K; i++) X *= 10;
    long R = (long)sqrtl((long double)(X + 2)) + 2;
    char *small = calloc(R + 1, 1);
    for (long i = 2; i * i <= R; i++) if (!small[i]) for (long j = i * i; j <= R; j += i) small[j] = 1;
    long *pr = malloc(sizeof(long) * (R + 1)); long np = 0;
    for (long i = 3; i <= R; i++) if (!small[i]) pr[np++] = i;
    long bucket[32] = {0};
    const long SEG = 1L << 24;
    char *seg = malloc(SEG);
    long lastp = 2;
    for (long lo = 3; lo <= X + 2; lo += 2 * SEG) {
        long hi = lo + 2 * SEG; if (hi > X + 3) hi = X + 3;           // odd numbers lo, lo+2, ..., < hi
        long cnt = (hi - lo + 1) / 2;
        memset(seg, 0, cnt);
        for (long t = 0; t < np; t++) {
            long q = pr[t]; if (q * q >= hi) break;
            long start = q * q;
            if (start < lo) start = ((lo + q - 1) / q) * q;
            if (start % 2 == 0) start += q;
            for (long j = start; j < hi; j += 2 * q) seg[(j - lo) >> 1] = 1;
        }
        for (long i = 0; i < cnt; i++) {
            if (seg[i]) continue;
            long p = lo + 2 * i;
            if (p - lastp == 2 && lastp <= X) {
                int b = 0; long c = 1; while (c < lastp) { c *= 10; b++; }   // smallest b with lastp <= 10^b
                bucket[b]++;
            }
            lastp = p;
        }
    }
    const long double g = 0.57721566490153286060651209L;
    const long double limit = 4.0L * expl(-2.0L * g);
    printf("X pi2(X) prod_{3<=q<=sqrtX}(1-2/q) rho(X) HL-Mertens-model limit\n");
    long cum = 0; long Xc = 1;
    for (int b = 0; b <= K; b++) {
        cum += bucket[b];
        if (b >= 3) {
            long double prod = 1.0L; long double s = sqrtl((long double)Xc);
            for (long t = 0; t < np && pr[t] <= s; t++) prod *= (1.0L - 2.0L / pr[t]);
            long double rho = ((long double)Xc / 2.0L) * prod / cum;
            // Li2(X) = int_2^X dt / log^2 t by Simpson on log scale
            long double a = logl(2.0L), bb = logl((long double)Xc); int N = 200000; long double h = (bb - a) / N, S = 0;
            for (int i = 0; i <= N; i++) { long double u = a + i * h; long double f = expl(u) / (u * u); S += f * ((i == 0 || i == N) ? 1 : (i % 2 ? 4 : 2)); }
            long double Li2 = S * h / 3.0L;
            long double model = limit * ((long double)Xc / (bb * bb)) / Li2;
            printf("10^%d %ld %.10Lf %.5Lf %.5Lf %.5Lf\n", b, cum, prod, rho, model, limit);
        }
        Xc *= 10;
    }
    return 0;
}
