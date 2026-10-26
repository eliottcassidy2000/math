// E(X)/sqrt(X) at K log-spaced heights in [10^4, Xmax] (E = S - 1 - N_U - N_D, the remainder after the two root cusps),
// for the logarithmic mean and sign statistics in HYP-9418 (2).  Branch parity via the regular-CF digit automaton
// (THM-4808 (b)); triples bucketed by binary search over the checkpoints.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static int parity_cf(long m, long n) {
    long d[128]; int L = 0; long a = m, b = n;
    while (b) { d[L++] = a / b; long r = a % b; a = b; b = r; }
    if (a != 1) return -1;
    if (L > 1 && d[L - 1] == 1) { d[L - 2] += 1; L--; }
    int count = 0, i = 0; long head = d[0];
    while (1) {
        int rem = L - i;
        if (rem == 1) return count & 1;
        if (head % 2 == 0) { count++; i += 1; head = d[i]; }
        else { if (rem == 2) return count & 1; i += 2; head = d[i] + 1; }
    }
}
int main(int argc, char **argv) {
    double Xmax = atof(argv[1]); int K = atoi(argv[2]);
    double *cp = malloc(sizeof(double) * K);
    for (int k = 0; k < K; k++) cp[k] = floor(pow(10.0, 4.0 + (log10(Xmax) - 4.0) * k / (K - 1)));
    long long *Rb = calloc(K + 1, sizeof(long long));
    long M = (long)sqrt(Xmax) + 2;
    for (long m = 2; m <= M; m++) for (long n = (m % 2 ? 2 : 1); n < m; n += 2) {
        double c = (double)m * m + (double)n * n; if (c > Xmax) break;
        int par = parity_cf(m, n); if (par < 0) continue;
        int lo = 0, hi = K;                       // first checkpoint index with cp >= c
        while (lo < hi) { int mid = (lo + hi) / 2; if (cp[mid] >= c) hi = mid; else lo = mid + 1; }
        int ray = (n == 1) || (n == m - 1);
        Rb[lo] += (par ? -1 : 1) - (ray ? 1 : 0);
    }
    long long R = 0; double sumlog = 0, wsum = 0; int pos = 0, neg = 0, changes = 0, last = 0; double mx = 0;
    printf("X E(X) E/sqrtX\n");
    for (int k = 0; k < K; k++) {
        R += Rb[k];
        double v = R / sqrt(cp[k]);
        if (k % (K / 20 > 0 ? K / 20 : 1) == 0 || k == K - 1) printf("%.0f %lld %.5f\n", cp[k], R, v);
        if (k > 0) { double w = log(cp[k] / cp[k - 1]); sumlog += v * w; wsum += w; }
        int s = (R > 0) - (R < 0); if (s > 0) pos++; else if (s < 0) neg++;
        if (s != 0 && last != 0 && s != last) changes++; if (s != 0) last = s;
        if (fabs(v) > mx) mx = fabs(v);
    }
    printf("log-spaced heights %d in [1e4, %.0e]: logarithmic mean of E/sqrtX = %.5f; positive %d, negative %d; sign changes %d; max |E|/sqrtX = %.4f\n",
           K, Xmax, sumlog / wsum, pos, neg, changes, mx);
    return 0;
}
