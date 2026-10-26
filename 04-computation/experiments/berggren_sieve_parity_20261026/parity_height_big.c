// Branch-character height sums S(X) = sum_{primitive triples, c <= X} (-1)^#A and the root-ray remainder R(X), for
// checkpoints up to Xmax without per-c arrays: each Euclid pair (m, n) is bucketed by checkpoint; #A comes from the
// regular continued fraction of m/n through the digit automaton of THM-4808 (P1'), which is cross-checked against the
// literal descent on every pair with m < 4000.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static int nA_descent(long a, long b) {
    int nA = 0;
    while (!(a == 2 && b == 1)) {
        if (a > 3 * b) { long k = (a - 2 * b - 1) / (2 * b); if (k < 1) k = 1; a -= 2 * k * b; }
        else if (a > 2 * b) { long t = a - 2 * b; a = b; b = t; nA++; }
        else { long t = 2 * b - a; a = b; b = t; }
    }
    return nA;
}
// returns -1 if gcd(m, n) != 1, else #A mod 2 via the CF digit automaton
static int parity_cf(long m, long n) {
    long d[128]; int L = 0; long a = m, b = n;
    while (b) { d[L++] = a / b; long r = a % b; a = b; b = r; }
    if (a != 1) return -1;
    if (L > 1 && d[L - 1] == 1) { d[L - 2] += 1; L--; }
    int count = 0, i = 0; long head = d[0];
    while (1) {
        int rem = L - i;                         // digits remaining including head
        if (rem == 1) return count & 1;          // even integer slope: only D steps remain (or the root)
        if (head % 2 == 0) { count++; i += 1; head = d[i]; }
        else {
            if (rem == 2) return count & 1;      // [h; c1] -> ... -> [1; 1] = root
            i += 2; head = d[i] + 1;
        }
    }
}
int main(int argc, char **argv) {
    double Xmax = atof(argv[1]);
    double cps[] = {1e8, 2e8, 5e8, 1e9, 2e9, 5e9, 1e10};
    int ncp = 0; while (ncp < 7 && cps[ncp] <= Xmax) ncp++;
    long long Sb[8] = {0}, Nb[8] = {0}, Rb[8] = {0};
    long M = (long)sqrt(Xmax) + 2; long mism = 0, checked = 0;
    for (long m = 2; m <= M; m++) {
        for (long n = (m % 2 ? 2 : 1); n < m; n += 2) {
            double c = (double)m * m + (double)n * n; if (c > Xmax) break;
            int par = parity_cf(m, n); if (par < 0) continue;
            if (m < 4000) { checked++; if ((nA_descent(m, n) & 1) != par) mism++; }
            int b = 0; while (b < ncp && c > cps[b]) b++;
            if (b == ncp) continue;
            int ray = (n == 1) || (n == m - 1);
            Sb[b] += par ? -1 : 1; Nb[b] += 1; Rb[b] += (par ? -1 : 1) - (ray ? 1 : 0);
        }
    }
    printf("automaton vs literal descent on %ld pairs (m < 4000): mismatches %ld\n", checked, mism);
    printf("X N(X) S(X) S/sqrtX R(X) R/sqrtX   [(1/2+1/sqrt2) = %.5f]\n", 0.5 + 1 / sqrt(2.0));
    long long S = 0, N = 0, R = 0;
    for (int b = 0; b < ncp; b++) {
        S += Sb[b]; N += Nb[b]; R += Rb[b];
        printf("%.0e %lld %lld %.5f %lld %.5f\n", cps[b], N, S, S / sqrt(cps[b]), R, R / sqrt(cps[b]));
    }
    return 0;
}
