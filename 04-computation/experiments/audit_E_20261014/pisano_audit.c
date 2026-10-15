/* Audit E (2026-10-14): Pisano facts of the results note section 3, by direct iteration of (F_k, F_(k+1)) mod n
   (C, independent of pisano_facts.py).  Checks for n <= NMAX: max pi(n)/n and where pi(n) = 6n; pi(p) | p-1 for
   p = +-1 mod 5 and pi(p) | 2(p+1) for p = +-2 mod 5; pi(p^2) in {pi(p), p pi(p)} for primes p < 142;
   the maximum of pi(p^k)/p^k over prime powers p^k <= NMAX with p != 5; the owner's examples. */
#include <stdio.h>
#include <stdlib.h>
static long pis(long n) {
    if (n == 1) return 1;
    long a = 0, b = 1 % n, k = 0;
    do { long c = (a + b) % n; a = b; b = c; k++; } while (!(a == 0 && b == 1));
    return k;
}
static int isprime(long n) { if (n < 2) return 0; for (long d = 2; d * d <= n; d++) if (n % d == 0) return 0; return 1; }
int main(int argc, char **argv) {
    long NMAX = argc > 1 ? atol(argv[1]) : 20000;
    long *P = malloc(sizeof(long) * (NMAX + 1));
    double mx = 0; long ex = 0;
    for (long n = 2; n <= NMAX; n++) { P[n] = pis(n); if ((double)P[n] / n > mx) mx = (double)P[n] / n; if (P[n] > 6 * n) ex++; }
    printf("NMAX=%ld: max pi(n)/n = %.6f; #n with pi(n) > 6n: %ld; n with pi(n) = 6n:", NMAX, mx, ex);
    for (long n = 2; n <= NMAX; n++) if (P[n] == 6 * n) printf(" %ld", n);
    printf("\n");
    long bad = 0;
    for (long p = 7; p <= NMAX; p++) if (isprime(p)) {
        if ((p % 5 == 1 || p % 5 == 4) && (p - 1) % P[p]) bad++;
        if ((p % 5 == 2 || p % 5 == 3) && (2 * (p + 1)) % P[p]) bad++;
        if ((p % 5 == 2 || p % 5 == 3) && (p + 1) % P[p] == 0) bad++;   /* phi^(p+1) = -1: pi never divides p+1 */
    }
    printf("split pi | p-1, inert pi | 2(p+1) and pi does not divide p+1 (primes 7..NMAX): violations %ld\n", bad);
    long bad2 = 0;
    for (long p = 2; p * p <= NMAX; p++) if (isprime(p)) { long a = P[p], b = P[p * p]; if (b != a && b != p * a) bad2++; }
    printf("pi(p^2) in {pi(p), p pi(p)} for p^2 <= NMAX: violations %ld\n", bad2);
    double best = 0; long bp = 0, bk = 0;
    for (long p = 2; p <= NMAX; p++) if (isprime(p) && p != 5) { long q = p; int k = 1; while (q <= NMAX) {
        double r = (double)P[q] / q; if (r > best + 1e-12) { best = r; bp = p; bk = k; } q *= p; k++; } }
    printf("max pi(p^k)/p^k over prime powers <= NMAX with p != 5: %.6f at p=%ld (k=%ld); pi(5^k)/5^k:", best, bp, bk);
    for (long q = 5; q <= NMAX; q *= 5) printf(" %.3f", (double)P[q] / q);
    printf("\nexamples: pi(10)=%ld pi(11)=%ld pi(12)=%ld pi(19)=%ld pi(20)=%ld pi(50)=%ld pi(5)=%ld pi(2)=%ld pi(3)=%ld pi(7)=%ld\n",
           P[10], P[11], P[12], P[19], P[20], P[50], P[5], P[2], P[3], P[7]);
    printf("split primes < 100 with pi = p-1:"); for (long p = 7; p < 100; p++) if (isprime(p) && (p % 5 == 1 || p % 5 == 4) && P[p] == p - 1) printf(" %ld", p);
    printf("\ninert primes < 100 with pi = 2(p+1):"); for (long p = 7; p < 100; p++) if (isprime(p) && (p % 5 == 2 || p % 5 == 3) && P[p] == 2 * (p + 1)) printf(" %ld", p);
    printf("\n");
    return 0;
}
