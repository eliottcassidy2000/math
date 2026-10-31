/* opus-2026-10-31-S70: Conway's PRIMEGAME and the owner's variant (15/14 -> 15/2, 1/7), exponent-vector simulation.

   Registers = exponents of the first ten primes 2,3,5,7,11,13,17,19,23,29.
   Checks, for both programs, up to a bound on n:
     (1) the powers of 2 reached are 2^p for exactly the primes p, in order (standard program);
     (2) phase law: every entry into a configuration with no state prime (11..29) is 2^n 7^(D(n)-1),
         D(n) = n / (least prime factor of n) the largest proper divisor (D(prime) = 1, D(1) = 1);
     (3) the variant reaches exactly the same configurations at every phase entry and the same powers of 2.
   usage: primegame_variant_phase_law NMAX
*/
#include <stdio.h>
#include <stdlib.h>

static const int P[10] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
typedef struct { int num[10], den[10]; } Frac;

static void factor(int x, int v[10]) {
    for (int i = 0; i < 10; i++) { v[i] = 0; while (x % P[i] == 0) { x /= P[i]; v[i]++; } }
    if (x != 1) { fprintf(stderr, "bad factor\n"); exit(1); }
}
static void mk(Frac *f, int a, int b) { factor(a, f->num); factor(b, f->den); }

static int spf(int n) { for (int p = 2; p * p <= n; p++) if (n % p == 0) return p; return n; }

static int nfr;
static Frac STD[14], OWN[14];

/* run until the configuration 2^nmax (or beyond) is reached; record phase entries and powers of 2 */
static int run(Frac *pr, long long nmax, int check_phase, long long *phases_out, long long *pow2_out, int *npow2) {
    int x[10] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int inphase = 0; long long nphase = 0; *npow2 = 0;
    for (long long t = 0;; t++) {
        int nostate = 1;
        for (int i = 4; i < 10; i++) if (x[i]) { nostate = 0; break; }
        if (nostate && !inphase) {
            nphase++;
            /* phase entry: expect 2^n 7^(D(n)-1) with r3 = r5 = 0 */
            int n = x[0];
            if (check_phase) {
                int D = (n == 1) ? 1 : n / spf(n);
                if (x[1] != 0 || x[2] != 0 || x[3] != D - 1) {
                    printf("PHASE LAW FAILS at phase %lld: r2=%d r3=%d r5=%d r7=%d (expected r7=%d)\n", nphase, x[0], x[1], x[2], x[3], D - 1);
                    return 0;
                }
            }
            if (x[3] > x[0]) { printf("INVARIANT r7 <= r2 FAILS at phase %lld\n", nphase); return 0; }
            if (n >= nmax) { *phases_out = nphase; return 1; }
            inphase = 1;
        }
        if (!nostate) inphase = 0;
        int k;
        for (k = 0; k < 14; k++) {
            int ok = 1;
            for (int i = 0; i < 10; i++) if (x[i] < pr[k].den[i]) { ok = 0; break; }
            if (ok) break;
        }
        if (k == 14) { printf("HALT\n"); return 0; }
        for (int i = 0; i < 10; i++) x[i] += pr[k].num[i] - pr[k].den[i];
        int pure2 = 1;
        for (int i = 1; i < 10; i++) if (x[i]) { pure2 = 0; break; }
        if (pure2) pow2_out[(*npow2)++] = x[0];
    }
}

int main(int argc, char **argv) {
    long long nmax = argc > 1 ? atoll(argv[1]) : 1000;
    int a[14] = {17, 78, 19, 23, 29, 77, 95, 77, 1, 11, 13, 15, 15, 55};
    int b[14] = {91, 85, 51, 38, 33, 29, 23, 19, 17, 13, 11, 14, 2, 1};
    int c[14] = {17, 78, 19, 23, 29, 77, 95, 77, 1, 11, 13, 15, 1, 55};
    int d[14] = {91, 85, 51, 38, 33, 29, 23, 19, 17, 13, 11, 2, 7, 1};
    for (int i = 0; i < 14; i++) { mk(&STD[i], a[i], b[i]); mk(&OWN[i], c[i], d[i]); }
    static long long p2s[200000], p2o[200000]; int ns, no; long long phs, pho;
    int ok1 = run(STD, nmax, 1, &phs, p2s, &ns);
    int ok2 = run(OWN, nmax, 1, &pho, p2o, &no);
    /* primes check on the standard program */
    int primes_ok = 1, cnt = 0;
    for (int i = 0; i < ns; i++) {
        int p = (int)p2s[i];
        if (spf(p) != p || p < 2) primes_ok = 0;
        cnt++;
    }
    /* every prime below nmax appears */
    int k = 0;
    for (int p = 2; p < nmax; p++) if (spf(p) == p) { if (k >= ns || p2s[k] != p) primes_ok = 0; k++; }
    int same = (ns == no);
    for (int i = 0; same && i < ns; i++) if (p2s[i] != p2o[i]) same = 0;
    printf("NMAX=%lld standard_ok=%d variant_ok=%d phases=%lld/%lld powers_of_2=%d/%d primes_exact=%d same_powers=%d\n",
           nmax, ok1, ok2, phs, pho, ns, no, primes_ok, same);
    return (ok1 && ok2 && primes_ok && same) ? 0 : 1;
}
