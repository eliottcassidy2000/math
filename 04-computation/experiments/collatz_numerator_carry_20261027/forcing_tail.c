/* forcing_tail.c -- mac-mini-2026-10-27-potentials
 * HYP-9421 (4), pre-registered: max over primitive odd seeds s < N (s != 5 mod 8) of the total native price
 * F(s) = sum lambda along the compressed ROOT orbit of s down to 1, lambda = 2^a/(6X + 2^a), X = 3s+1.
 * Prediction: max F(s) < 0.1127503 for N = 10^8.  Long double accumulation; 128-bit integers for the orbit.
 * Build: cc -O2 -o forcing_tail forcing_tail.c ; run: ./forcing_tail 100000000 */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned __int128 u128;
int main(int argc, char **argv) {
    unsigned long long N = strtoull(argv[1], 0, 10);
    long double best = 0; unsigned long long arg = 0, maxsteps = 0;
    for (unsigned long long s0 = 3; s0 < N; s0 += 2) {
        if (s0 % 8 == 5) continue;
        u128 s = s0; long double F = 0; unsigned long long steps = 0;
        while (s != 1) {
            int a = (s % 8 == 1) ? 2 : 1;
            u128 X = 3 * s + 1;
            F += (long double)(1 << a) / (6.0L * (long double)X + (long double)(1 << a));
            u128 z = X >> a;
            while (z % 8 == 5) z = (z - 1) / 4;
            s = z; steps++;
            if (steps > 100000) { printf("orbit too long at %llu\n", s0); return 1; }
        }
        if (steps > maxsteps) maxsteps = steps;
        if (F > best) { best = F; arg = s0; }
    }
    printf("N = %llu: max F(s) = %.12Lf at s = %llu; longest native orbit %llu edges; prediction max F < 0.1127503: %s\n",
           N, best, arg, maxsteps, best < 0.1127503L ? "CONFIRMED" : "REFUTED");
    return 0;
}
