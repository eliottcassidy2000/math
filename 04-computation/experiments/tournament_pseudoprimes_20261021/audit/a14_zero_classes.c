/* AUDIT a14 (independent auditor, 2026-10-21): zero classes Z_p = {1 <= r <= p-1 : p | D(r)} for odd primes p <= PMAX
 * via the EGF-inversion method (a6-style; factorials invertible since r < p).  Prints "p r1 r2 ..." per prime. */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64;
static u64 mulm(u64 a, u64 b, u64 m) { return a * b % m; }
static u64 powm(u64 b, u64 e, u64 m) { u64 r = 1 % m; b %= m; while (e) { if (e & 1) r = mulm(r, b, m); b = mulm(b, b, m); e >>= 1; } return r; }
static int isprime(int n) { if (n < 2) return 0; for (int d = 2; d * d <= n; d++) if (n % d == 0) return 0; return 1; }
int main(int argc, char **argv) {
    int PMAX = atoi(argv[1]);
    for (int p = 3; p <= PMAX; p += 2) {
        if (!isprime(p)) continue;
        int N = p - 1; u64 *u = malloc((N + 1) * 8), *v = malloc((N + 1) * 8), *f = malloc((N + 1) * 8);
        f[0] = 1; for (int m = 1; m <= N; m++) f[m] = mulm(f[m - 1], m, p);
        u64 T = 1; for (int m = 0; m <= N; m++) { if (m) T = mulm(T, powm(2, m - 1, p), p); u[m] = mulm(T, powm(f[m], p - 2, p), p); }
        v[0] = 1; for (int m = 1; m <= N; m++) { u64 a = 0; for (int k = 1; k <= m; k++) a = (a + mulm(u[k], v[m - k], p)) % p; v[m] = (p - a) % p; }
        printf("%d", p);
        for (int r = 1; r <= N; r++) { u64 s = (p - mulm(f[r], v[r], p)) % p; u64 c = (r == 1) ? 1 : ((r & 1) ? powm(2, (r - 1) / 2, p) : 0);
            if (s == c) printf(" %d", r); }
        printf("\n"); free(u); free(v); free(f);
    }
    return 0; }
