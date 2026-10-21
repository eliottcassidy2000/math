/* AUDIT a6 (independent auditor, 2026-10-21): second, structurally different computation of s(m) mod p for a PRIME p
 * and m < p: ordinary power-series inversion of the EGF  u(x) = sum_m 2^C(m,2) x^m / m!  over F_p (factorials are
 * invertible for m < p), v = 1/u, and s(m) = -m! v_m.  No binomial coefficients are used.
 * Usage: ./a6_egf_prime m p      prints s(m) mod p, c(m) mod p and D(m) mod p.
 */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64;
typedef unsigned __int128 u128;
static u64 mulm(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }
static u64 powm(u64 b, u64 e, u64 m) { u64 r = 1 % m; b %= m; while (e) { if (e & 1) r = mulm(r, b, m); b = mulm(b, b, m); e >>= 1; } return r; }
int main(int argc, char **argv) {
    long N = atol(argv[1]); u64 p = strtoull(argv[2], 0, 10);
    if ((u64)N >= p) { fprintf(stderr, "need m < p\n"); return 1; }
    u64 *u = malloc((N + 1) * 8), *v = malloc((N + 1) * 8), *fact = malloc((N + 1) * 8);
    fact[0] = 1; for (long m = 1; m <= N; m++) fact[m] = mulm(fact[m - 1], m, p);
    u64 Tm = 1; /* 2^C(m,2) */
    for (long m = 0; m <= N; m++) { if (m > 0) Tm = mulm(Tm, powm(2, m - 1, p), p); u[m] = mulm(Tm, powm(fact[m], p - 2, p), p); }
    v[0] = 1;
    for (long m = 1; m <= N; m++) { u64 acc = 0; for (long k = 1; k <= m; k++) { acc += mulm(u[k], v[m - k], p); if (acc >= p) acc -= p; } v[m] = (p - acc) % p; }
    u64 sm = (p - mulm(fact[N], v[N], p)) % p;
    u64 cm = (N == 1) ? 1 : ((N & 1) ? powm(2, (N - 1) / 2, p) : 0);
    printf("m=%ld p=%llu s(m) mod p = %llu  c(m) mod p = %llu  D(m) mod p = %llu\n", N, p, sm, cm, (sm + p - cm) % p);
    return 0;
}
