/* Audit B (2026-10-14): independent Wilson-quotient scan for THM-4651 (1).
 *
 * Method (different from the session's PARI scan, which multiplies 1..p-1 modulo p^2):
 *   (p-1)! = prod_{l prime < p} l^{e_l},  e_l = sum_i floor((p-1)/l^i)  (Legendre),
 * each prime power taken modulo p^2 by binary exponentiation with unsigned __int128 products.
 * Then W = ((p-1)! + 1)/p mod p  (exact: (p-1)! == -1 mod p by Wilson, asserted).
 *
 * Second, for p < PDIRECT, the Giuga sum S(p) = sum_{k<p} k^(p-1) mod p^2 is computed directly
 * (k^(p-1) by binary exponentiation) and compared with p - 1 + p*w_p and with p + (p-1)!.
 *
 * Usage: wilson_legendre N PDIRECT
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;

static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)(((u128)a * b) % m); }
static u64 powmod(u64 b, u64 e, u64 m) {
    u64 r = 1 % m; b %= m;
    while (e) { if (e & 1) r = mulmod(r, b, m); b = mulmod(b, b, m); e >>= 1; }
    return r;
}

int main(int argc, char **argv) {
    u64 N = argc > 1 ? strtoull(argv[1], 0, 10) : 200000;
    u64 PD = argc > 2 ? strtoull(argv[2], 0, 10) : 20000;
    char *comp = calloc(N + 1, 1);
    u64 *pr = malloc(sizeof(u64) * (N + 1));
    u64 np = 0;
    for (u64 i = 2; i <= N; i++) {
        if (!comp[i]) { pr[np++] = i; for (u64 j = i * i; j <= N; j += i) comp[j] = 1; }
    }
    u64 cnt0 = 0, cntm1 = 0, cntp1 = 0, nchecked = 0, ndirect = 0;
    long long hist[5] = {0, 0, 0, 0, 0}; /* w = -2,-1,0,1,2 */
    printf("primes p with 3 <= p < %llu\n", (unsigned long long)N);
    for (u64 ip = 1; ip < np && pr[ip] < N; ip++) {
        u64 p = pr[ip], p2 = p * p, f = 1;
        for (u64 il = 0; il < ip; il++) {
            u64 l = pr[il], e = 0, q = p - 1;
            while (q) { q /= l; e += q; }
            f = mulmod(f, powmod(l, e, p2), p2);
        }
        if ((f + 1) % p != 0) { printf("Wilson theorem FAILED at %llu\n", (unsigned long long)p); return 1; }
        u64 w = ((f + 1) / p) % p;
        nchecked++;
        if (w == 0) { printf("  w_p == 0  : %llu\n", (unsigned long long)p); cnt0++; }
        if (w == p - 1) { printf("  w_p == -1 : %llu\n", (unsigned long long)p); cntm1++; }
        if (w == 1) { printf("  w_p == +1 : %llu\n", (unsigned long long)p); cntp1++; }
        if (w == p - 2) { printf("  w_p == -2 : %llu\n", (unsigned long long)p); }
        if (w == 2) { printf("  w_p == +2 : %llu\n", (unsigned long long)p); }
        if (p < PD) {
            u64 S = 0;
            for (u64 k = 1; k < p; k++) { S += powmod(k, p - 1, p2); if (S >= p2) S -= p2; }
            u64 rhs1 = (p - 1 + mulmod(p, w, p2)) % p2;
            u64 rhs2 = (p + f) % p2;
            if (S != rhs1 || S != rhs2) { printf("congruence FAILED at %llu\n", (unsigned long long)p); return 1; }
            ndirect++;
        }
    }
    printf("primes scanned: %llu; w==0: %llu, w==-1: %llu, w==+1: %llu\n", (unsigned long long)nchecked,
           (unsigned long long)cnt0, (unsigned long long)cntm1, (unsigned long long)cntp1);
    printf("direct check S(p) == p-1+p*w_p == p+(p-1)! (mod p^2) for all %llu primes 3 <= p < %llu: OK\n",
           (unsigned long long)ndirect, (unsigned long long)PD);
    (void)hist;
    return 0;
}
