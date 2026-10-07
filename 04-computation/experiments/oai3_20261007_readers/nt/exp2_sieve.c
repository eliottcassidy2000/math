/*
 * exp2_sieve.c -- unconditional enumeration of negative discriminants D (fundamental
 * or not) whose form class group C(D) has exponent <= 2, for |D| = N <= X1.
 *
 * Criterion (Gauss): C(D) has exponent <= 2  <=>  every reduced form of discriminant D
 * (primitive or not) is ambiguous (b = 0, |b| = a, or a = c).  Non-primitive
 * non-ambiguous forms only occur when some C(D/g^2) -- a quotient of C(D) -- already has
 * exponent > 2, so marking them is harmless.
 *
 * Rejection lemma (EKN Lemma 1 / Weinberger Lemma 5 for c = 2, valid for orders):
 *   if an odd prime p has (D/p) = +1 and 4p^2 < N, then (p, b, (b^2+N)/(4p)) with
 *   0 < |b| < p, b == D mod 2, b^2 == D mod p, is a primitive reduced non-ambiguous form.
 *   If D == 1 mod 8 and N > 15, then (2, 1, (N+1)/8) is one.
 *
 * Part A: every N in [3, X0], N == 0,3 mod 4: prime test, then full reduced-form check.
 * Part B: N in (X0, X1]: CRT wheel mod M = 16*3*5*7*11*13*17*19*23 + 64-bit-word bitmask
 *         tables for primes 29..P2 (valid since 4*P2^2 < X0), then direct Legendre tests
 *         for primes up to sqrt(N)/2, then the full check for any survivor.
 *
 * Written for this lane (mac-mini-2026-10-07-oaimath3, lane nt); not openai/math code.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef unsigned long long u64;
typedef __uint128_t u128;

static int *primes; static int nprimes;

static void sieve_primes(int lim) {
    char *c = calloc(lim + 1, 1); primes = malloc(sizeof(int) * (lim / 2 + 10)); nprimes = 0;
    for (int i = 2; i <= lim; i++) { if (!c[i]) { primes[nprimes++] = i; for (long j = (long)i * i; j <= lim; j += i) c[j] = 1; } }
    free(c);
}

/* Jacobi symbol (a/n), n odd positive */
static int jacobi(u64 a, u64 n) {
    int t = 1; a %= n;
    while (a) {
        while ((a & 1) == 0) { a >>= 1; u64 r = n & 7; if (r == 3 || r == 5) t = -t; }
        u64 tmp = a; a = n; n = tmp;
        if ((a & 3) == 3 && (n & 3) == 3) t = -t;
        a %= n;
    }
    return n == 1 ? t : 0;
}

/* full check: return 1 if some reduced form (a,b,c) with 0<b<a<c has 4ac-b^2 = N */
static int has_nonambiguous(u64 N) {
    u64 amax = (u64)sqrtl((long double)N / 3.0L) + 2;
    for (u64 a = 2; a <= amax; a++) {
        u64 fa = 4 * a;
        for (u64 b = (N & 1) ? 1 : 2; b < a; b += 2) {
            u64 t = N + b * b;
            if (t % fa == 0) { u64 c = t / fa; if (c > a) return 1; }
        }
    }
    return 0;
}

/* prime pre-test: returns 1 if rejected by a split prime p (4p^2 < N) or by the 2-adic test */
static int split_reject(u64 N, int startidx) {
    if ((N & 7) == 7 && N > 15) return 1;
    for (int i = startidx; i < nprimes; i++) {
        u64 p = primes[i];
        if (4 * p * p >= N) break;
        u64 r = (p - N % p) % p;           /* r = D mod p */
        if (r && jacobi(r, p) == 1) return 1;
    }
    return 0;
}

static int cmpu64(const void *x, const void *y) { u64 a = *(const u64 *)x, b = *(const u64 *)y; return a < b ? -1 : a > b; }

int main(int argc, char **argv) {
    u64 X0 = argc > 1 ? strtoull(argv[1], 0, 10) : 100000000ULL;
    int KW = argc > 2 ? atoi(argv[2]) : 2;          /* words per bitmask: K = 64*KW */
    int P2 = argc > 3 ? atoi(argv[3]) : 211;
    sieve_primes(2000000);
    clock_t t0 = clock();

    /* ---------------- Part A ---------------- */
    u64 *found = malloc(sizeof(u64) * 100000); int nf = 0;
    u64 survA = 0;
    for (u64 N = 3; N <= X0; N++) {
        if ((N & 3) == 1 || (N & 3) == 2) continue;
        if (split_reject(N, 1)) continue;
        survA++;
        if (!has_nonambiguous(N)) found[nf++] = N;
    }
    fprintf(stderr, "Part A: N<=%llu: %llu prime-test survivors, %d with exponent<=2 (%.1fs)\n",
            X0, survA, nf, (double)(clock() - t0) / CLOCKS_PER_SEC);

    /* ---------------- Part B ---------------- */
    const int W[8] = {3, 5, 7, 11, 13, 17, 19, 23};
    u64 M = 16; for (int i = 0; i < 8; i++) M *= W[i];
    int K = 64 * KW;
    u64 X1 = (u64)K * M;
    if ((u64)4 * P2 * P2 >= X0) { fprintf(stderr, "P2 too large for X0\n"); return 1; }
    /* allowed residues */
    int L16[6] = {0, 3, 4, 8, 11, 12};
    int nL[8]; int *L[8];
    for (int i = 0; i < 8; i++) {
        int p = W[i]; L[i] = malloc(sizeof(int) * p); nL[i] = 0;
        for (int r = 0; r < p; r++) { int d = (p - r) % p; if (d == 0 || jacobi(d, p) == -1) L[i][nL[i]++] = r; }
    }
    /* CRT idempotents for moduli 16, W[0..7] */
    u64 mods[9]; mods[0] = 16; for (int i = 0; i < 8; i++) mods[i + 1] = W[i];
    u64 E[9];
    for (int i = 0; i < 9; i++) {
        u64 Mi = M / mods[i]; u64 inv = 0;
        for (u64 t = 1; t < mods[i]; t++) if ((Mi % mods[i]) * t % mods[i] == 1) { inv = t; break; }
        if (mods[i] == 1) inv = 0;
        E[i] = (u64)(((u128)Mi * inv) % M);
    }
    /* tier-2 primes */
    int t2[200]; int nt2 = 0;
    for (int i = 0; i < nprimes && primes[i] <= P2; i++) if (primes[i] > 23) t2[nt2++] = primes[i];
    /* bit tables */
    u64 **T = malloc(sizeof(u64 *) * nt2);
    for (int j = 0; j < nt2; j++) {
        int p = t2[j]; T[j] = calloc((size_t)p * KW, sizeof(u64));
        int mp = M % p;
        char *ok = malloc(p);
        for (int r = 0; r < p; r++) { int d = (p - r) % p; ok[r] = (d == 0 || jacobi(d, p) == -1); }
        for (int i = 0; i < p; i++)
            for (int k = 0; k < K; k++) if (ok[(i + (long)k * mp) % p]) T[j][(size_t)i * KW + k / 64] |= 1ULL << (k % 64);
        free(ok);
    }
    int idx3 = 0; while (primes[idx3] <= P2) idx3++;
    u64 nres = 0, surv2 = 0, surv3 = 0;
    u64 *mask = malloc(sizeof(u64) * KW);
    int c[9];
    /* nested CRT enumeration */
    for (c[0] = 0; c[0] < 6; c[0]++)
    for (c[1] = 0; c[1] < nL[0]; c[1]++)
    for (c[2] = 0; c[2] < nL[1]; c[2]++)
    for (c[3] = 0; c[3] < nL[2]; c[3]++)
    for (c[4] = 0; c[4] < nL[3]; c[4]++)
    for (c[5] = 0; c[5] < nL[4]; c[5]++)
    for (c[6] = 0; c[6] < nL[5]; c[6]++)
    for (c[7] = 0; c[7] < nL[6]; c[7]++)
    for (c[8] = 0; c[8] < nL[7]; c[8]++) {
        u128 acc = (u128)E[0] * L16[c[0]];
        for (int i = 1; i < 9; i++) acc += (u128)E[i] * L[i - 1][c[i]];
        u64 r = (u64)(acc % M);
        nres++;
        for (int w = 0; w < KW; w++) mask[w] = ~0ULL;
        int alive = 1;
        for (int j = 0; j < nt2 && alive; j++) {
            u64 *row = T[j] + (size_t)(r % t2[j]) * KW; alive = 0;
            for (int w = 0; w < KW; w++) { mask[w] &= row[w]; alive |= (mask[w] != 0); }
        }
        if (!alive) continue;
        for (int w = 0; w < KW; w++) {
            u64 m = mask[w];
            while (m) {
                int bit = __builtin_ctzll(m); m &= m - 1;
                u64 N = r + (u64)(w * 64 + bit) * M;
                if (N <= X0) continue;
                surv2++;
                if (split_reject(N, idx3)) continue;
                surv3++;
                fprintf(stderr, "tier-3 survivor N=%llu -- running full check\n", N);
                if (!has_nonambiguous(N)) { found[nf++] = N; fprintf(stderr, "  EXPONENT<=2: N=%llu\n", N); }
            }
        }
    }
    fprintf(stderr, "Part B: (%llu, %llu]: residues %llu, tier-2 survivors %llu, tier-3 survivors %llu (%.1fs)\n",
            X0, X1, nres, surv2, surv3, (double)(clock() - t0) / CLOCKS_PER_SEC);
    qsort(found, nf, sizeof(u64), cmpu64);
    printf("# negative discriminants -N with C(-N) of exponent <= 2, N <= %llu: %d values\n", X1, nf);
    for (int i = 0; i < nf; i++) printf("%llu\n", found[i]);
    return 0;
}
