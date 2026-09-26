/*
 * procgen_sporadic_20260926_sweep.c  (session collatz-procgen-20260922, lane "sporadic", 2026-09-26)
 *
 * Sweep of the (3m+d)-maps T_d(y) = y/2 (y even), (3y + d)/2 (y odd) on the positive integers over all
 * d in [dlo, dhi] with gcd(d, 6) = 1 (optionally one residue class d = stripe (mod nstripe) of them).
 * For each d every odd y in [ylo, X_d] (default [1, factor * d]) is tested as the least element of a cycle as in
 * procgen_sporadic_20260926_traj.c (odd-value orbit; stop when it drops below y or returns to y; Brent
 * after `steplim` odd steps).  This finds every cycle of T_d whose least element is <= X_d.
 * A cycle is PRIMITIVE (Lagarias 1990; Belaga-Mignotte 2006, Def. 3) iff gcd(m, d) = 1 for its least
 * element m (gcd(T_d(y), d) = gcd(y, d) along orbits because gcd(3, d) = 1).
 *
 * Output:
 *   c <d> <m> <p> <a> <prim>     one line per cycle (least element m, period p, odd steps a, prim 0/1)
 *   e <d> <y> <m> <p> <a>        orbit of y <= X_d entered a cycle with least element m > X_d (m < 2^64)
 *   D <d> <ylo> <yhi> <nprim> <nnonprim> <unresolved> <overflow>   (per d; cycles with least element in [ylo, yhi])
 * Usage: sweep dlo dhi factor [nstripe stripe] [steplim] [brentlim]      (ylo = 1, yhi = factor*d)
 *    or: sweep -f rangefile [nstripe stripe] [steplim] [brentlim]         (lines "d ylo yhi"; stripes by line)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

typedef unsigned __int128 u128;

static inline int ctz128(u128 x) {
    uint64_t lo = (uint64_t)x;
    if (lo) return __builtin_ctzll(lo);
    return 64 + __builtin_ctzll((uint64_t)(x >> 64));
}

static uint64_t gcd64(uint64_t a, uint64_t b) {
    while (b) { uint64_t t = a % b; a = b; b = t; }
    return a;
}

static const u128 LIMIT = ((u128)1) << 120;

static void run_one(uint64_t d, uint64_t ylo, uint64_t X, uint64_t steplim, uint64_t brentlim);

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: sweep dlo dhi factor [nstripe stripe] [steplim] [brentlim]\n"
                        "   or: sweep -f rangefile [nstripe stripe] [steplim] [brentlim]  (lines: d ylo yhi)\n");
        return 2;
    }
    if (argv[1][0] == '-' && argv[1][1] == 'f') {
        FILE *fp = fopen(argv[2], "r");
        if (!fp) { perror("rangefile"); return 2; }
        uint64_t nstripe = argc > 3 ? strtoull(argv[3], 0, 10) : 1;
        uint64_t stripe = argc > 4 ? strtoull(argv[4], 0, 10) : 0;
        uint64_t steplim = argc > 5 ? strtoull(argv[5], 0, 10) : 3000;
        uint64_t brentlim = argc > 6 ? strtoull(argv[6], 0, 10) : 50000000;
        uint64_t d, ylo, yhi, idx = 0;
        while (fscanf(fp, "%" SCNu64 " %" SCNu64 " %" SCNu64, &d, &ylo, &yhi) == 3) {
            if (nstripe > 1 && (idx++ % nstripe) != stripe) continue;
            run_one(d, ylo, yhi, steplim, brentlim);
        }
        fclose(fp);
        return 0;
    }
    if (argc < 4) { fprintf(stderr, "usage: sweep dlo dhi factor [nstripe stripe] [steplim] [brentlim]\n"); return 2; }
    uint64_t dlo = strtoull(argv[1], 0, 10), dhi = strtoull(argv[2], 0, 10);
    uint64_t factor = strtoull(argv[3], 0, 10);
    uint64_t nstripe = argc > 4 ? strtoull(argv[4], 0, 10) : 1;
    uint64_t stripe = argc > 5 ? strtoull(argv[5], 0, 10) : 0;
    uint64_t steplim = argc > 6 ? strtoull(argv[6], 0, 10) : 3000;
    uint64_t brentlim = argc > 7 ? strtoull(argv[7], 0, 10) : 50000000;
    for (uint64_t d = dlo; d <= dhi; d++) {
        if (d % 2 == 0 || d % 3 == 0) continue;
        if (nstripe > 1 && (d % nstripe) != stripe) continue;
        run_one(d, 1, factor * d, steplim, brentlim);
    }
    return 0;
}

static void run_one(uint64_t d, uint64_t ylo, uint64_t X, uint64_t steplim, uint64_t brentlim) {
    {
        uint64_t nprim = 0, nnon = 0, nunres = 0, nover = 0;
        const u128 D = (u128)d;
        if (ylo % 2 == 0) ylo++;
        for (uint64_t y = ylo; y <= X; y += 2) {
            u128 x = y;
            uint64_t k = 0, p = 0;
            int status = 0;
            for (;;) {
                u128 u = 3 * x + D;
                int v = ctz128(u);
                x = u >> v;
                k++;
                p += (uint64_t)v;
                if (x < y) { status = 2; break; }
                if (x == y) { status = 1; break; }
                if (x >= LIMIT) { status = 3; break; }
                if (k >= steplim) { status = 4; break; }
            }
            uint64_t cp = 0, ca = 0;
            if (status == 1) { cp = p; ca = k; }
            else if (status == 3) { nover++; printf("O %" PRIu64 " %" PRIu64 "\n", d, y); continue; }
            else if (status == 4) {
                u128 power = 1, lam = 1, tort = x, hare = (3 * x + D) >> ctz128(3 * x + D);
                uint64_t cnt = 0;
                int ok = 1;
                while (tort != hare) {
                    if (power == lam) { tort = hare; power <<= 1; lam = 0; }
                    if (hare >= LIMIT) { ok = 0; break; }
                    hare = (3 * hare + D) >> ctz128(3 * hare + D);
                    lam++;
                    if (++cnt > brentlim) { ok = 0; break; }
                }
                if (!ok) { nunres++; printf("U %" PRIu64 " %" PRIu64 "\n", d, y); continue; }
                u128 m = hare, z = hare;
                uint64_t pp = 0, aa = 0;
                do {
                    u128 u = 3 * z + D;
                    int v = ctz128(u);
                    z = u >> v;
                    pp += (uint64_t)v; aa++;
                    if (z < m) m = z;
                } while (z != hare);
                if (m > (u128)X) {            /* a cycle beyond the systematic range, met incidentally */
                    printf("e %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", d, y, (uint64_t)m, pp, aa);
                    continue;
                }
                if (m != (u128)y) continue;   /* m > y: found later as least element; m < y: not least */
                cp = pp; ca = aa;
            } else continue;
            int prim = gcd64(y, d) == 1;
            if (prim) nprim++; else nnon++;
            printf("c %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %d\n", d, y, cp, ca, prim);
        }
        printf("D %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", d, ylo, X, nprim, nnon, nunres, nover);
        fflush(stdout);
    }
}
