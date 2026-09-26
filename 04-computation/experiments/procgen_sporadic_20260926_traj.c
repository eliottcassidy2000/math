/*
 * procgen_sporadic_20260926_traj.c  (session collatz-procgen-20260922, lane "sporadic", 2026-09-26)
 *
 * Exhaustive least-element search for cycles of T(y) = y/2 (y even), (q*y + d)/2 (y odd)
 * on the positive integers, q odd >= 1, d odd (d may be negative as long as q*y + d > 0 for y >= 1).
 *
 * For every odd y in [ylo, yhi] (a residue class modulo 2*nstripe if nstripe > 1) the orbit of y is
 * followed on odd values only, x -> (q x + d) / 2^v, v = v_2(q x + d).  Since the intermediate even
 * values x' 2^j (j >= 1) are even and exceed the next odd value x', the orbit of y
 *   - drops below y          iff some odd value is < y          -> y is not the least element of a cycle;
 *   - returns to y           iff some odd value equals y        -> y is the least element of a cycle;
 * (an odd y can only equal an odd value).  Every cycle of positive integers has an odd least element
 * (an even least element x would have T(x) = x/2 < x in the cycle), so this finds EVERY cycle whose least
 * element lies in [ylo, yhi], with its period p (total steps) and a (odd steps).
 *
 * If neither happens within `steplim` odd steps, Brent's algorithm is run on the orbit (capped at
 * `brentlim` odd steps) and the cycle it lands on is walked once: least element m == y -> reported as a
 * cycle (long cycle); m > y -> reported as an "entered" cycle (it is found again, as a least element, if
 * m <= yhi); m < y -> the orbit drops below y later, nothing to report.  Brent cap reached -> UNRESOLVED.
 * Arithmetic is unsigned 128-bit; a value >= 2^120 is reported as OVERFLOW (never happens in the runs
 * recorded in the .out; the runner asserts that).
 *
 * Output lines (stdout):
 *   C <y> <p> <a> <maxodd>        least element y of a cycle, period p, odd steps a, largest odd value
 *   E <y> <m> <p> <a>             orbit of y entered a cycle with least element m > y (period p, a)
 *   U <y>                         unresolved (Brent cap reached)
 *   O <y>                         overflow
 *   S <nodd> <oddsteps>           summary: number of starts processed and total odd steps
 *
 * Usage: traj q d ylo yhi [nstripe stripe] [steplim] [brentlim]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

typedef unsigned __int128 u128;
typedef __int128 s128;

static inline int ctz128(u128 x) {
    uint64_t lo = (uint64_t)x;
    if (lo) return __builtin_ctzll(lo);
    return 64 + __builtin_ctzll((uint64_t)(x >> 64));
}

static void print_u128(u128 x) {
    char buf[64];
    int i = 63;
    buf[i] = 0;
    if (x == 0) { putchar('0'); return; }
    while (x) { buf[--i] = (char)('0' + (int)(x % 10)); x /= 10; }
    fputs(buf + i, stdout);
}

static const u128 LIMIT = ((u128)1) << 120;

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: traj q d ylo yhi [nstripe stripe] [steplim] [brentlim]\n");
        return 2;
    }
    int64_t q = atoll(argv[1]);
    int64_t d = atoll(argv[2]);
    uint64_t ylo = strtoull(argv[3], 0, 10), yhi = strtoull(argv[4], 0, 10);
    uint64_t nstripe = argc > 5 ? strtoull(argv[5], 0, 10) : 1;
    uint64_t stripe = argc > 6 ? strtoull(argv[6], 0, 10) : 0;
    uint64_t steplim = argc > 7 ? strtoull(argv[7], 0, 10) : 20000;
    uint64_t brentlim = argc > 8 ? strtoull(argv[8], 0, 10) : 50000000;
    if (q < 1 || (q % 2) == 0 || (d % 2) == 0) { fprintf(stderr, "q, d must be odd, q >= 1\n"); return 2; }
    if (ylo % 2 == 0) ylo++;
    const u128 Q = (u128)q;
    const s128 D = (s128)d;
    uint64_t nproc = 0;
    u128 totsteps = 0;
    for (uint64_t y = ylo; y <= yhi; y += 2) {
        if (nstripe > 1 && ((y >> 1) % nstripe) != stripe) continue;
        nproc++;
        u128 x = y;
        uint64_t k = 0, p = 0;
        u128 mx = x;
        int status = 0; /* 1 = cycle, 2 = below, 3 = overflow, 4 = limit */
        for (;;) {
            s128 t = (s128)(Q * x) + D;
            if (t <= 0) { status = 2; break; }       /* left the positive integers: not a positive cycle */
            u128 u = (u128)t;
            int v = ctz128(u);
            x = u >> v;
            k++;
            p += (uint64_t)v;
            if (x > mx) mx = x;
            if (x < y) { status = 2; break; }
            if (x == y) { status = 1; break; }
            if (x >= LIMIT) { status = 3; break; }
            if (k >= steplim) { status = 4; break; }
        }
        totsteps += k;
        if (status == 1) {
            printf("C %" PRIu64 " %" PRIu64 " %" PRIu64 " ", y, p, k);
            print_u128(mx);
            putchar('\n');
        } else if (status == 3) {
            printf("O %" PRIu64 "\n", y);
        } else if (status == 4) {
            /* Brent cycle detection on the odd-value map F(x) = (q x + d)/2^v */
            u128 power = 1, lam = 1, tort = x, hare;
            s128 t = (s128)(Q * x) + D;
            hare = ((u128)t) >> ctz128((u128)t);
            uint64_t cnt = 0;
            int ok = 1;
            while (tort != hare) {
                if (power == lam) { tort = hare; power <<= 1; lam = 0; }
                t = (s128)(Q * hare) + D;
                if (t <= 0 || hare >= LIMIT) { ok = 0; break; }
                hare = ((u128)t) >> ctz128((u128)t);
                lam++;
                if (++cnt > brentlim) { ok = 0; break; }
            }
            if (!ok) { printf("U %" PRIu64 "\n", y); fflush(stdout); continue; }
            /* hare is on the cycle; walk it once to get least element, p, a */
            u128 m = hare, z = hare;
            uint64_t pp = 0, aa = 0;
            do {
                t = (s128)(Q * z) + D;
                int v = ctz128((u128)t);
                z = ((u128)t) >> v;
                pp += (uint64_t)v; aa++;
                if (z < m) m = z;
            } while (z != hare);
            if (m > y) {
                printf("E %" PRIu64 " ", y);
                print_u128(m);
                printf(" %" PRIu64 " %" PRIu64 "\n", pp, aa);
            } else if (m == y) {
                /* y is the least element of a cycle longer than steplim odd steps */
                u128 mx2 = hare;
                z = hare;
                do {
                    t = (s128)(Q * z) + D;
                    z = ((u128)t) >> ctz128((u128)t);
                    if (z > mx2) mx2 = z;
                } while (z != hare);
                printf("C %" PRIu64 " %" PRIu64 " %" PRIu64 " ", y, pp, aa);
                print_u128(mx2);
                putchar('\n');
            }
            /* m < y: the orbit drops below y after the step limit, so y is not a least element */
            fflush(stdout);
        }
    }
    printf("S %" PRIu64 " ", nproc);
    print_u128(totsteps);
    putchar('\n');
    return 0;
}
