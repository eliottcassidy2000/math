#!/usr/bin/env python3
"""collatz_landing_20260926_probe2.py -- how does the landing multiplicity scale with L on long segments?
(session collatz-landing-20260926, opus; HYP-9161 probe)

For L in {20, 30, 40, 60, 80} and theta_X = 1.05 log_2(L)/L (the bootstrap's choice), X = 2^L: on long 3n+1
segments (orbits of 2^L - 1, of random odd starts near 2^L, of 2^(L-1)+1, and of 27*2^(L-5)+...), and on 5n+1
orbits (which keep climbing), tabulate: number of indices <= X, #ND, #D, number of landing points, the maximum
and mean multiplicity, and mean/(theta L) and max/L.
Usage: python3 collatz_landing_20260926_probe2.py
"""
import math, random


def T(n, q=3):
    return n // 2 if n % 2 == 0 else (q * n + 1) // 2


def orbit(n, steps, q=3):
    ys = [n]
    for _ in range(steps):
        n = T(n, q)
        ys.append(n)
        if n == 1:
            break
    return ys


def landing_stats(ys, L, theta):
    X = 2 ** L
    thr = 2.0 ** (-theta * L)
    mult = {}
    nd = 0
    total = 0
    for i, y in enumerate(ys):
        if y > X or i + L >= len(ys):
            continue
        total += 1
        land = None
        for s in range(1, L + 1):
            if ys[i + s] < y * thr:
                land = i + s
                break
        if land is None:
            nd += 1
        else:
            mult[land] = mult.get(land, 0) + 1
    ms = sorted(mult.values(), reverse=True)
    D = sum(ms)
    return total, nd, D, len(ms), (max(ms) if ms else 0), (D / len(ms) if ms else 0.0)


def main():
    random.seed(20260926)
    print("L  theta   thetaL | segment                         total    ND     D   landing  maxmult  mean  mean/(thetaL)  max/L")
    for L in (20, 30, 40, 60, 80):
        theta = 1.05 * math.log2(L) / L
        segs = [("3n+1 from 2^L-1", orbit(2 ** L - 1, 40 * L)),
                ("3n+1 from 2^(L-1)+1", orbit(2 ** (L - 1) + 1, 40 * L)),
                ("3n+1 from random odd ~2^L", orbit(random.randrange(2 ** (L - 1), 2 ** L) | 1, 40 * L)),
                ("3n+1 from random odd ~2^(L-3)", orbit(random.randrange(2 ** (L - 4), 2 ** (L - 3)) | 1, 40 * L)),
                ("5n+1 from 7 (climbs)", orbit(7, 40 * L, q=5)),
                ("5n+1 from random odd ~2^(L/2)", orbit(random.randrange(2 ** (L // 2 - 1), 2 ** (L // 2)) | 1, 40 * L, q=5))]
        for name, ys in segs:
            tot, nd, D, lp, mx, mean = landing_stats(ys, L, theta)
            print("%2d  %.3f  %5.2f | %-30s %6d %5d %5d %7d %7d  %5.2f  %8.2f  %6.3f" % (L, theta, theta * L, name, tot, nd, D, lp, mx, mean, mean / (theta * L), mx / L))


if __name__ == '__main__':
    main()
