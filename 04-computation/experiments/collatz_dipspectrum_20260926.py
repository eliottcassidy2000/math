#!/usr/bin/env python3
"""collatz_dipspectrum_20260926.py -- FINITE-EXACT control for the dip-spectrum theorem
(session collatz-exponent-atlas-20260926, opus; THM-4487).

Claim: for gamma in (log_4 3, 1] and b = +-1,
   D_b(X, gamma) := #{ n <= X : T_b^i(n) >= n^gamma for all 0 <= i <= floor(log_2 n) }
                  = X^(h(gamma / log_2 3) + o(1)),
where h is the binary entropy. At gamma = 1 this is the Terras undecided count X^(0.95);
at gamma = log_4 3 = 0.7925 (Korec's exponent) the predicted exponent reaches 1.
Every n <= X is iterated for floor(log_2 n) steps, the minimum ratio
min_i log T^i(n) / log n is recorded, and for a grid of gamma the exact counts D_b(2^t, gamma),
log D / log X, and the slopes over four doublings are printed. The finite-size counts carry a
polynomial prefactor (the prefix/ballot condition costs about 1/k at word length k), so the local
slope is expected near h(rho) - 1/(k ln 2), which is also printed.
(An earlier version of this script counted n < 2^(t+1) into the row for 2^t; fixed 2026-09-26.)
Usage: python3 collatz_dipspectrum_20260926.py [tmax]
"""
import math, sys

ALPHA = math.log2(3.0)


def h(p):
    return 0.0 if p <= 0 or p >= 1 else -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def pred(g):
    return h(max(0.5, g / ALPHA))


def min_exponent(n, b):
    if n < 2:
        return float('inf')
    k = n.bit_length() - 1
    x = n
    mn = n
    for _ in range(k):
        x = x // 2 if x % 2 == 0 else (3 * x + b) // 2
        if x < mn:
            mn = x
            if mn <= 1:
                break
    return math.log(mn) / math.log(n) if mn >= 1 else -1.0


def run(tmax, gammas):
    for b in (1, -1):
        print("== sheet 3n%+d, n <= 2^%d ==" % (b, tmax))
        cum = {t: [0] * len(gammas) for t in range(8, tmax + 1)}
        Xmax = 1 << tmax
        for n in range(2, Xmax + 1):
            e = min_exponent(n, b)
            # n <= 2^tt  iff  tt >= ceil(log2 n)
            tmin = n.bit_length() - 1 if (n & (n - 1)) == 0 else n.bit_length()
            for gi, g in enumerate(gammas):
                if e >= g:
                    for tt in range(max(tmin, 8), tmax + 1):
                        cum[tt][gi] += 1
        print("   gamma:        " + "  ".join("%7.3f" % g for g in gammas))
        print("   predicted:    " + "  ".join("%7.4f" % pred(g) for g in gammas))
        for t in range(10, tmax + 1, 2):
            print("   X=2^%2d D:      " % t + "  ".join("%7d" % c for c in cum[t]))
        for t in range(10, tmax + 1, 2):
            print("   X=2^%2d logD/logX:" % t + "  ".join("%7.4f" % (math.log(c) / math.log(1 << t) if c > 0 else float('nan')) for c in cum[t]))
        for t0 in range(10, tmax - 3, 2):
            t1 = t0 + 4
            print("   slope 2^%2d->2^%2d: " % (t0, t1) + "  ".join("%7.4f" % (math.log2(cum[t1][gi] / cum[t0][gi]) / 4 if cum[t0][gi] > 0 else float('nan')) for gi in range(len(gammas))))
        kmid = tmax - 2
        print("   predicted - 1/(k ln 2) at k=%d: " % kmid + "  ".join("%7.4f" % (pred(g) - 1 / (kmid * math.log(2))) for g in gammas))


if __name__ == '__main__':
    tmax = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    gammas = [0.7925, 0.82, 0.85, 0.88, 0.91, 0.94, 0.97, 1.0]
    run(tmax, gammas)
