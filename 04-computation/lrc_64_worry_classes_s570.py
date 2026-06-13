#!/usr/bin/env python3
"""Continue the worry framework with the 64 self-converse round classes (n=14).
opus-2026-06-03-S570. Census the transversals (realizing speed sets of the worry-
container): all lonely via an unblocked small pair (S569); AP the unique floor-tight."""
from fractions import Fraction
from math import gcd
from itertools import combinations, product
N = 14; THR = Fraction(1, N); C = 2 * N - 1

def dist(x):
    x = x % 1
    return min(x, 1 - x)

def positive_measure(V):
    eps = set()
    for v in V:
        for k in range(v + 1):
            for s in (-1, 1):
                eps.add(Fraction(k * N + s, N * v) % 1)
    pts = sorted(eps); L = len(pts)
    for i in range(L):
        a, b = pts[i], pts[(i + 1) % L]; ln = (b - a) if b > a else (b - a + 1)
        if all(dist(v * ((a + ln / 2) % 1)) > THR for v in V):
            return True
    return False

def floor_tight(V):
    return any(all(dist(v * Fraction(j, N)) >= THR for v in V) for j in range(1, N))

def unblocked_small(V):
    for a, b in [(x, y) for x, y in combinations(sorted(V), 2) if (x + y) // gcd(x, y) <= N]:
        D = a + b
        for m in range(1, D):
            if dist(Fraction(a * m, D)) >= THR and dist(Fraction(b * m, D)) >= THR:
                if all(dist(Fraction(c * m, D)) >= THR for c in V if c not in (a, b)):
                    return True
    return False

if __name__ == "__main__":
    pairs = [(a, C - a) for a in range(1, N)]
    tot = lonely = 0; tight = []; no_ub = []
    for choice in product(*[(lo, hi) for lo, hi in pairs]):
        V = tuple(sorted(choice)); g = 0
        for v in V:
            g = gcd(g, v)
        if g != 1:
            continue
        tot += 1
        pm = positive_measure(V); ft = floor_tight(V)
        if pm or ft:
            lonely += 1
        if (not pm) and ft:
            tight.append(V)
        if not unblocked_small(V):
            no_ub.append(V)
    print(f"transversals (gcd-1) n=14: {tot}; lonely={lonely}; floor-tight(worry)={tight}; "
          f"no-unblocked-small-pair={len(no_ub)}")
