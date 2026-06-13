#!/usr/bin/env python3
"""Proof-program check: across even n, is every measure-0 (worry-set) config exactly
FLOOR-TIGHT (M=1/n), never below, and a perfect antipodal transversal mod 2n-1?
opus-2026-06-02-S568. LRC(n) <=> every measure-0 set is floor-tight."""
from fractions import Fraction
from math import gcd
from itertools import combinations

def positive_measure(V, n):
    thr = Fraction(1, n); eps = set()
    for v in V:
        for k in range(v + 1):
            for s in (-1, 1):
                eps.add(Fraction(k * n + s, n * v) % 1)
    pts = sorted(eps); L = len(pts)
    for i in range(L):
        a, b = pts[i], pts[(i + 1) % L]; ln = (b - a) if b > a else (b - a + 1)
        mid = (a + ln / 2) % 1
        if all(min((v * mid) % 1, 1 - (v * mid) % 1) > thr for v in V):
            return True
    return False

def nonempty_floor(V, n):
    thr = Fraction(1, n); eps = set()
    for v in V:
        for k in range(v + 1):
            for s in (-1, 1):
                eps.add(Fraction(k * n + s, n * v) % 1)
    return any(all(min((v * t) % 1, 1 - (v * t) % 1) >= thr for v in V) for t in sorted(eps))

def is_transversal(V, n):
    M = 2 * n - 1; res = [v % M for v in V]
    if 0 in res:
        return False
    return all(sum(1 for r in res if r in (a, M - a)) == 1 for a in range(1, n))

def prim(V):
    g = 0
    for v in V:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in V))

if __name__ == "__main__":
    boxes = {4: 9, 6: 11, 8: 13, 10: 15, 12: 17, 14: 18}
    for n in [4, 6, 8, 10, 12, 14]:
        m = n - 1; B = boxes[n]
        tested = meas0 = cex = transv = 0; spor = []
        for combo in combinations(range(1, B + 1), m):
            if prim(combo) != tuple(combo):
                continue
            tested += 1
            if positive_measure(combo, n):
                continue
            meas0 += 1
            if not nonempty_floor(combo, n):
                cex += 1
            if is_transversal(combo, n):
                transv += 1
            elif nonempty_floor(combo, n):
                spor.append(combo)
        print(f"n={n:2d} box[1,{B}] ({tested} sets): measure-0={meas0}, "
              f"COUNTEREXAMPLES={cex}, transversal={transv}/{meas0}"
              + (f", non-transv sporadic={len(spor)} e.g.{spor[0]}" if spor else ""))
    print("\nLRC(n) <=> every measure-0 set is floor-tight (counterexamples must be 0).")
