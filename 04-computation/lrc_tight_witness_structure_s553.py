#!/usr/bin/env python3
"""
Structure of the n=14 tight family: where do the measure-zero witnesses sit?
opus-2026-06-01-S553 (remote-control).

For each TIGHT n=14 set (M=1/14, safe set = finitely many points) we list the
exact witness times and test the hypothesis that every tight witness is a
half-division point t = j/(2n).  A clean "yes" gives a canonical witness lattice
for the whole tight family -- the object a proof of "every tight config is
lonely" would use.
"""

from fractions import Fraction
from math import gcd
from itertools import combinations


def primitive(V):
    g = 0
    for v in V:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in V))


def witnesses(V, n):
    """All measure-zero safe points of a tight set (exact)."""
    thr = Fraction(1, n)
    endpoints = set()
    for v in V:
        for k in range(v + 1):
            for s in (-1, 1):
                endpoints.add(Fraction(k * n + s, v * n) % 1)
    out = []
    for t in sorted(endpoints):
        if all(min((v * t) % 1, 1 - (v * t) % 1) >= thr for v in V):
            out.append(t)
    return out


def is_tight(V, n):
    """measure 0 and nonempty."""
    thr = Fraction(1, n)
    endpoints = set()
    for v in V:
        for k in range(v + 1):
            for s in (-1, 1):
                endpoints.add(Fraction(k * n + s, v * n) % 1)
    pts = sorted(endpoints)
    m = len(pts)
    for i in range(m):
        a = pts[i]
        b = pts[(i + 1) % m]
        length = (b - a) if b > a else (b - a + 1)
        mid = (a + length / 2) % 1
        if all(min((v * mid) % 1, 1 - (v * mid) % 1) >= thr for v in V):
            return False                      # positive measure
    return len(witnesses(V, n)) > 0


def run(n=14, B=18):
    m = n - 1
    tight = []
    for combo in combinations(range(1, B + 1), m):
        if primitive(combo) != tuple(combo):
            continue
        if is_tight(combo, n):
            tight.append(combo)
    print(f"== n={n}: {len(tight)} tight sets among primitive {m}-subsets of "
          f"[1,{B}] ==\n")
    all_half_div = True
    for V in tight:
        ws = witnesses(V, n)
        half_div = all((w * 2 * n).denominator == 1 for w in ws)
        all_half_div = all_half_div and half_div
        tag = "AP" if primitive(V) == tuple(range(1, m + 1)) else "non-AP"
        wlist = ", ".join(f"{w.numerator}/{w.denominator}" for w in ws)
        print(f"  {tag:6s} {V}")
        print(f"         witnesses ({len(ws)}): {wlist}")
        print(f"         all at j/(2n)={2*n}-division points? {half_div}")
    print(f"\n  => EVERY tight witness is a j/(2n) division point: {all_half_div}")


if __name__ == "__main__":
    run(14, 18)
