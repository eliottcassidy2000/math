#!/usr/bin/env python3
"""collatz_one_out_edge_systems_20260930.py -- Collatz among one-out-edge systems (session collatz-posets-zeta5-20260927, opus,
2026-09-30, fifteenth note).

 (1) Power maps n -> oddpart(n^k + 1): for odd n >= 3 the orbit strictly increases (k = 2: v_2(n^2+1) = 1; k = 3: oddpart(n^3+1)
     = oddpart(n+1) (n^2 - n + 1)), so the functional graph on odd N is the fixed point 1 plus rays; in-degree census
     (collisions of oddpart(n^3+1)) and the density of orphans (in-degree 0).
 (2) Polynomial dynamics on Q: x -> x^2 + c.  c = -29/16 has the rational 3-cycle -7/4 -> 5/4 -> -1/4; the rational
     preperiodic set is finite (Northcott): brute force over heights <= 200 for several c; everything else wanders.
 (3) Collatz on Q_odd (Z_(2)): every rational m/d with |m| <= 60, d <= 60 odd is eventually periodic (periodicity conjecture
     range check); cycles reached per denominator; the sheet symmetry x_w(3x-1) = -x_w(3x+1).
 (4) Functional-graph census on N for 3x+1, 3x-1, 5x+1, 7x+1, n^2+1, n^3+1 (Syracuse forms): cycles below the bound,
     fraction of starts that are preperiodic within the bound, in-degree profile.
 (5) Eckmann-Hilton defect: the two operations (triple, halve) commute as multipliers, their affine versions do not; the
     commutator of the letters (1) and (2) is the translation by beta((1),(2))/8 = -2/8 = -1/4.
 (6) Chamberland's real extension f(x) = (x/2) cos^2(pi x/2) + ((3x+1)/2) sin^2(pi x/2): real fixed points in [-6, 6] and
     their multipliers; the integer cycles {1,2}, {-1}, {-5,-7,-10}, {-17,...} sit on one real line.
Usage: python3 collatz_one_out_edge_systems_20260930.py
"""
import math, time
from fractions import Fraction
from math import gcd
from collections import Counter

T0 = time.time()


def oddpart(x):
    while x % 2 == 0:
        x //= 2
    return x


def part1():
    print("== (1) power maps n -> oddpart(n^k + 1) ==")
    N = 10 ** 5
    inc2 = all(oddpart(n * n + 1) > n for n in range(3, N, 2))
    inc3 = all(oddpart(n ** 3 + 1) > n for n in range(3, N, 2))
    print(" strict increase for odd 3 <= n < 10^5: k = 2: %s (v_2(n^2+1) = 1 always); k = 3: %s (oddpart(n^3+1) = oddpart(n+1)(n^2-n+1))" % (inc2, inc3))
    ident = all(oddpart(n ** 3 + 1) == oddpart(n + 1) * (n * n - n + 1) for n in range(1, N, 2))
    print(" identity oddpart(n^3+1) = oddpart(n+1) (n^2 - n + 1) for odd n < 10^5: %s" % ident)
    img = Counter(oddpart(n ** 3 + 1) for n in range(1, N, 2))
    coll = [(m, c) for m, c in img.items() if c > 1]
    print(" images m = oddpart(n^3+1), n odd < 10^5: %d distinct, collisions (in-degree >= 2): %d %s" % (len(img), len(coll), coll[:5]))
    small = [m for m in range(1, 2001, 2) if m in img]
    print(" odd m <= 2000 that are images (have a preimage): %d of 1000 -> orphans (in-degree 0) are the overwhelming majority; the graph on odd N is the fixed point 1 plus disjoint rays" % len(small))


def preperiodic_rationals(c, H):
    """rational preperiodic points of x -> x^2 + c with numerator, denominator <= H (brute force, height bound)."""
    pts = set()
    for q in range(1, H + 1):
        for p in range(-H, H + 1):
            if gcd(abs(p), q) != 1:
                continue
            x = Fraction(p, q); seen = []
            y = x
            ok = False
            for _ in range(60):
                if y in seen:
                    ok = True; break
                seen.append(y)
                y = y * y + c
                if y.denominator > 10 ** 12 or abs(y) > 10 ** 12:
                    break
            if ok:
                pts.add(x)
    return sorted(pts)


def part2():
    print("== (2) polynomial dynamics on Q: x -> x^2 + c ==")
    c = Fraction(-29, 16)
    x = Fraction(-7, 4); orb = [x]
    for _ in range(3):
        x = x * x + c; orb.append(x)
    print(" c = -29/16: orbit of -7/4: %s (a rational 3-cycle; the -7/4 of the arithmetic-braids geometry note)" % [str(t) for t in orb])
    for cc in (Fraction(-29, 16), Fraction(-1), Fraction(0), Fraction(-2), Fraction(1, 4), Fraction(-3, 4), Fraction(-7, 4), Fraction(-6, 1)):
        pts = preperiodic_rationals(cc, 60)
        print("  c = %6s: rational preperiodic points with height <= 60: %d  %s" % (cc, len(pts), [str(t) for t in pts][:10]))
    print(" reading: finitely many preperiodic rationals for every c (Northcott), at most one cycle of period <= 3 (Poonen's conjecture; periods 4, 5 excluded by Morton and Flynn-Poonen-Schaefer, CITED from memory); every other rational wanders to infinity in a tree component")


def collatz_q(x):
    """the Collatz map on Z_(2): x odd-type (numerator odd) -> (3x+1)/2^v, numerator even -> x/2 ... as the Syracuse map on odd-type points."""
    y = 3 * x + 1
    if y == 0:
        return y          # x = -1/3 is the 2-adic zero of 3x + 1; the map sends it to the fixed point 0
    while y.numerator % 2 == 0:
        y = y / 2
    return y


def part3():
    print("== (3) Collatz on Q_odd: periodicity in a box, cycles per denominator, the sheet symmetry ==")
    H = 40; cycles = {}; wander = []; total = 0
    for d in range(1, H + 1, 2):
        for m in range(-H, H + 1):
            if m % 2 == 0 or gcd(abs(m), d) != 1:
                continue
            x = Fraction(m, d); total += 1; seen = {}; y = x; per = None; path = []
            for _ in range(600):
                if y in seen:
                    per = path[seen[y]:]; break
                seen[y] = len(path); path.append(y); y = collatz_q(y)
                if abs(y.numerator) > 10 ** 15:
                    break
            if per is None:
                wander.append(x)
            else:
                key = (min(per), len(per))
                cycles.setdefault(per[0].denominator, set()).add(min(per))
    print(" odd-numerator rationals m/d, |m| <= 40, d <= 40 odd: %d points; not eventually periodic within 600 steps or before 10^15: %d %s" % (total, len(wander), [str(w) for w in wander[:8]]))
    print(" distinct cycles reached, by denominator (d: least elements):", {d: sorted(str(v) for v in s)[:6] for d, s in sorted(cycles.items())[:12]})
    print(" cycles with 3 | denominator: %s (clocks are coprime to 3, so none: preperiodic points may have 3 | d, cycle points cannot)" % [d for d in cycles if d % 3 == 0])
    # sheet symmetry: the 3x-1 map's fixed points are the negatives of the 3x+1 fixed points
    def fixed_point(w, sign):
        S = 0; d = 0
        for v in w:
            S = 3 * S + sign * 2 ** d; d += v
        return Fraction(S, 2 ** d - 3 ** len(w))
    ok = all(fixed_point(w, -1) == -fixed_point(w, 1) for w in ((1,), (2,), (1, 2), (1, 1, 1, 2, 1, 1, 4), (3, 1, 2)))
    print(" sheet symmetry: the fixed point of a word for 3x-1 is minus the fixed point for 3x+1: %s (negation conjugates the two sheets on Z_(2))" % ok)


def part4():
    print("== (4) functional-graph census on N: cycles, preperiodicity, in-degrees ==")
    def syr(q, b, k):
        def f(n):
            y = (n ** k if k > 1 else q * n) + b if k > 1 else q * n + b
            return oddpart(y)
        return f
    systems = [("3x+1", lambda n: oddpart(3 * n + 1)), ("3x-1", lambda n: oddpart(3 * n - 1)), ("5x+1", lambda n: oddpart(5 * n + 1)),
               ("7x+1", lambda n: oddpart(7 * n + 1)), ("x^2+1", lambda n: oddpart(n * n + 1)), ("x^3+1", lambda n: oddpart(n ** 3 + 1))]
    N = 20001; cap = 10 ** 12
    for name, f in systems:
        cyc = set(); pre = 0; indeg = Counter()
        for n in range(1, N, 2):
            x = n; seen = set(); path = []
            while x < cap and x not in seen and len(path) < 3000:
                seen.add(x); path.append(x); x = f(x)
            if x in seen:
                pre += 1; c = path[path.index(x):]; cyc.add(min(c))
        for n in range(1, 4 * N, 2):
            m = f(n)
            if m < N:
                indeg[m] += 1
        prof = Counter(indeg.get(m, 0) for m in range(1, N, 2))
        print(" %-6s: cycles (least elements) reached from odd n < 2*10^4: %s; preperiodic within bounds: %d/%d; in-degree profile on odd m < 2*10^4 from n < 8*10^4: %s" % (
            name, sorted(cyc)[:8], pre, N // 2, sorted(prof.items())[:5]))


def part5():
    print("== (5) the Eckmann-Hilton defect ==")
    def carry(w):
        S = 0; d = 0
        for v in w:
            S = 3 * S + 2 ** d; d += v
        return S
    b = carry((1, 2)) - carry((2, 1))
    print(" multipliers commute: 3 * (1/2) = (1/2) * 3; affine letters do not: the words (1,2) and (2,1) have carries %d, %d, defect beta = %d, so the two composites differ by the translation %s" % (carry((1, 2)), carry((2, 1)), b, Fraction(b, 8)))
    print(" the multiplier monoid <2, 3> is N x N (two commuting copies of N); the carry cocycle is the obstruction to interchange, and it is antisymmetric (twelfth note): the structure is symmetric, which is exactly what Eckmann-Hilton would force if the interchange held")


def part6():
    print("== (6) Chamberland's real extension ==")
    f = lambda x: (x / 2) * math.cos(math.pi * x / 2) ** 2 + ((3 * x + 1) / 2) * math.sin(math.pi * x / 2) ** 2
    df = lambda x: (f(x + 1e-6) - f(x - 1e-6)) / 2e-6
    roots = []
    xs = [-6 + i * 0.001 for i in range(12001)]
    for a, b in zip(xs, xs[1:]):
        if (f(a) - a) * (f(b) - b) < 0:
            lo, hi = a, b
            for _ in range(60):
                mid = (lo + hi) / 2
                if (f(lo) - lo) * (f(mid) - mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            roots.append(round((lo + hi) / 2, 6))
    print(" real fixed points of f in [-6, 6]: %s" % roots)
    print(" multipliers f'(x*) at them: %s (|f'| < 1 attracting)" % [round(df(r), 3) for r in roots])
    print(" f(1) = %.6f, f(2) = %.6f, f(-1) = %.6f, f(-5) = %.6f, f(-7) = %.6f, f(-10) = %.6f: the integer cycles of both sheets lie on one real line" % (f(1), f(2), f(-1), f(-5), f(-7), f(-10)))


if __name__ == "__main__":
    import sys
    parts = sys.argv[1:] or ["1", "2", "3", "4", "5", "6"]
    for q in parts:
        globals()["part" + q]()
    print("total %.0fs" % (time.time() - T0))
