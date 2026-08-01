#!/usr/bin/env python3
"""The Factorial Conjecture FC(2) and FC(3): moment-determinacy framing + census.

FC(n): for f in C[x_1..x_n] with L(x^alpha) = alpha!, if L(f^m) = 0 for every
m >= 1 then f = 0.

THE FRAMING.  L is the moment functional of the PRODUCT EXPONENTIAL measure:

    L(g) = int_{[0,inf)^n} g(x) e^{-(x_1+...+x_n)} dx,     L(x^alpha) = alpha!.

So L(f^m) is the m-th moment of the pushforward nu = f_* mu.  FC says: all
moments (m >= 1) of nu vanish => f = 0.  If nu were DETERMINATE this is
immediate (vanishing moments + total mass force nu = delta_0, so f = 0 a.e.,
so f = 0 as a polynomial).  Hence FC is exactly a moment-determinacy question,
and the classical Krein/Carleman criterion locates the difficulty by DEGREE,
not by dimension: a density decaying like exp(-t^rho) is determinate iff
rho >= 1/2, and the pushforward of e^{-|x|} under a degree-d polynomial has
tail exp(-t^{1/d}).  So

    deg f <= 2  ->  rho = 1/d >= 1/2  ->  determinate  ->  the argument closes
    deg f >= 3  ->  rho < 1/2         ->  INDETERMINATE -> the argument fails

which is why FC is easy in low degree and hard from cubics on, in EVERY
number of variables.  (Caveat: for complex f the pushforward is a complex
measure and the criterion needs the real/imaginary decomposition; the clean
statement above is for real f.)

THE CENSUS.  Independently, search small-support integer f for the first few
vanishing moments, which is the finite shadow of the conjecture.
"""
import sys
from itertools import product
from math import factorial


def L_of_power(coeffs, m, nvars):
    """L(f^m) exactly, f = sum coeffs[alpha] x^alpha."""
    poly = {(0,) * nvars: 1}
    for _ in range(m):
        new = {}
        for a, ca in poly.items():
            for b, cb in coeffs.items():
                k = tuple(a[i] + b[i] for i in range(nvars))
                new[k] = new.get(k, 0) + ca * cb
        poly = {k: v for k, v in new.items() if v}
    tot = 0
    for a, c in poly.items():
        t = c
        for e in a:
            t *= factorial(e)
        tot += t
    return tot


def census(nvars, maxdeg, crange, mmax, report=6):
    """All integer f on the monomials of total degree <= maxdeg with
       coefficients in crange, having L(f^m) = 0 for m = 1..mmax."""
    mons = [a for a in product(range(maxdeg + 1), repeat=nvars)
            if 0 < sum(a) <= maxdeg]
    survivors, total = [], 0
    for vals in product(crange, repeat=len(mons)):
        if not any(vals):
            continue
        total += 1
        f = {mons[i]: vals[i] for i in range(len(mons)) if vals[i]}
        ok = True
        for m in range(1, mmax + 1):
            if L_of_power(f, m, nvars) != 0:
                ok = False
                break
        if ok:
            survivors.append(f)
    return survivors, total, mons


if __name__ == "__main__":
    print("FACTORIAL CONJECTURE -- census of small integer f with L(f^m)=0, m=1..M")
    print("(constant term excluded: L(f^m)=0 for all m forces it anyway at m->0)")
    print()
    for nvars, maxdeg, cr, mmax in ((2, 2, range(-2, 3), 3),
                                    (2, 2, range(-3, 4), 4),
                                    (3, 1, range(-3, 4), 4),
                                    (3, 2, range(-1, 2), 3)):
        s, tot, mons = census(nvars, maxdeg, cr, mmax)
        print(f"  n={nvars}, deg<={maxdeg}, coeffs in [{min(cr)},{max(cr)}], "
              f"m=1..{mmax}: {len(s)} survivors out of {tot}")
        print(f"     monomials: {mons}")
        for f in s[:4]:
            print(f"       {f}")
        if len(s) > 4:
            print(f"       ... and {len(s)-4} more")
    print()
    print("""READING.  Any survivor is a polynomial whose first M moments vanish;
FC asserts none survives ALL m.  A nonempty census at small M is expected and
is not evidence against FC -- it is the finite window.  What matters is
whether survivors persist as M grows, and whether they are supported in
degree >= 3, which is exactly where the determinacy argument above stops
working.""")
