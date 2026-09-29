#!/usr/bin/env python3
"""The seed-1 test of positive-density log-time convergence.

For an odd x whose Syracuse orbit x = x_0, ..., x_n = 1 has word (a_1..a_n), A = sum a_j:
    1/x = 3^n 2^-A Pi(x),   Pi(x) = prod_(j<n) (1 + 1/(3 x_j))  (exact telescoping identity),
and, the x_j (j < n) being distinct odd integers for a first arrival, Pi(x) <= exp((1/3) sum_(k<=n) 1/(2k-1)) <= e^(1/3) n^(1/6).
Hence the undiscounted harmonic mass G_n = sum_(first arrivals at depth n) 1/x satisfies G_n <= 1.3956 n^(1/6) H_n
with H_n = 3^n mu_n(1 mod 3^n) the harmonic mass (Theorem A of the session note).  Positive lower density of
{x : tau(x) <= C ln x} forces sum_(n <= C ln X) G_n >= (c/2) ln(X/X_0), hence limsup n^(1/6) H_n > 0.

This script: (i) Pi(x) over all odd x <= 2 10^6 (max, argmax, mean), (ii) first-arrival versus all-arrival masses
H_n^first, H_n, and G_n^first at small depth (truncated words), (iii) the split of the depth-n layer of the tree of 1
by residue mod 3 (the exact recursion H_(n+1) = H_n^(1) + 2 H_n^(2)).
Session: opus, collatz-poset-dag-20260927 (S19), 2026-09-28.
Run: python 04-computation/experiments/mazur_seed1_test_20260928.py
"""
from __future__ import annotations

import math
import sys
from fractions import Fraction

sys.path.insert(0, "04-computation/experiments")
from mazur_harmonic_mass_20260928 import mu_exact  # noqa: E402

LN2, LN3 = math.log(2), math.log(3)


def orbit_stats(x: int):
    n = A = 0
    s = 0.0
    y = x
    while y != 1:
        s += 1.0 / y
        t = 3 * y + 1
        a = (t & -t).bit_length() - 1
        y = t >> a
        n += 1
        A += a
    return n, A, s


if __name__ == "__main__":
    print("== (i) the discount product Pi(x) = 2^A/(3^n x) over odd x <= 2 10^6 ==")
    best = (0.0, 0)
    tot = 0.0
    cnt = 0
    worst_bound = 0.0
    for x in range(1, 2 * 10 ** 6 + 1, 2):
        n, A, s = orbit_stats(x)
        if n == 0:
            continue
        lp = A * LN2 - n * LN3 - math.log(x)
        pi = math.exp(lp)
        # the bound e^(1/3) n^(1/6) and the tighter exp(s/3)
        assert lp <= s / 3 + 1e-9 and s <= 0.5 * math.log(n) + 1 + 1e-12, (x, n, s)
        worst_bound = max(worst_bound, pi / (math.exp(1 / 3) * n ** (1 / 6)))
        tot += pi
        cnt += 1
        if pi > best[0]:
            best = (pi, x)
    print(f"   max Pi = {best[0]:.5f} at x = {best[1]}; mean Pi = {tot / cnt:.5f}; max of Pi / (e^(1/3) n^(1/6)) = {worst_bound:.4f} (the bound holds with room)")
    for x in (3, 5, 9, 27, 97, 871):
        n, A, s = orbit_stats(x)
        print(f"   x={x}: depth n={n}, A={A}, Pi = {math.exp(A * LN2 - n * LN3 - math.log(x)):.5f}, sum 1/x_j = {s:.4f}")

    print("== (ii) first-arrival versus all-arrival harmonic masses of the tree of 1 (words with valuations <= 16) ==")
    AMAX = 16
    for n in range(1, 7):
        Hfirst = Fraction(0)
        Hall = Fraction(0)
        Gfirst = Fraction(0)
        split = [Fraction(0)] * 3
        stack = [(1, 0, 0, False)]  # value, depth, valuation sum, revisited-1 flag
        while stack:
            x, d, A, rev = stack.pop()
            if d == n:
                w = Fraction(3 ** n, 2 ** A)
                Hall += w
                if not rev:
                    Hfirst += w
                    Gfirst += Fraction(1, x)
                    split[x % 3] += w
                continue
            for a in range(1, AMAX + 1):
                y = 2 ** a * x - 1
                if y % 3 == 0:
                    z = y // 3
                    stack.append((z, d + 1, A + a, rev or z == 1))
        exact = 3 ** n * mu_exact(n).get(1 % 3 ** n, Fraction(0))
        print(f"   n={n}: H_n (all) truncated {float(Hall):.5f} vs exact 3^n mu_n(1) = {float(exact):.5f}; H_n^first {float(Hfirst):.5f}; "
              f"G_n^first = sum 1/x = {float(Gfirst):.5f}; G/H^first = {float(Gfirst / Hfirst):.4f}; layer split by class mod 3 (0,1,2): "
              f"{[round(float(v / Hfirst), 4) for v in split]}; predicted H_(n+1)^first = H^(1) + 2 H^(2) = {float(split[1] + 2 * split[2]):.5f}")
    print("   (H_1 exact: split 1/21, 16/21, 4/21 gives H_2 = 8/7; the all-arrival mass adds the loop's (3/4)^n-weighted earlier layers)")
    print("DONE")
