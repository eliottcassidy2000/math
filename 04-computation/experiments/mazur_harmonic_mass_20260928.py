#!/usr/bin/env python3
"""The harmonic mass of the Syracuse tree of 1, layer by layer, equals 3^n times the
probability that the Syracuse random variable Y_n (Tao; Mazur's mu_n) is 1 mod 3^n:

    H_n := sum over odd x of Syracuse depth exactly n (to 1) of 3^n / 2^(A(x))  =  3^n mu_n(1 mod 3^n)  =  (3/2) rho_n(1),

where A(x) is the total valuation along the orbit of x to 1.  Direct enumeration of the tree layers
(valuations truncated) against the exact residue recursion, and the values of rho_n(1) to n = 9.

Session: opus, collatz-poset-dag-20260927 (S19), 2026-09-28.
Run: python 04-computation/experiments/mazur_harmonic_mass_20260928.py
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction


def v2(x: int) -> int:
    return (x & -x).bit_length() - 1


def mu_exact(q: int) -> dict:
    """Exact law of Y_q on Z/3^q: Y_0 = 0, Y_(k+1) = 2^-A (3 Y_k + 1), A geometric(1/2) on {1, 2, ...}."""
    mu = {0: Fraction(1)}
    for lev in range(1, q + 1):
        mod = 3 ** lev
        L = 2 * 3 ** (lev - 1)  # order of 2 modulo 3^lev
        pr = {r: Fraction(2 ** (L - r), 2 ** L - 1) for r in range(1, L + 1)}  # P(A = r mod L), r = 1..L
        inv = {r: pow(2, -r, mod) for r in range(1, L + 1)}
        new = defaultdict(Fraction)
        for y, p in mu.items():
            base = (3 * y + 1) % mod
            for r in range(1, L + 1):
                new[(inv[r] * base) % mod] += p * pr[r]
        mu = dict(new)
    return mu


def mu_float(q: int) -> dict:
    mu = {0: 1.0}
    for lev in range(1, q + 1):
        mod = 3 ** lev
        L = 2 * 3 ** (lev - 1)
        corr = 1.0 - 2.0 ** (-L)  # P(A = r mod L) = 2^-r / (1 - 2^-L); the denominator matters for small L (L = 2: 2/3, 1/3)
        pr = {r: 2.0 ** (-r) / corr for r in range(1, L + 1)}
        inv = {r: pow(2, -r, mod) for r in range(1, L + 1)}
        new = defaultdict(float)
        for y, p in mu.items():
            base = (3 * y + 1) % mod
            for r in range(1, L + 1):
                new[(inv[r] * base) % mod] += p * pr[r]
        mu = dict(new)
    return mu


def tree_layer_mass(n: int, amax: int) -> Fraction:
    """Direct: sum over words (a_1..a_n) with a_i <= amax admissible at the seed 1 (inverse orbit
    integral, all intermediate values odd positive; the seed 1 is reached with the given valuations)
    of 3^n 2^-A.  Admissibility is decided by actually forming the inverse orbit."""
    total = Fraction(0)
    # depth-first over inverse steps from 1
    stack = [(1, 0, 0)]  # (current value, depth, valuation sum)
    while stack:
        x, d, A = stack.pop()
        if d == n:
            total += Fraction(3 ** n, 2 ** A)
            continue
        for a in range(1, amax + 1):
            y = 2 ** a * x - 1
            if y % 3 == 0:
                z = y // 3
                # z is odd and positive; its Syracuse image is x with valuation a exactly
                stack.append((z, d + 1, A + a))
    return total


if __name__ == "__main__":
    print("== harmonic mass of the Syracuse tree of 1: H_n = sum_(depth n) 3^n 2^-A  versus 3^n mu_n(1) ==")
    for n in range(1, 7):
        mu = mu_exact(n)
        pred = 3 ** n * mu.get(1 % 3 ** n, Fraction(0))
        # truncated direct sums converge to the exact value as amax grows (geometric tails)
        direct = [tree_layer_mass(n, amax) for amax in (8, 14, 20)] if n <= 4 else [tree_layer_mass(n, 12)]
        print(f"   n={n}: 3^n mu_n(1) = {float(pred):.6f} (exact {pred if pred.denominator < 10**12 else float(pred)}); direct truncated layer sums {[round(float(v), 6) for v in direct]}")
        if n <= 4:
            assert abs(float(direct[-1]) - float(pred)) < 1e-4, (n, direct[-1], pred)
    print("   (the truncated sums approach the exact value: the identity H_n = 3^n mu_n(1) holds; the residue 1 is an ordinary unit of Z/3^n)")
    print("   rho_n(1) = (2/3) 3^n mu_n(1) and the units' mean of rho_n is 1:")
    vals = []
    for n in range(1, 9):
        mu = mu_float(n)
        r1 = (2 / 3) * 3 ** n * mu.get(1, 0.0)
        units = [(2 / 3) * 3 ** n * mu.get(y, 0.0) for y in range(3 ** n) if y % 3]
        vals.append(r1)
        mod = 3 ** n
        argmax = max((y for y in range(mod) if y % 3), key=lambda y: mu.get(y, 0.0))
        rm1 = (2 / 3) * mod * mu.get((-1) % mod, 0.0)
        rm5 = (2 / 3) * mod * mu.get((-5) % mod, 0.0)
        rm17 = (2 / 3) * mod * mu.get((-17) % mod, 0.0)
        print(f"   n={n}: rho_n(1) = {r1:.5f}; H_n = {1.5 * r1:.5f}; rho_n on units: min {min(units):.4f}, max {max(units):.4f} at y = {argmax} (= -1 mod 3^n: {argmax == (-1) % mod}); "
              f"rho_n(-1) = {rm1:.4f} vs (2/3)(3/2)^n = {(2/3)*1.5**n:.4f}; rho_n(-5) = {rm5:.4f} vs (2/3)(9/8)^(n/2) = {(2/3)*(9/8)**(n/2):.4f}; rho_n(-17) = {rm17:.4f} vs (2/3)(2187/2048)^(n/7) = {(2/3)*(2187/2048)**(n/7):.4f}")
    print("   Mazur chooses a seed by averaging over residue classes; the seed 1 itself has the masses H_n above, which the mixing estimate says converge; positivity of the limit for the seed 1 is not claimed by the paper.")
    print("DONE")
