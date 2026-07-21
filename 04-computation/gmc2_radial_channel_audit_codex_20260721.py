#!/usr/bin/env python3
"""Exact checks for THM-2014, MISTAKE-211, and HYP-8765.

Methodology / Tournament Analysis
---------------------------------
The observables here are scalar Gaussian moment equations with many-to-one
balanced relations.  There is no canonical antisymmetric pair observable, so a
tournament on monomials or charges would discard the cancellation being tested.
Candidate vertex sets explicitly considered: charges, monomials, primitive
circuits, radial channels, Newton-face defects, and proof obligations.  The
faithful state is at least (charge, radial height); charge alone preserves
angular balance and destroys the factorial weight.  Tournament Analysis is
therefore not used as mathematical evidence.  The only useful orientation is
the proof-obligation path

    moment identity -> circuit localization -> channel determinant
    -> pair radical -> one-sided nullcone -> GMC(2).

All algebra below is exact SymPy arithmetic; there is no floating-point test
except decimal display of exact rational ratios in the final stress table.
"""

from __future__ import annotations

from itertools import product
from math import factorial

import sympy as sp


u = sp.symbols("u")


def laplace_poly(expr: sp.Expr) -> sp.Expr:
    """L(sum c_j u^j) = sum c_j j! exactly."""
    poly = sp.Poly(sp.expand(expr), u)
    return sp.expand(
        sum(coeff * factorial(power[0]) for power, coeff in poly.terms())
    )


def compositions(total: int, parts: int):
    """Yield weak compositions of total into parts slots."""
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first,) + tail


def wick_moment(terms, m: int) -> sp.Expr:
    """Moment of sum_i coeff_i Z^a_i W^b_i by exact Wick balance."""
    ans = 0
    for counts in compositions(m, len(terms)):
        zdeg = sum(r * term[0] for r, term in zip(counts, terms))
        wdeg = sum(r * term[1] for r, term in zip(counts, terms))
        if zdeg != wdeg:
            continue
        multinomial = factorial(m)
        monomial = 1
        for r, (_, _, coeff) in zip(counts, terms):
            multinomial //= factorial(r)
            monomial *= coeff**r
        ans += multinomial * factorial(zdeg) * monomial
    return sp.expand(ans)


def slice_moment(bpoly: sp.Expr, gamma: sp.Expr, m: int) -> sp.Expr:
    """Equation (1) of THM-2014."""
    ans = 0
    for k in range(m // 2 + 1):
        n = m - 2 * k
        ans += (
            factorial(m)
            * gamma**k
            * laplace_poly(u**k * bpoly**n)
            / (factorial(k) ** 2 * factorial(n))
        )
    return sp.expand(ans)


def check_first_return_collision():
    a, b, c = sp.symbols("a b c", nonzero=True)
    terms = [(6, 0, a), (0, 2, b), (0, 18, c)]
    m4 = wick_moment(terms, 4)
    expected = 4 * factorial(6) * a * b**3 + 4 * factorial(18) * a**3 * c
    assert sp.expand(m4 - expected) == 0
    c_cancel = -sp.Rational(factorial(6), factorial(18)) * b**3 / a**2
    assert sp.simplify(m4.subs(c, c_cancel)) == 0
    return m4, c_cancel


def check_return_two_faker():
    A = -sp.Rational(15, 2) + sp.Rational(3, 2) * sp.sqrt(3) * sp.I
    B = 9 - 3 * sp.sqrt(3) * sp.I
    H = u * (u**2 + A * u + B)
    values = [sp.simplify(laplace_poly(H**j)) for j in range(1, 4)]
    expected3 = 3888 * sp.I * (3 * sp.sqrt(3) - 2 * sp.I)
    assert values[:2] == [0, 0]
    assert sp.simplify(values[2] - expected3) == 0
    return values


def check_return_three_faker():
    a = sp.symbols("a")
    b = sp.symbols("b")
    terms = [(0, 1, a), (0, 4, b), (1, 2, 1), (2, 0, 1)]
    b_sub = -a**2 / 12 - a / 2 - 1
    f = a**4 - 15 * a**3 - 240 * a**2 - 945 * a - 1260

    m3 = sp.expand(wick_moment(terms, 3).subs(b, b_sub))
    m6 = sp.expand(wick_moment(terms, 6).subs(b, b_sub))
    m9 = sp.expand(wick_moment(terms, 9).subs(b, b_sub))
    assert m3 == 0
    assert sp.rem(sp.Poly(m6, a), sp.Poly(f, a)).as_expr() == 0

    rem9 = sp.rem(sp.Poly(sp.cancel(m9 / 60480), a), sp.Poly(f, a)).as_expr()
    expected = 360 * (
        3560 * a**3 + 38530 * a**2 + 137613 * a + 174552
    )
    assert sp.expand(rem9 - expected) == 0
    resultant = sp.resultant(f, expected, a)
    expected_resultant = 8335804094167484712960000
    assert abs(int(resultant)) == expected_resultant
    return f, rem9, resultant


def check_slice_formula_and_stress():
    # Direct Wick vs THM-2014 formula on a nontrivial quadratic radial middle.
    z, w = sp.symbols("z w")
    del z, w  # exponents, not symbolic expansion, drive wick_moment.
    bpoly = u**2 - 4 * u + 2
    gamma = sp.Integer(1)
    direct_terms = [
        (1, 0, 1),
        (0, 1, 1),
        (2, 2, 1),
        (1, 1, -4),
        (0, 0, 2),
    ]
    for m in range(1, 9):
        assert sp.expand(wick_moment(direct_terms, m) - slice_moment(bpoly, gamma, m)) == 0

    ratios = []
    for m in (10, 20, 40, 80):
        radial = laplace_poly(bpoly**m)
        total = slice_moment(bpoly, gamma, m)
        ratios.append((m, sp.Rational(total, radial)))
    return ratios


def main():
    m4, c_cancel = check_first_return_collision()
    values2 = check_return_two_faker()
    f, rem9, resultant = check_return_three_faker()
    ratios = check_slice_formula_and_stress()

    print("GMC2 RADIAL-CHANNEL AUDIT — exact gates")
    print("gate first-return star collision: PASS")
    print(f"  E4 = {m4}")
    print(f"  cancelling c = {c_cancel}")
    print("gate R=2 / death-at-6 faker: PASS")
    print(f"  L(H),L(H^2),L(H^3) = {values2}")
    print("gate R=3 / death-at-9 faker: PASS")
    print(f"  f(a) = {f}")
    print(f"  E9/60480 mod f = {rem9}")
    print(f"  resultant = {resultant}")
    print("gate THM-2014 charge formula vs direct Wick, m=1..8: PASS")
    print("quadratic stress b=u^2-4u+2, gamma=1: M_m/L(b^m)")
    for m, ratio in ratios:
        print(f"  m={m:2d}: {sp.N(ratio, 12)}")
    print("assumption challenge: charge-only quotient REJECTED; retain (charge,height)")
    print("Tournament Analysis: not mathematically applicable; no canonical pair orientation")
    print("ALL EXACT GATES PASS")


if __name__ == "__main__":
    main()
