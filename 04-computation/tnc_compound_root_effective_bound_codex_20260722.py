#!/usr/bin/env python3
"""Lightweight exact referee for THM-2101.

Reproduction:
    python3 04-computation/tnc_compound_root_effective_bound_codex_20260722.py
    python3 -O 04-computation/tnc_compound_root_effective_bound_codex_20260722.py

The symbolic half constructs the compound polynomial from labelled roots and
then rewrites its coefficients in elementary symmetric functions.  The finite
half exhausts a stated small integer-coefficient universe and computes literal
Laurent constant terms by exact convolution.  Neither half is used as a proof
of the all-coefficient theorem; they referee its two load-bearing formulas and
its principal hostile/equality controls.
"""

from __future__ import annotations

from itertools import combinations, product
from math import comb, floor

import sympy as sp


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def bound(m: int, n: int) -> tuple[int, int, int]:
    a, b = sorted((m, n))
    c = comb(a + b, a)
    # Uniform incidence of one root among all a-subsets: C*a/(a+b).
    k = comb(a + b - 1, a - 1)
    return c, k, c


def multiply(p: dict[int, int], q: dict[int, int]) -> dict[int, int]:
    out: dict[int, int] = {}
    for i, x in p.items():
        for j, y in q.items():
            out[i + j] = out.get(i + j, 0) + x * y
    return {e: x for e, x in out.items() if x}


def first_nonzero_ct(f: dict[int, int], cutoff: int) -> tuple[int, int] | None:
    power = {0: 1}
    for m in range(1, cutoff + 1):
        power = multiply(power, f)
        value = power.get(0, 0)
        if value:
            return m, value
    return None


def exponent_of(poly: sp.Expr, generator: sp.Symbol, generators: list[sp.Symbol]) -> int:
    if poly == 0:
        return -1
    pp = sp.Poly(poly, *generators)
    idx = generators.index(generator)
    return max(monom[idx] for monom in pp.monoms())


def symbolic_compound_check(d: int, a: int) -> tuple[int, int, int, list[int], int]:
    b = d - a
    roots = sp.symbols(f"x1:{d + 1}")
    y, t, c, z = sp.symbols("Y t c z")
    subset_products = [sp.prod(roots[i] for i in s) for s in combinations(range(d), a)]
    compound = sp.Poly(sp.prod(y - p for p in subset_products), y)
    size = len(subset_products)
    require(size == comb(d, a), ("compound size", d, a, size))

    coefficient_forms: list[sp.Expr] = []
    formal_generators: list[sp.Symbol] | None = None
    moving: sp.Symbol | None = None
    actual_powers: list[int] = []
    for j, coeff in enumerate(compound.all_coeffs()):
        symmetric, remainder, mapping = sp.symmetrize(coeff, roots, formal=True)
        require(remainder == 0, ("nonsymmetric remainder", d, a, j, remainder))
        generators = [pair[0] for pair in mapping]
        if formal_generators is None:
            formal_generators = generators
            moving = generators[b - 1]
        require(generators == formal_generators, ("generator mismatch", d, a, j))
        coefficient_forms.append(symmetric)
        power_moving = exponent_of(symmetric, moving, generators)
        actual_powers.append(max(power_moving, 0))
        one_sided_bound = floor(j * a / b)
        complement_bound = size - j
        require(
            power_moving <= min(one_sided_bound, complement_bound),
            (
                "dual pole bound",
                d,
                a,
                j,
                power_moving,
                one_sided_bound,
                complement_bound,
            ),
        )

    require(formal_generators is not None and moving is not None, ("missing generators", d, a))
    assert formal_generators is not None and moving is not None  # typing only
    c_size, k, theorem_bound = bound(a, b)
    require(c_size == size, ("bound size", d, a, c_size, size))

    aa, bb = sp.symbols("A B", nonzero=True)
    substitution: dict[sp.Symbol, sp.Expr] = {}
    for i, generator in enumerate(formal_generators, start=1):
        substitution[generator] = aa + bb / t if i == b else sp.Symbol(f"E{i}")
    compound_in_e = sum(
        coefficient_forms[j].subs(substitution) * y ** (size - j)
        for j in range(size + 1)
    )
    # Complement duality aligns pole_order(h_j) <= C-j with Y^(C-j).
    # Therefore Y=c*t clears every coefficient without a uniform t^K factor.
    direct_line = sp.cancel(compound_in_e.subs(y, c * t))
    numerator, denominator = sp.fraction(direct_line)
    require(sp.expand(denominator) == 1, ("line denominator", d, a, denominator))
    line_polynomial = sp.Poly(sp.expand(numerator), t)
    require(not line_polynomial.is_zero, ("zero line polynomial", d, a))
    require(
        line_polynomial.degree() == size,
        ("line degree", d, a, line_polynomial.degree(), size),
    )

    # The divided difference along two order-one lines has at worst a simple
    # pole. Multiplying it by t must therefore be polynomial in t.
    divided = sp.cancel(
        t * (compound_in_e.subs(y, c * t) - compound_in_e.subs(y, z * t))
        / ((c - z) * t)
    )
    divided_numerator, divided_denominator = sp.fraction(divided)
    require(
        sp.expand(divided_denominator) == 1,
        ("divided-difference pole worse than simple", d, a, divided_denominator),
    )
    require(divided_numerator != 0, ("zero divided difference", d, a))
    return size, k, theorem_bound, actual_powers, line_polynomial.degree()


def exhaustive_integer_check(m: int, n: int) -> tuple[int, int, dict[int, int]]:
    _, _, theorem_bound = bound(m, n)
    interior = [q for q in range(-m + 1, n) if q != 0]
    checked = 0
    latest = -1
    latest_witness: dict[int, int] = {}
    for left, right in product((-1, 1), repeat=2):
        for values in product((-1, 0, 1), repeat=len(interior)):
            f = {-m: left, n: right}
            f.update({q: x for q, x in zip(interior, values) if x})
            checked += 1
            answer = first_nonzero_ct(f, theorem_bound)
            require(answer is not None, ("missing return", m, n, f, theorem_bound))
            assert answer is not None  # typing only
            first, _ = answer
            if first > latest:
                latest = first
                latest_witness = dict(sorted(f.items()))
    return checked, latest, latest_witness


def main() -> None:
    print("THM-2101 compound-root effective-bound exact referee")
    print("symbolic compound checks")
    for d, a in ((2, 1), (3, 1), (4, 2), (5, 2)):
        size, k, theorem_bound, powers, line_degree = symbolic_compound_check(d, a)
        print(
            f"  (a,b)=({a},{d-a}) C={size} eta={k} bound={theorem_bound} "
            f"deg Q={line_degree} actual e_b powers by coefficient={powers}; "
            "direct line polynomial and simple-pole divided difference: PASS"
        )

    print("exact exhaustive integer-coefficient universe")
    total = 0
    for m, n in ((1, 1), (1, 2), (1, 3), (1, 4), (2, 2), (2, 3), (3, 3)):
        checked, latest, witness = exhaustive_integer_check(m, n)
        total += checked
        _, _, theorem_bound = bound(m, n)
        print(
            f"  (M,N)=({m},{n}) checked={checked} max_first={latest} "
            f"C-bound={theorem_bound} witness={witness}: PASS"
        )
    print(f"  total exact Laurent polynomials checked={total}")

    print("sharp min(M,N)=1 binomial controls")
    for n in range(1, 9):
        f = {-1: 1, n: 1}
        _, _, theorem_bound = bound(1, n)
        answer = first_nonzero_ct(f, theorem_bound)
        require(answer is not None, ("binomial missing return", n, theorem_bound))
        assert answer is not None  # typing only
        first, value = answer
        require(first == n + 1 == theorem_bound, ("binomial first", n, first, theorem_bound))
        require(value == n + 1, ("binomial value", n, value))
        print(f"  u^-1+u^{n}: first={first}=C, CT={value}: PASS")

    hostile = {-2: -1, -1: 1, 1: 1, 2: 1}
    # u^2+u+u^-1-u^-2 from THM-2070: CT at 1,2,3 vanishes; CT at 4 is -12.
    answer = first_nonzero_ct(hostile, bound(2, 2)[2])
    require(answer == (4, -12), ("dihedral hostile", answer))
    print("THM-2070 dihedral hostile control: first=4, CT=-12: PASS")
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
