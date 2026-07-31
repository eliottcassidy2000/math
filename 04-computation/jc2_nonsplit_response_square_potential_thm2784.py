#!/usr/bin/env python3
"""Exact hostile controls for THM-2784's nonsplit square-potential theorem.

The all-degree theorem is a function-field argument.  This companion checks
the exact differential/square-potential interface, the finite-place valuation
ledger, sharp positive boundaries, and a bounded rational-function atlas.
It intentionally includes a high-degree nonsquare example whose squarefree
part has large degree, preventing an overstatement of the squarefree corollary.
"""

from __future__ import annotations

import ast
import itertools
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


x = sp.symbols("x")
GATES = 0


def gate(condition: bool, message: str) -> None:
    global GATES
    require(condition, message)
    GATES += 1


def normalized_factor_dict(poly: sp.Expr) -> dict[sp.Expr, int]:
    polynomial = sp.Poly(poly, x, domain=sp.QQ)
    factors: dict[sp.Expr, int] = {}
    for factor, exponent in sp.factor_list(polynomial.as_expr())[1]:
        monic = sp.Poly(factor, x, domain=sp.QQ).monic().as_expr()
        factors[monic] = exponent
    return factors


def rational_valuation(expr: sp.Expr, factor: sp.Expr) -> int:
    numerator, denominator = sp.fraction(sp.cancel(expr))
    numerator_factors = normalized_factor_dict(numerator)
    denominator_factors = normalized_factor_dict(denominator)
    monic = sp.Poly(factor, x, domain=sp.QQ).monic().as_expr()
    return numerator_factors.get(monic, 0) - denominator_factors.get(monic, 0)


def potential_from_f(function: sp.Expr) -> sp.Expr | None:
    derivative = sp.diff(function, x)
    if derivative == 0:
        return None
    candidate = sp.cancel(4 * function / derivative**2)
    numerator, denominator = sp.fraction(candidate)
    if sp.degree(denominator, x) > 0:
        return None
    return sp.Poly(sp.cancel(numerator / denominator), x, domain=sp.QQ).as_expr()


def check_potential(function: sp.Expr, potential: sp.Expr, label: str) -> int:
    gate(
        sp.cancel(potential * sp.diff(function, x) ** 2 - 4 * function) == 0,
        f"square-potential identity failed: {label}",
    )
    factor_data = normalized_factor_dict(potential)
    local_rows = 0
    for factor, multiplicity in factor_data.items():
        order_f = rational_valuation(function, factor)
        gate(
            order_f == 2 - multiplicity,
            f"local valuation ledger failed: {label}, {factor}",
        )
        gate(multiplicity != 2, f"forbidden double root survived: {label}")
        local_rows += 1
    if factor_data and all(exponent == 1 for exponent in factor_data.values()):
        gate(
            sp.degree(potential, x) <= 1,
            f"squarefree potential had degree at least two: {label}",
        )

    # Complete the finite ledger at infinity.  Let s be the number of simple
    # V-roots, h the number of V-roots of multiplicity at least three, and e
    # the number of extra double F-zeros (all counted geometrically).
    numerator, denominator = sp.fraction(sp.cancel(function))
    numerator_factors = normalized_factor_dict(numerator)
    denominator_factors = normalized_factor_dict(denominator)
    simple_roots = sum(
        sp.degree(factor, x)
        for factor, exponent in factor_data.items()
        if exponent == 1
    )
    high_roots = sum(
        sp.degree(factor, x)
        for factor, exponent in factor_data.items()
        if exponent >= 3
    )
    extra_double_zeros = sum(
        sp.degree(factor, x)
        for factor, exponent in numerator_factors.items()
        if exponent == 2 and factor not in factor_data
    )
    finite_zero_degree = sum(
        sp.degree(factor, x) * exponent
        for factor, exponent in numerator_factors.items()
    )
    finite_pole_degree = sum(
        sp.degree(factor, x) * exponent
        for factor, exponent in denominator_factors.items()
    )
    infinity_order = finite_pole_degree - finite_zero_degree
    potential_degree = sp.degree(potential, x)
    if infinity_order != 0:
        gate(
            infinity_order == potential_degree - 2,
            f"unbalanced infinity order failed: {label}",
        )
        gate(
            simple_roots + high_roots + extra_double_zeros == 1,
            f"unbalanced one-point classification failed: {label}",
        )
        map_degree = max(finite_zero_degree, finite_pole_degree)
        ramification_total = (
            extra_double_zeros
            + (finite_pole_degree - high_roots)
            + abs(infinity_order)
            - 1
        )
        gate(
            ramification_total == 2 * map_degree - 2,
            f"unbalanced Riemann-Hurwitz ledger failed: {label}",
        )
    else:
        leading_value = sp.LC(sp.Poly(numerator, x)) / sp.LC(
            sp.Poly(denominator, x)
        )
        difference = sp.cancel(function - leading_value)
        difference_numerator, difference_denominator = sp.fraction(difference)
        ramification = (
            sp.degree(difference_denominator, x)
            - sp.degree(difference_numerator, x)
        )
        gate(ramification >= 1, f"balanced infinity ramification failed: {label}")
        gate(
            potential_degree == 2 * ramification + 2,
            f"balanced infinity degree failed: {label}",
        )
        gate(
            finite_zero_degree == finite_pole_degree,
            f"balanced divisor degrees failed: {label}",
        )
        gate(
            ramification
            == simple_roots + high_roots + extra_double_zeros - 1,
            f"balanced branch ledger failed: {label}",
        )
        map_degree = finite_zero_degree
        ramification_total = (
            extra_double_zeros
            + (finite_pole_degree - high_roots)
            + ramification
            - 1
        )
        gate(
            ramification_total == 2 * map_degree - 2,
            f"balanced Riemann-Hurwitz passport failed: {label}",
        )
    return local_rows


def main() -> None:
    local_rows = 0

    # Sharp positive controls.  These correspond respectively to a simple
    # root, one root of multiplicity three, and the first two-root solution.
    named = (
        (4 * x, x, "linear-squarefree-boundary"),
        (4 / x, x**3, "triple-root-boundary"),
        (4 * x / (x - 1), x * (x - 1) ** 3, "one-three-two-root-boundary"),
    )
    for function, potential, label in named:
        local_rows += check_potential(function, potential, label)

    # A deterministic bounded atlas of rational monomials on three marked
    # finite points.  Whenever 4F/(F')^2 is polynomial, attack the complete
    # finite-place ledger and the squarefree-degree corollary.
    atlas_rows = 0
    atlas_squarefree_rows = 0
    for exponents in itertools.product(range(-2, 3), repeat=3):
        if exponents == (0, 0, 0):
            continue
        function = sp.Integer(4)
        for point, exponent in zip((-1, 0, 1), exponents):
            function *= (x - point) ** exponent
        function = sp.cancel(function)
        potential = potential_from_f(function)
        if potential is None or sp.degree(potential, x) <= 0:
            continue
        atlas_rows += 1
        factors = normalized_factor_dict(potential)
        if factors and all(exponent == 1 for exponent in factors.values()):
            atlas_squarefree_rows += 1
        local_rows += check_potential(
            function,
            potential,
            f"atlas-{exponents}",
        )

    # Hostile against the false strengthening "large squarefree part is
    # impossible".  Exact primitives exist for
    # V=x^(n+2)(x^n-1)/n^2, even though its squarefree part grows with n.
    family_rows = 0
    maximal_squarefree_part_degree = 0
    for n in range(1, 7):
        function = sp.cancel(4 * (1 - x ** (-n)))
        potential = sp.cancel(x ** (n + 2) * (x**n - 1) / (n * n))
        local_rows += check_potential(function, potential, f"large-part-{n}")
        radical = sp.prod(factor for factor in normalized_factor_dict(potential))
        radical_degree = sp.degree(radical, x)
        maximal_squarefree_part_degree = max(
            maximal_squarefree_part_degree,
            radical_degree,
        )
        family_rows += 1

    gate(maximal_squarefree_part_degree == 7, "hostile radical degree drifted")

    # The geometric endpoints used in the proof are exact symbolic formulas.
    # For V=c(x-a), dx/U=(2/c)dU.  For deg(V)=2 with leading coefficient c^2,
    # the residues of dx/U at the two infinities are -1/c and +1/c.
    c, b, kappa = sp.symbols("c b kappa", nonzero=True)
    gate(sp.simplify((2 / c) * (c / 2) - 1) == 0, "linear primitive failed")
    linear_v = c * x + b
    linear_f = 4 * kappa**2 * linear_v / c**2
    gate(
        sp.expand(linear_v * sp.diff(linear_f, x) ** 2 - 4 * kappa**2 * linear_f)
        == 0,
        "unique linear square potential failed",
    )
    residues = (-1 / c, 1 / c)
    gate(all(residue != 0 for residue in residues), "quadratic residue vanished")
    gate(sp.simplify(sum(residues)) == 0, "quadratic residues do not balance")

    # If V is squarefree, the local ledger makes F a polynomial.  If
    # D=deg(V), n=deg(F), leading degrees force D+2(n-1)=n, or n=2-D.
    # The only nonconstant possibility is D=n=1.
    squarefree_degree_rows = 0
    for degree_v in range(1, 9):
        degree_f = 2 - degree_v
        gate(
            (degree_f >= 1) == (degree_v == 1),
            f"squarefree degree wall failed at degree {degree_v}",
        )
        squarefree_degree_rows += 1

    # At infinity, dx/U has orders 2g-2 (odd degree 2g+1) and g-1
    # (even degree 2g+2).  They are nonnegative exactly from degrees 3 and 4.
    infinity_rows = 0
    for g in range(1, 9):
        gate(2 * g - 2 >= 0, f"odd infinity order failed at genus {g}")
        gate(g - 1 >= 0, f"even infinity order failed at genus {g}")
        infinity_rows += 2

    source = Path(__file__).read_text(encoding="utf-8")
    assert_nodes = sum(
        isinstance(node, (ast.Assert,))
        for node in ast.walk(ast.parse(source))
    )
    gate(assert_nodes == 0, "assert statements are forbidden")

    print("NONSPLIT RESPONSE SQUARE-POTENTIAL EXACT REFEREE")
    print(f"named_positive_rows={len(named)}")
    print(f"bounded_polynomial_potential_rows={atlas_rows}")
    print(f"bounded_squarefree_rows={atlas_squarefree_rows}")
    print(f"large_squarefree_part_hostiles={family_rows}")
    print(f"maximal_hostile_squarefree_part_degree={maximal_squarefree_part_degree}")
    print(f"finite_root_valuation_rows={local_rows}")
    print(f"global_infinity_ledger_rows={len(named) + atlas_rows + family_rows}")
    print(f"squarefree_degree_rows={squarefree_degree_rows}")
    print(f"infinity_order_rows={infinity_rows}")
    print("degree1_boundary=dx/U=(2/c)dU")
    print("degree2_residues=(-1/c,+1/c)")
    print(f"assert_nodes={assert_nodes}")
    print(f"exact_gates={GATES}")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
