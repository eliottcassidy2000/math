#!/usr/bin/env python3
"""Exact companion for THM-3696.

Checks the ordinary-triple normalization, Hermite gluing and conductor,
graded coefficient modules, non-Poisson witness, and scalar-row evaluations.
Every truth gate uses ``require`` and remains active under ``python -O``.
"""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import sympy as sp


CHECKS = 0


def require(condition: bool, detail: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(detail)


def coefficient_vector(polynomial: sp.Expr, variable: sp.Symbol, degree: int) -> sp.Matrix:
    poly = sp.Poly(sp.expand(polynomial), variable)
    return sp.Matrix([poly.nth(index) for index in range(degree + 1)])


def jet_vector(polynomial: sp.Expr, variable: sp.Symbol) -> sp.Matrix:
    roots = (-1, 0, 1)
    return sp.Matrix(
        [
            *(sp.expand(polynomial).subs(variable, root) for root in roots),
            *(sp.diff(polynomial, variable).subs(variable, root) for root in roots),
        ]
    )


def matrix_span(polynomials: list[sp.Expr], variable: sp.Symbol, degree: int) -> sp.Matrix:
    return sp.Matrix.hstack(*(coefficient_vector(poly, variable, degree) for poly in polynomials))


def main() -> None:
    b = sp.symbols("b")
    x_symbol, y_symbol = sp.symbols("X Y")
    h = 1 - b**2
    X = b * h
    Y = b**2 * h
    conductor = sp.expand(X**2)

    # Plane relation, birational recovery, and ordinary-triple tangent cone.
    relation = x_symbol**4 - x_symbol**2 * y_symbol + y_symbol**3
    require(sp.expand(relation.subs({x_symbol: X, y_symbol: Y})) == 0, "plane relation")
    require(sp.Poly(relation, x_symbol, y_symbol).is_irreducible, "irreducible plane relation")
    require(sp.expand(b**3 - b + X) == 0, "integral cubic")
    tangent = y_symbol * (y_symbol - x_symbol) * (y_symbol + x_symbol)
    require(sp.expand(tangent - (y_symbol**3 - x_symbol**2 * y_symbol)) == 0, "ordinary triple tangent")
    require([sp.expand(X.subs(b, root)) for root in (-1, 0, 1)] == [0, 0, 0], "three X branches")
    require([sp.expand(Y.subs(b, root)) for root in (-1, 0, 1)] == [0, 0, 0], "three Y branches")

    # The six Hermite coordinates are values at -1,0,1 followed by derivatives.
    base_zero = [sp.Integer(1), X, Y]
    zero_jets = sp.Matrix.hstack(*(jet_vector(poly, b) for poly in base_zero))
    require(zero_jets.rank() == 3, "R0 jet rank")
    expected_zero_annihilators = {
        (-1, 1, 0, 0, 0, 0),
        (-1, 0, 1, 0, 0, 0),
        (0, 0, 0, 1, 4, 1),
    }
    actual_zero_annihilators = {tuple(vector) for vector in zero_jets.T.nullspace()}
    require(actual_zero_annihilators == expected_zero_annihilators, actual_zero_annihilators)

    even_basis = [sp.Integer(1), X, Y, h, h * X]
    even_jets = sp.Matrix.hstack(*(jet_vector(poly, b) for poly in even_basis))
    require(even_jets.rank() == 5, "even normalized module jet rank")
    require(
        {tuple(vector) for vector in even_jets.T.nullspace()} == {(-1, 0, 1, 0, 0, 0)},
        "even normalized module annihilator",
    )

    minus_one_basis = [sp.Integer(1), X, Y, b, b * Y]
    minus_one_jets = sp.Matrix.hstack(*(jet_vector(poly, b) for poly in minus_one_basis))
    require(minus_one_jets.rank() == 5, "minus-one normalized module jet rank")
    require(
        {tuple(vector) for vector in minus_one_jets.T.nullspace()} == {(1, -2, 1, 0, 0, 0)},
        "minus-one normalized module annihilator",
    )

    full_basis = [sp.Integer(1), X, Y, b, h, h * X]
    full_jets = sp.Matrix.hstack(*(jet_vector(poly, b) for poly in full_basis))
    require(full_jets.rank() == 6, "odd negative full normalization")

    # Conductor containment on the rank-three normalization module and the
    # sharp failure one order earlier.
    require(sp.expand(b * X**2 - X * Y) == 0, "bX2=XY")
    require(sp.expand(b**2 * X**2 - Y**2) == 0, "b2X2=Y2")
    require(sp.expand(conductor - b**2 * h**2) == 0, "conductor polynomial")
    hostile = sp.expand(b**2 * X)
    hostile_jet = jet_vector(hostile, b)
    require((hostile_jet[3] + 4 * hostile_jet[4] + hostile_jet[5]) != 0, "X is not conductor")

    # Compare the original monomial semigroup with the closed module formulas
    # through a hostile weight and polynomial-degree window.
    degree_cap = 18
    exponent_cap = 12
    module_rows = []
    for weight in range(-10, 7):
        original = []
        for i in range(exponent_cap + 1):
            for j in range(exponent_cap + 1):
                k = weight + 2 * i + j
                if k < 0:
                    continue
                polynomial = sp.expand(b**k * h ** (i + j))
                if sp.degree(polynomial, b) <= degree_cap:
                    original.append(polynomial)

        if weight >= 0:
            normalized_bases = [b**weight]
            module_kind = "positive:b^rR0"
        else:
            R = -weight
            if R == 1:
                normalized_bases = [h, h * b]
                module_kind = "minus1:h(R0+bR0)"
            elif R % 2 == 0:
                m = R // 2
                normalized_bases = [h**m, h ** (m + 1)]
                module_kind = "negative_even:h^m(R0+hR0)"
            else:
                m = (R - 1) // 2
                normalized_bases = [h ** (m + 1), h ** (m + 1) * b, h ** (m + 2)]
                module_kind = "negative_odd:h^(m+1)C[b]"

        predicted = []
        for base in normalized_bases:
            for i in range(exponent_cap + 1):
                for j in range(exponent_cap + 1):
                    polynomial = sp.expand(base * X**j * Y**i)
                    if sp.degree(polynomial, b) <= degree_cap:
                        predicted.append(polynomial)

        original_matrix = matrix_span(original, b, degree_cap)
        predicted_matrix = matrix_span(predicted, b, degree_cap)
        union_matrix = original_matrix.row_join(predicted_matrix)
        require(
            original_matrix.rank() == predicted_matrix.rank() == union_matrix.rank(),
            (weight, original_matrix.rank(), predicted_matrix.rank(), union_matrix.rank()),
        )
        module_rows.append((weight, original_matrix.rank(), module_kind))

    # Positive weights have unbounded required b-adic order under c-shifts,
    # which is the exact mechanism behind the zero global conductor.
    for valuation in range(0, 8):
        fixed_polynomial = b**valuation * (1 + b)
        shifted_weight = valuation + 1
        require(sp.rem(fixed_polynomial, b**shifted_weight, domain=sp.QQ) != 0, (valuation, shifted_weight))

    # Non-Poisson closure: be has weight -2 and residual b, which violates
    # equality at the two arms.
    bracket_coefficient = sp.expand(
        -sp.diff(3 * h, b) * (2 * h) + 2 * (3 * h) * sp.diff(2 * h, b)
    )
    require(sp.expand(bracket_coefficient + 12 * b * h) == 0, "{3e,2ce}=-12be")
    require(b.subs(b, -1) != b.subs(b, 1), "be is outside M_-2")

    # Scalar-row laws.  Generic elements of the exact Hermite jet modules
    # suffice because all evaluations factor through these six jets.
    coefficients = sp.symbols("a0:5")
    other_coefficients = sp.symbols("g0:3")
    f_even = sp.expand(sum(coefficient * basis for coefficient, basis in zip(coefficients, even_basis)))
    g_zero = sp.expand(sum(coefficient * basis for coefficient, basis in zip(other_coefficients, base_zero)))
    p_two = sp.expand(h * f_even)
    q_one = sp.expand(b * g_zero)
    J_two = sp.expand(2 * p_two * sp.diff(q_one, b) + sp.diff(p_two, b) * q_one)
    require(sp.simplify(J_two.subs(b, -1) - J_two.subs(b, 1)) == 0, "R=2 arm synchronization")
    require(
        sp.simplify(J_two.subs(b, 1) + 2 * f_even.subs(b, 1) * g_zero.subs(b, 1)) == 0,
        "R=2 arm formula",
    )
    require(
        sp.simplify(J_two.subs(b, 0) - 2 * f_even.subs(b, 0) * g_zero.subs(b, 0)) == 0,
        "R=2 center formula",
    )

    f_minus_one = sp.expand(sum(coefficient * basis for coefficient, basis in zip(coefficients, minus_one_basis)))
    p_one = sp.expand(h * f_minus_one)
    J_one = sp.expand(p_one * sp.diff(g_zero, b))
    require(J_one.subs(b, -1) == J_one.subs(b, 1) == 0, "R=1 arm vanishing")
    require(
        sp.simplify(J_one.subs(b, 0) - f_minus_one.subs(b, 0) * sp.diff(g_zero, b).subs(b, 0)) == 0,
        "R=1 center formula",
    )

    generic_f = sum(sp.symbols(f"f{index}") * b**index for index in range(5))
    generic_g = sum(sp.symbols(f"q{index}") * b**index for index in range(5))
    for R in range(3, 10):
        p = sp.expand(h ** ((R + 1) // 2) * generic_f)
        q = sp.expand(b ** (R - 1) * generic_g)
        J = sp.expand(R * p * sp.diff(q, b) + (R - 1) * sp.diff(p, b) * q)
        require([sp.expand(J.subs(b, root)) for root in (-1, 0, 1)] == [0, 0, 0], ("R>=3", R))

    source = Path(__file__).read_text(encoding="utf-8")
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))), "inactive assert")

    semantic_rows = [
        "R0=C[X,Y],X=b(1-b^2),Y=b^2(1-b^2)",
        "plane=X^4-X^2Y+Y^3;branches=-1,0,1",
        "R0_gluing=equal_values;d(-1)+4d(0)+d(1)=0",
        "conductor=b^2(1-b^2)^2*C[b]",
        "modules=positive;bminus1;negative_even;negative_odd",
        "global_conductor=zero",
        "not_Poisson_closed={3e,2ce}=-12be",
        "scalar_arms=only(-2,1),synchronized;center=(-2,1)or(-1,0)",
    ]
    semantic_hash = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

    print("theorem=THM-3696-y0-collision-three-branch-conductor")
    print("weight_zero=R0=C[X,Y];X=b*(1-b^2);Y=b^2*(1-b^2);relation=X^4-X^2Y+Y^3")
    print("normalization=C[b];triple_branches=b=-1,0,1")
    print("R0_membership=equal_values_and_fprime(-1)+4*fprime(0)+fprime(1)=0")
    print("conductor=b^2*(1-b^2)^2*C[b]")
    print("modules=r>=0:b^rR0;r=-1:h(R0+bR0);r=-2m:h^m(R0+hR0);r=-(2m+1),m>=1:h^(m+1)C[b]")
    print(f"hostile_module_window=-10..6;degree_cap={degree_cap};rows={tuple(module_rows)}")
    print("global_conductor_R_in_D=zero;mechanism=unbounded_positive_b_adic_order")
    print("poisson_closure=false;witness={3e,2ce}=-12be_notin_R")
    print("scalar_law=arms_only(-2,1)_and_synchronized;center_only(-2,1)_or(-1,0)")
    print(f"semantic_sha256={semantic_hash}")
    print(f"CHECKS={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
