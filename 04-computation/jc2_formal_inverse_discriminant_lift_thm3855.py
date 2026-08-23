#!/usr/bin/env python3
"""Exact companion for THM-3855's formal inverse-discriminant lift."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.expand(expression) == 0, message)


A, C = sp.symbols("A C")
a, b, c, d = sp.symbols("a b c d")


def discriminant(aa: sp.Expr, bb: sp.Expr, cc: sp.Expr, dd: sp.Expr) -> sp.Expr:
    return sp.expand(
        bb**2 * cc**2
        - 4 * aa * cc**3
        - 4 * bb**3 * dd
        - 27 * aa**2 * dd**2
        + 18 * aa * bb * cc * dd
    )


def homogeneous(expression: sp.Expr, degree: int) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), A, C)
    return sp.Add(
        *[
            coefficient * A**i * C**j
            for (i, j), coefficient in polynomial.terms()
            if i + j == degree
        ]
    )


base = (A, C, 7 * A, -3 * A)
Delta0 = A * (C + 5 * A) * (4 * C + 19 * A) * (3 * C - 17 * A)
zero(discriminant(*base) - Delta0, "four-ray discriminant")

D_symbolic = discriminant(a, b, c, d)
gradients = [
    sp.expand(sp.diff(D_symbolic, variable).subs(dict(zip((a, b, c, d), base))))
    for variable in (a, b, c, d)
]
expected_gradients = [
    -2 * A**2 * (929 * A + 189 * C),
    -2 * A * (189 * A**2 - 49 * A * C - 18 * C**2),
    -2 * A * (294 * A**2 + 27 * A * C - 7 * C**2),
    2 * (81 * A**3 + 63 * A**2 * C - 2 * C**3),
]
for index, (actual, expected) in enumerate(zip(gradients, expected_gradients)):
    zero(actual - expected, f"coefficient gradient {index}")

cubic_monomials = (A**3, A**2 * C, A * C**2, C**3)
gradient_matrix = sp.Matrix(
    [
        [sp.Poly(gradient, A, C).coeff_monomial(monomial) for gradient in gradients]
        for monomial in cubic_monomials
    ]
)
expected_matrix = sp.Matrix(
    [
        [-1858, -378, -588, 162],
        [-378, 98, -54, 126],
        [0, 36, 14, 0],
        [0, 0, 0, -4],
    ]
)
gate(gradient_matrix == expected_matrix, "gradient coefficient matrix")
gate(gradient_matrix.det() == 640000, "gradient determinant")
gradient_inverse = gradient_matrix.inv()
gate(gradient_matrix * gradient_inverse == sp.eye(4), "gradient right inverse")

groebner = sp.groebner(gradients, A, C, order="grevlex")
gate(
    [sp.Poly(polynomial, A, C).monic().as_expr() for polynomial in groebner.polys]
    == [A**3, A**2 * C, A * C**2, C**3],
    "gradient ideal is the cube of the maximal ideal",
)


def right_inverse(form: sp.Expr, degree: int) -> list[sp.Expr]:
    """Lift a homogeneous degree-`degree` form through the cubic gradients."""

    corrections = [sp.S.Zero] * 4
    for (i, j), coefficient in sp.Poly(sp.expand(form), A, C).terms():
        gate(i + j == degree, "right-inverse input is homogeneous")
        if i >= 3:
            basis_index = 0
            multiplier = A ** (i - 3) * C**j
        elif i == 2:
            basis_index = 1
            multiplier = C ** (j - 1)
        elif i == 1:
            basis_index = 2
            multiplier = C ** (j - 2)
        else:
            basis_index = 3
            multiplier = C ** (j - 3)
        for column in range(4):
            corrections[column] += (
                coefficient * gradient_inverse[column, basis_index] * multiplier
            )
    zero(
        sum(gradient * correction for gradient, correction in zip(gradients, corrections))
        - form,
        "homogeneous gradient right inverse",
    )
    return [sp.expand(correction) for correction in corrections]


# The one-place target of THM-3853.  The recursion is performed at lambda=1;
# homogeneity gives the displayed scalar first correction for arbitrary lambda.
target = sp.expand(Delta0 + C**5)
coefficients = list(base)
first_correction = None
for total_degree in range(5, 13):
    current = homogeneous(discriminant(*coefficients), total_degree)
    wanted = C**5 if total_degree == 5 else sp.S.Zero
    correction = right_inverse(wanted - current, total_degree)
    gate(
        all(
            homogeneous(entry, total_degree - 3) == entry
            for entry in correction
        ),
        f"correction degree {total_degree - 3}",
    )
    coefficients = [
        sp.expand(entry + update) for entry, update in zip(coefficients, correction)
    ]
    zero(
        homogeneous(discriminant(*coefficients), total_degree) - wanted,
        f"inverse-discriminant recursion through degree {total_degree}",
    )
    if total_degree == 5:
        first_correction = correction

expected_first = [
    sp.Rational(91449, 40000) * C**2,
    sp.Rational(151263, 40000) * C**2,
    -sp.Rational(194481, 20000) * C**2,
    -sp.Rational(1, 4) * C**2,
]
gate(first_correction == expected_first, "explicit first quadratic correction")

error = sp.expand(discriminant(*coefficients) - target)
gate(
    all(homogeneous(error, degree) == 0 for degree in range(0, 13)),
    "truncated lift has no error through degree twelve",
)
gate(homogeneous(error, 13) != 0, "degree-thirteen hostile remains after truncation")
gate(
    all(entry.subs({A: 0, C: 0}) == 0 for entry in coefficients),
    "all lifted index coefficients remain in the maximal ideal",
)

semantic = {
    "base_packet": "(A,C,7A,-3A)",
    "gradient_matrix_det": 640000,
    "gradient_ideal": "(A,C)^3",
    "target": "Delta0+C^5",
    "formal_lift": "all m-adic orders; corrections begin in degree 2",
    "finite_replay": "degrees 5 through 12; first residual degree 13",
    "index_gate": "all coefficients remain in (A,C), so no unit value",
    "scope": "formal completion only; polynomial algebraization and connected Keller surface open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3855-formal-inverse-discriminant-lift-and-algebraization-gate")
print("base=(A,C,7A,-3A);gradient_det=640000;gradient_ideal=(A,C)^3")
print("target=Delta0+C^5;formal_lift=ALL_ORDERS")
print("first_correction=C^2*(91449/40000,151263/40000,-194481/20000,-1/4)")
print("finite_replay=through_degree_12;first_truncation_residual_degree=13")
print("index=coefficients_remain_in_(A,C);unit_value=NONE")
print("scope=formal_only;polynomial_algebraization_and_connected_Keller_surface_OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
