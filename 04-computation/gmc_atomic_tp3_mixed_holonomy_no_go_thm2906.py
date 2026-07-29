#!/usr/bin/env python3
"""Exact companion for THM-2906.

The two cells below have positive adjacent-difference tensor data and
strictly positive local response TP3 determinants, but their mixed
endpoint-holonomy determinants have opposite signs.

All truth-bearing gates use explicit exceptions so optimized Python runs
the same verification.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations_with_replacement, product
from math import factorial

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def adjacent_tensor(indices):
    """Return L(prod_j d_(indices[j])) by literal 2^k expansion."""

    indices = tuple(sorted(indices))
    value = Fraction(0)
    for mask in range(1 << len(indices)):
        degrees = []
        sign = 1
        for position, index in enumerate(indices):
            if mask & (1 << position):
                degrees.append(index + 1)
            else:
                degrees.append(index)
                sign = -sign
        term = Fraction(factorial(sum(degrees)), 1)
        for degree in degrees:
            term /= factorial(degree)
        value += sign * term
    return sp.Rational(value.numerator, value.denominator)


def moment(*vectors):
    """Multilinear adjacent-difference moment of sparse coefficient maps."""

    return sp.expand(
        sum(
            sp.prod(coefficient for _, coefficient in terms)
            * adjacent_tensor(tuple(index for index, _ in terms))
            for terms in product(*(tuple(vector.items()) for vector in vectors))
        )
    )


def determinant_three(matrix):
    (a, b, c), (d, e, f), (g, h, i) = matrix
    return (
        a * (e * i - f * h)
        - b * (d * i - f * g)
        + c * (d * h - e * g)
    )


def local_response_tp3(vector):
    return determinant_three(
        [
            [
                moment({1 + row: sp.Integer(1)}, *([vector] * power))
                for power in (1, 2, 3)
            ]
            for row in range(3)
        ]
    )


def cell(x_value, y_value):
    U = {1: sp.Integer(1), 3: sp.Integer(x_value)}
    V = {2: sp.Integer(1), 3: sp.Integer(y_value)}

    g_0 = moment(U, U)
    g_1 = moment(U, V)
    g_2 = moment(V, V)

    t_0 = moment(U, U, U)
    t_1 = moment(U, U, V)
    t_2 = moment(U, V, V)
    t_3 = moment(V, V, V)

    A = [
        moment(*([U] * (4 - upper_count) + [V] * upper_count))
        for upper_count in range(5)
    ]

    r_left = sp.cancel((2 * A[1] * g_0 - A[0] * g_1) / g_0**2)
    r_right = sp.cancel((2 * A[3] * g_2 - A[4] * g_1) / g_2**2)
    holonomy = sp.expand(g_0**2 * g_2**2 * (r_left - r_right))

    invariant_one = sp.expand(
        3 * t_1 * g_0 * g_2
        - t_3 * g_0**2
        - 2 * t_0 * g_1 * g_2
    )
    invariant_two = sp.expand(
        3 * t_2 * g_0 * g_2
        - 2 * t_3 * g_1 * g_0
        - t_0 * g_2**2
    )

    return {
        "g": (g_0, g_1, g_2),
        "gram": sp.expand(g_0 * g_2 - g_1**2),
        "r_left": r_left,
        "r_right": r_right,
        "holonomy": holonomy,
        "invariants": (invariant_one, invariant_two),
        "delta_U": local_response_tp3(U),
        "delta_V": local_response_tp3(V),
    }


def polynomial_digest(polynomial):
    records = [
        f"{monomial[0]}:{monomial[1]}:{coefficient}"
        for monomial, coefficient in polynomial.terms()
    ]
    return sha256(("\n".join(records) + "\n").encode()).hexdigest()


def main():
    print("THM-2906 exact mixed-endpoint holonomy no-go")

    atomic_cells = 0
    for order in (2, 3, 4):
        for indices in combinations_with_replacement(range(1, 4), order):
            require(
                adjacent_tensor(indices) > 0,
                f"nonpositive atomic tensor: {indices}",
            )
            atomic_cells += 1

    x, y, z, w = sp.symbols("x y z w")
    U = {1: sp.Integer(1), 3: x}
    V = {2: sp.Integer(1), 3: y}
    g_0 = moment(U, U)
    g_1 = moment(U, V)
    g_2 = moment(V, V)
    A = [
        moment(*([U] * (4 - upper_count) + [V] * upper_count))
        for upper_count in range(5)
    ]

    J = sp.Poly(
        sp.expand(
            (2 * A[1] * g_0 - A[0] * g_1) * g_2**2
            - (2 * A[3] * g_2 - A[4] * g_1) * g_0**2
        ),
        x,
        y,
    )
    positive_terms = sum(1 for _, coefficient in J.terms() if coefficient > 0)
    negative_terms = sum(1 for _, coefficient in J.terms() if coefficient < 0)
    require(
        (J.degree(x), J.degree(y), len(J.terms())) == (5, 5, 30),
        "endpoint polynomial profile changed",
    )
    require(
        (positive_terms, negative_terms) == (13, 17),
        "endpoint polynomial sign census changed",
    )
    selected = (
        J.coeff_monomial(1),
        J.coeff_monomial(x),
        J.coeff_monomial(y),
        J.coeff_monomial(x * y),
    )
    require(
        selected == (273888, -5839200, 5559264, -82033920),
        "selected endpoint coefficients changed",
    )

    q = g_0 + 2 * g_1 * z + g_2 * z**2
    f = sum(sp.binomial(4, index) * A[index] * z**index for index in range(5))
    q_right = sp.expand(w**2 * q.subs(z, 1 / w))
    f_right = sp.expand(w**4 * f.subs(z, 1 / w))
    r_left = sp.cancel((2 * A[1] * g_0 - A[0] * g_1) / g_0**2)
    r_right = sp.cancel((2 * A[3] * g_2 - A[4] * g_1) / g_2**2)
    require(
        sp.cancel(sp.diff(f / q, z).subs(z, 0) - 2 * r_left) == 0,
        "left quotient derivative identity changed",
    )
    require(
        sp.cancel(sp.diff(f_right / q_right, w).subs(w, 0) - 2 * r_right)
        == 0,
        "right quotient derivative identity changed",
    )

    first = cell(1, 1)
    second = cell(1, 2)
    require(
        first
        == {
            "g": (30, 37, 46),
            "gram": 11,
            "r_left": sp.Rational(73097746, 75),
            "r_right": sp.Rational(515413080, 529),
            "holonomy": 610878432,
            "invariants": (141000, 169320),
            "delta_U": 36467293512,
            "delta_V": 86548240800,
        },
        "first hostile cell changed",
    )
    require(
        second
        == {
            "g": (30, 61, 126),
            "gram": 59,
            "r_left": sp.Rational(142883488, 75),
            "r_right": sp.Rational(120168140, 63),
            "holonomy": -33115086144,
            "invariants": (957960, 2436120),
            "delta_U": 36467293512,
            "delta_V": 2255452579800,
        },
        "second hostile cell changed",
    )
    require(
        first["holonomy"] > 0 > second["holonomy"],
        "mixed holonomy no longer reverses sign",
    )
    require(
        all(value > 0 for value in first["invariants"] + second["invariants"]),
        "a hostile cell accidentally entered the cubic-null ideal",
    )
    require(
        min(
            first["delta_U"],
            first["delta_V"],
            second["delta_U"],
            second["delta_V"],
        )
        > 0,
        "a local TP3 control lost positivity",
    )

    print(f"atomic_positive_cells={atomic_cells}")
    print(
        "endpoint_polynomial="
        f"degree_5_5;terms={len(J.terms())};"
        f"positive={positive_terms};negative={negative_terms}"
    )
    print("selected_coefficients=273888,-5839200,5559264,-82033920")
    print("endpoint_polynomial_sha256=" + polynomial_digest(J))
    print(
        "cell_1_1="
        "g(30,37,46);"
        "rL=73097746/75;rR=515413080/529;"
        "J=610878432"
    )
    print(
        "cell_1_1="
        "I(141000,169320);"
        "TP3(36467293512,86548240800)"
    )
    print(
        "cell_1_2="
        "g(30,61,126);"
        "rL=142883488/75;rR=120168140/63;"
        "J=-33115086144"
    )
    print(
        "cell_1_2="
        "I(957960,2436120);"
        "TP3(36467293512,2255452579800)"
    )
    print("quotient_endpoint_derivatives=2rL,2rR")
    print("missing_predicate=simultaneous_cubic_nullity")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
