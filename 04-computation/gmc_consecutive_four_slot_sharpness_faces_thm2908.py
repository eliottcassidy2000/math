#!/usr/bin/env python3
"""Exact coordinate-face and sharpness companion for THM-2908.

The two nonconsecutive three-slot faces of
``span{s^n,s^(n+1),s^(n+2),s^(n+3)}`` are reduced by the first
factorial moment to a quadratic and a cubic in one projective coordinate.
Their resultants have coefficient-positive nonlinear factors.  Together
with THM-2812 on the two consecutive faces, this makes every
Krull witness for the first three moments full-support.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


n, a, b, c, t = sp.symbols("n a b c t")


def normalized_moment(offsets: tuple[int, int, int], power: int) -> sp.Expr:
    """Return L(H^power)/(power*n)! for H on the given offset support."""

    variables = (a, b, c)
    total = sp.Integer(0)
    for i in range(power + 1):
        for j in range(power - i + 1):
            k = power - i - j
            counts = (i, j, k)
            offset = sum(
                count * support_offset
                for count, support_offset in zip(counts, offsets)
            )
            multinomial = (
                sp.factorial(power)
                / (
                    sp.factorial(i)
                    * sp.factorial(j)
                    * sp.factorial(k)
                )
            )
            total += (
                multinomial
                * sp.prod(
                    variable**count
                    for variable, count in zip(variables, counts)
                )
                * sp.rf(power * n + 1, offset)
            )
    return sp.expand(total)


expected = {
    (0, 1, 3): (
        16 * (n + 1) ** 7 * (n + 2) * (n + 3) ** 6,
        sp.Poly(
            186624 * n**10
            + 5443200 * n**9
            + 42920064 * n**8
            + 165194816 * n**7
            + 370368960 * n**6
            + 519725672 * n**5
            + 467821669 * n**4
            + 268012211 * n**3
            + 93650626 * n**2
            + 18023940 * n
            + 1453192,
            n,
        ),
        (2, 6, 15, 3, 9, 28),
    ),
    (0, 2, 3): (
        4 * (n + 1) ** 5 * (n + 2) ** 7 * (n + 3) ** 6,
        sp.Poly(
            28344976 * n**12
            + 855385784 * n**11
            + 7207831697 * n**10
            + 30256509235 * n**9
            + 76068084600 * n**8
            + 124601091222 * n**7
            + 138707331317 * n**6
            + 106863625003 * n**5
            + 56960915338 * n**4
            + 20591978516 * n**3
            + 4810575912 * n**2
            + 653870592 * n
            + 39189888,
            n,
        ),
        (2, 6, 18, 3, 9, 34),
    ),
}

records: list[str] = []
for offsets, (linear_factors, nonlinear, profile) in expected.items():
    first = normalized_moment(offsets, 1)
    first_a = sp.Poly(first, a).coeff_monomial(a)
    first_b = sp.Poly(first, b).coeff_monomial(b)
    first_c = sp.Poly(first, c).coeff_monomial(c)
    require(first_a == 1, f"first-moment a coefficient changed: {offsets}")

    a_value = -(first_b * t + first_c)
    quadratic = sp.Poly(
        normalized_moment(offsets, 2).subs({a: a_value, b: t, c: 1}),
        t,
        n,
    )
    cubic = sp.Poly(
        normalized_moment(offsets, 3).subs({a: a_value, b: t, c: 1}),
        t,
        n,
    )
    actual_profile = (
        quadratic.degree(t),
        quadratic.degree(n),
        len(quadratic.terms()),
        cubic.degree(t),
        cubic.degree(n),
        len(cubic.terms()),
    )
    require(actual_profile == profile, f"moment profile changed: {offsets}")

    resultant = sp.Poly(
        sp.resultant(quadratic.as_expr(), cubic.as_expr(), t),
        n,
    )
    expected_resultant = sp.Poly(
        sp.expand(linear_factors * nonlinear.as_expr()),
        n,
    )
    require(
        resultant == expected_resultant,
        f"three-slot face resultant changed: {offsets}",
    )
    require(
        all(coefficient > 0 for coefficient in nonlinear.all_coeffs()),
        f"nonlinear face factor lost coefficient positivity: {offsets}",
    )
    records.append(
        f"offsets={offsets};profiles={profile};"
        f"nonlinear_degree={nonlinear.degree()};"
        f"nonlinear_constant={nonlinear.eval(0)};"
        "all_coefficients_positive=true"
    )

print("THM-2908 SHARPNESS FACE CERTIFICATES")
for record in records:
    print(record)
print("coordinate_faces=two_consecutive_by_THM2812+two_exact_resultants")
print("krull_witness=full_support;factorial_depth=4;gaussian_depth=8")
print("all_exact_checks=PASS")
