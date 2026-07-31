#!/usr/bin/env python3
"""Exact companion for THM-2917.

For each translation-normalized three-subset of ``{0,1,2,3,4}``, eliminate
the first factorial moment and compute the resultant of the remaining
quadratic and cubic.  Every irreducible factor has strictly positive
coefficients.
"""

from __future__ import annotations

from hashlib import sha256

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


n, a, b, c, t = sp.symbols("n a b c t")
SUPPORTS = (
    (0, 1, 2),
    (0, 1, 3),
    (0, 2, 3),
    (0, 1, 4),
    (0, 2, 4),
    (0, 3, 4),
)

EXPECTED_PROFILES = {
    (0, 1, 2): (2, 4, 12, 3, 6, 22, 18, 19),
    (0, 1, 3): (2, 6, 15, 3, 9, 28, 24, 25),
    (0, 2, 3): (2, 6, 18, 3, 9, 34, 30, 31),
    (0, 1, 4): (2, 8, 18, 3, 12, 34, 30, 31),
    (0, 2, 4): (2, 8, 21, 3, 12, 40, 36, 37),
    (0, 3, 4): (2, 8, 24, 3, 12, 46, 42, 43),
}

EXPECTED_NONLINEAR = {
    (0, 1, 2): (
        5,
        29,
        "132895f955d191bc3822c732bbc37b43331b43f8fd9bcf677fa9e3282cfd146c",
    ),
    (0, 1, 3): (
        10,
        1453192,
        "6fa8012c7b9d57e35886cbbecc477a43f7ffe42168c3798b24b55f8da35c3685",
    ),
    (0, 2, 3): (
        12,
        39189888,
        "cd6fdf453309116c9620eff9eadd17b51b25bbbccd60f597ff1f82b99c33fb12",
    ),
    (0, 1, 4): (
        14,
        765149683104,
        "cf35e711e297e259fdf0358655e349b4e1ed0443fc09a876b513d94b34d9b43f",
    ),
    (0, 2, 4): (
        17,
        16074139599576,
        "98b87bf47164bb6cb73c5f141ec76be103e7642fbcb6099677c508974cc03727",
    ),
    (0, 3, 4): (
        19,
        744201921369600,
        "820ecaabdbf14eb2b21e7cb56d5235cc167bf49d1e44361ffd93296d8cb73db1",
    ),
}


def normalized_moment(
    offsets: tuple[int, int, int], power: int
) -> sp.Expr:
    """Return ``L(H^power)/(power*n)!`` on one normalized support."""

    answer = sp.Integer(0)
    for i in range(power + 1):
        for j in range(power - i + 1):
            k = power - i - j
            offset = i * offsets[0] + j * offsets[1] + k * offsets[2]
            answer += (
                sp.factorial(power)
                / (
                    sp.factorial(i)
                    * sp.factorial(j)
                    * sp.factorial(k)
                )
                * a**i
                * b**j
                * c**k
                * sp.rf(power * n + 1, offset)
            )
    return sp.expand(answer)


def coefficient_digest(polynomial: sp.Poly) -> str:
    payload = ",".join(str(value) for value in polynomial.all_coeffs()) + "\n"
    return sha256(payload.encode()).hexdigest()


records: list[str] = []
for offsets in SUPPORTS:
    first = sp.Poly(normalized_moment(offsets, 1), a, b, c)
    first_a = first.coeff_monomial(a)
    first_b = first.coeff_monomial(b)
    first_c = first.coeff_monomial(c)
    require(first_a == 1, f"first-moment normalization changed: {offsets}")
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
    resultant = sp.Poly(
        sp.resultant(quadratic.as_expr(), cubic.as_expr(), t),
        n,
    )
    profile = (
        quadratic.degree(t),
        quadratic.degree(n),
        len(quadratic.terms()),
        cubic.degree(t),
        cubic.degree(n),
        len(cubic.terms()),
        resultant.degree(),
        len(resultant.terms()),
    )
    require(
        profile == EXPECTED_PROFILES[offsets],
        f"resultant profile changed: {offsets}",
    )

    unit, factors = sp.factor_list(resultant.as_expr(), n)
    normalized_factors: list[tuple[sp.Poly, int]] = []
    for factor, exponent in factors:
        factor_poly = sp.Poly(factor, n)
        if factor_poly.LC() < 0:
            factor_poly = -factor_poly
        require(
            all(value > 0 for value in factor_poly.all_coeffs()),
            f"factor lost coefficient positivity: {offsets}",
        )
        normalized_factors.append((factor_poly, exponent))

    nonlinear = max(normalized_factors, key=lambda item: item[0].degree())[0]
    expected_degree, expected_constant, expected_digest = (
        EXPECTED_NONLINEAR[offsets]
    )
    require(
        (
            nonlinear.degree(),
            nonlinear.eval(0),
            coefficient_digest(nonlinear),
        )
        == (expected_degree, expected_constant, expected_digest),
        f"nonlinear factor changed: {offsets}",
    )
    factor_profile = ",".join(
        f"{factor.degree()}^{exponent}"
        for factor, exponent in sorted(
            normalized_factors,
            key=lambda item: (item[0].degree(), item[0].eval(0)),
        )
    )
    records.append(
        f"support={offsets};moment_profile={profile[:6]};"
        f"resultant={profile[6]},{profile[7]};unit={unit};"
        f"factor_profile={factor_profile};"
        f"nonlinear={expected_degree},{expected_constant},{expected_digest}"
    )

print("THM-2917 ALL THREE-SLOT DIAMETER-FOUR DETECTION")
for record in records:
    print(record)
print("all_resultant_factors_have_positive_coefficients=true")
print("two_term_boundary=positive_L2_after_L1")
print("scope=first-window SFC(3), every translated three-subset of five consecutive exponents")
print("all_exact_checks=PASS")
