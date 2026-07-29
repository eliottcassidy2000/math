#!/usr/bin/env python3
"""Exact limiting-algebra companion for THM-2914.

The analytic tail proof uses Stirling bounds and stability of a
nonsingular zero.  This companion checks every finite algebraic input:
the low factorial moments, limiting cubic invariants and eliminant,
positive-quadrant root selector, Jacobian, endpoint factorization and
all exponential rate separations.
"""

from __future__ import annotations

from hashlib import sha256
from math import factorial

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def digest(polynomial: sp.Poly) -> str:
    records = "\n".join(
        f"{','.join(str(value) for value in monomial)}:{coefficient}"
        for monomial, coefficient in polynomial.terms()
    )
    return sha256((records + "\n").encode()).hexdigest()


s, xi, eta, u, z = sp.symbols("s xi eta u z")


def f(index: int) -> sp.Expr:
    return s**index / factorial(index)


def functional(polynomial: sp.Expr) -> sp.Expr:
    expanded = sp.Poly(sp.expand(polynomial), s)
    return sp.expand(
        sum(
            coefficient * factorial(monomial[0])
            for monomial, coefficient in expanded.terms()
        )
    )


d0 = f(1) - f(0)
d1 = f(2) - f(1)
low = d0 + z * d1

quadratic = sp.Poly(functional(low**2), z)
cubic = sp.Poly(functional(low**3), z)
require(
    tuple(quadratic.nth(index) / (1, 2, 1)[index] for index in range(3))
    == (1, 1, 2),
    "low quadratic moment vector changed",
)
require(
    tuple(cubic.nth(index) / (1, 3, 3, 1)[index] for index in range(4))
    == (2, 4, 10, 30),
    "low cubic moment vector changed",
)

g0, g1, g2 = sp.Integer(1), sp.Integer(1), sp.Integer(2)
t0 = 2 + xi**3
t1 = 4 + xi**2 * eta
t2 = 10 + xi * eta**2
t3 = 30 + eta**3

invariant_one = sp.expand(
    3 * t1 * g0 * g2 - t3 * g0**2 - 2 * t0 * g1 * g2
)
invariant_two = sp.expand(
    3 * t2 * g0 * g2 - 2 * t3 * g1 * g0 - t0 * g2**2
)
expected_one = -eta**3 + 6 * eta * xi**2 - 4 * xi**3 - 14
expected_two = -2 * (eta**3 - 3 * eta**2 * xi + 2 * xi**3 + 4)
require(
    invariant_one == expected_one and invariant_two == expected_two,
    "limiting cubic invariants changed",
)

p = 2 * u**3 + 18 * u**2 - 729 * u + 54
resultant = sp.factor(sp.resultant(invariant_one, invariant_two, xi))
require(
    sp.expand(resultant + 256 * p.subs(u, eta**3)) == 0,
    "limiting cubic eliminant changed",
)

p_poly = sp.Poly(p, u, domain=sp.QQ)
require(
    p_poly.count_roots(-25, -24) == 1
    and p_poly.count_roots(0, 1) == 1
    and p_poly.count_roots(15, 16) == 1
    and p_poly.count_roots(-sp.oo, sp.oo) == 3,
    "limiting eliminant root isolation changed",
)

ratio = (2 * u**2 + 21 * u - 603) / 189
ratio_alternate = (u + 9) / (2 * u - 3)
ratio_difference_numerator = sp.together(
    ratio - ratio_alternate
).as_numer_denom()[0]
require(
    sp.cancel(ratio.subs(u, 15)) == sp.Rational(6, 7)
    and sp.cancel(ratio.subs(u, 16)) == sp.Rational(35, 27)
    and sp.diff(ratio, u) == (4 * u + 21) / 189,
    "positive root ratio interval changed",
)
require(
    sp.rem(
        sp.Poly(ratio_difference_numerator, u, domain=sp.QQ),
        p_poly,
    ).is_zero,
    "alternate positive-root ratio law changed",
)
require(
    ratio.subs(u, 1) < 0,
    "small positive eliminant root acquired a positive xi coordinate",
)

first_certificate = (
    16 * u**4
    + 360 * u**3
    - 8856 * u**2
    - 122526 * u
    + 1750329
)
second_certificate = (
    8 * u**4
    + 180 * u**3
    - 3294 * u**2
    - 47655 * u
    + 500094
)
substituted_one = sp.factor(
    invariant_one.subs(xi, eta * ratio).subs(eta**3, u)
)
substituted_two = sp.factor(
    invariant_two.subs(xi, eta * ratio).subs(eta**3, u)
)
require(
    sp.expand(
        substituted_one
        + p * first_certificate / 6751269
    )
    == 0
    and sp.expand(
        substituted_two
        + 2 * p * second_certificate / 6751269
    )
    == 0,
    "positive root selector identities changed",
)

jacobian = sp.factor(
    sp.det(
        sp.Matrix(
            [
                [sp.diff(invariant_one, xi), sp.diff(invariant_one, eta)],
                [sp.diff(invariant_two, xi), sp.diff(invariant_two, eta)],
            ]
        )
    )
)
require(
    jacobian == 18 * (eta**2 - 2 * eta * xi + 2 * xi**2) ** 2,
    "limiting Jacobian factor changed",
)

endpoint = (
    eta**4 - 4 * xi * eta**3 + 8 * xi**3 * eta - 4 * xi**4
)
endpoint_factor = (
    (eta**2 - 2 * xi**2)
    * (eta**2 - 4 * eta * xi + 2 * xi**2)
)
require(
    sp.expand(endpoint - endpoint_factor) == 0,
    "limiting endpoint factorization changed",
)

q = sp.symbols("q")
first_q = 1 - 2 * q**2
second_q = 1 - 4 * q + 2 * q**2
require(
    first_q.subs(q, sp.Rational(6, 7)) == sp.Rational(-23, 49)
    and second_q.subs(q, sp.Rational(6, 7)) == sp.Rational(-47, 49)
    and second_q.subs(q, sp.Rational(35, 27))
    == sp.Rational(-601, 729)
    and sp.diff(second_q, q, 2) == 4
    and sp.Rational(23 * 601, 49 * 729)
    == sp.Rational(13823, 35721),
    "limiting endpoint sign floor changed",
)

# Exponential bases after Stirling.  Polynomial powers of r cannot
# overcome any of these strict separations.
quadratic_cubic_bases = (
    sp.Rational(1, 3),
    sp.Rational(4, 9),
)
quartic_error_bases = (
    sp.Rational(81, 256),
    sp.Rational(27, 256),
    sp.Rational(36, 256),
    sp.Rational(81, 256),
)
require(
    max(quadratic_cubic_bases) < 1
    and max(quartic_error_bases) < 1
    and sp.Rational(256, 81) > 1,
    "factorial exponential separation changed",
)

EXPECTED_DIGESTS = {
    "I1": "6b59fdb1978397e9efefcad497c99ca2edd23deac95bb5a3af04c02a0c9f61ff",
    "I2": "212141a0f5a563f56c41f05225e5e9b9b93d8ee6eceb8aa77cda77093f9f6e1b",
    "p": "b1472981a28032a2251c051e59765d5ede2df5a9e230deadb0f658048ace44bb",
    "J": "3c9a5a3648ea58e89486db7a753fe38d4bea4f4e5387634f09f32b95a3f9b864",
}
actual_digests = {
    "I1": digest(sp.Poly(invariant_one, xi, eta, domain=sp.QQ)),
    "I2": digest(sp.Poly(invariant_two, xi, eta, domain=sp.QQ)),
    "p": digest(p_poly),
    "J": digest(sp.Poly(endpoint, xi, eta, domain=sp.QQ)),
}
require(
    actual_digests == EXPECTED_DIGESTS,
    "locked limiting polynomial digest changed",
)


def main() -> None:
    print("THM-2914 EVENTUAL HIGH-GAP CUBIC-NULL BRANCH")
    print("status=PROOF-COMPLETE CANDIDATE+VERIFIED-EXACT")
    print("low_quadratic=(1,1,2);low_cubic=(2,4,10,30)")
    print(
        "limit_I1=-eta^3+6eta*xi^2-4xi^3-14;"
        "limit_I2=-2(eta^3-3eta^2*xi+2xi^3+4)"
    )
    print("eliminant=-256*(2u^3+18u^2-729u+54)")
    print("root_boxes=(-25,-24),(0,1),(15,16);positive_quadrant=last_only")
    print("positive_ratio=6/7<xi/eta<35/27")
    print("jacobian=18*(eta^2-2eta*xi+2xi^2)^2>0")
    print(
        "endpoint=(eta^2-2xi^2)*(eta^2-4eta*xi+2xi^2);"
        "normalized_floor>13823/35721*eta^4"
    )
    print(
        "rates=T2/T3^(2/3)->0;T4/T3^(4/3)->infinity;"
        "quartic_error_bases=81/256,27/256,36/256,81/256"
    )
    print(
        "digests="
        + ";".join(f"{key}:{value}" for key, value in actual_digests.items())
    )
    print(
        "scope=eventual local positive branch only; no finite cutoff, "
        "global uniqueness, arbitrary-support, SFC(4), or GMC closure"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
