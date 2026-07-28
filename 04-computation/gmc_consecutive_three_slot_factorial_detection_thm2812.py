from __future__ import annotations

from math import factorial

import sympy as sp


s, t, n = sp.symbols("s t n", integer=True, nonnegative=True)
x, y, z = sp.symbols("x y z")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def rising(base: sp.Expr, length: int) -> sp.Expr:
    return sp.prod(base + j for j in range(1, length + 1))


def normalized_factorial_moment(power: int) -> sp.Expr:
    """Return L((s^n(x+y*s+z*s^2))^power)/(power*n)!."""
    poly = sp.Poly(sp.expand((x + y * s + z * s**2) ** power), s)
    return sp.expand(
        sum(
            coefficient * rising(power * n, degree[0])
            for degree, coefficient in poly.terms()
        )
    )


def factorial_functional(poly: sp.Expr) -> sp.Expr:
    expanded = sp.Poly(sp.expand(poly), s)
    return sp.simplify(
        sum(
            coefficient * factorial(degree[0])
            for degree, coefficient in expanded.terms()
        )
    )


F1 = normalized_factorial_moment(1)
F2 = normalized_factorial_moment(2)
F3 = normalized_factorial_moment(3)

require(
    sp.factor(F1 - (x + (n + 1) * y + (n + 1) * (n + 2) * z)) == 0,
    "first normalized moment identity failed",
)

z_from_F1 = -(1 + (n + 1) * t) / ((n + 1) * (n + 2))
second_reduced = sp.factor(F2.subs({x: 1, y: t, z: z_from_F1}))
third_reduced = sp.factor(F3.subs({x: 1, y: t, z: z_from_F1}))

Q = (
    2 * (n + 1) ** 2 * (2 * n + 1) * t**2
    + 6 * (n + 1) * (2 * n + 1) * t
    + 9 * n
    + 5
)
require(
    sp.factor(second_reduced - Q / (n + 1)) == 0,
    "second-moment conic identity failed",
)

disc_Q = sp.factor(sp.discriminant(Q, t))
require(
    disc_Q == -4 * (n + 1) ** 2 * (2 * n + 1),
    "second-moment discriminant identity failed",
)

A = 36 * n**4 + 57 * n**3 + 15 * n**2 - 9 * n - 3
B = 52 * n**3 + 49 * n**2 + 8 * n - 1
expected_remainder = -(A * t + B) / ((n + 1) ** 2 * (2 * n + 1))

third_num, third_den = sp.together(third_reduced).as_numer_denom()
second_num, _ = sp.together(second_reduced).as_numer_denom()
raw_remainder = (
    sp.Poly(third_num, t).rem(sp.Poly(second_num, t)).as_expr() / third_den
)
require(
    sp.factor(raw_remainder - expected_remainder) == 0,
    "third-moment Euclidean remainder identity failed",
)

require(
    (A.subs(n, 0), B.subs(n, 0)) == (-3, -1),
    "n=0 remainder boundary failed",
)
for n0 in range(1, 201):
    require(int(A.subs(n, n0)) > 0, f"A positivity control failed at n={n0}")
    require(int(B.subs(n, n0)) > 0, f"B positivity control failed at n={n0}")

# A sharp full-support witness at n=1.
sqrt3 = sp.sqrt(3)
t_star = -sp.Rational(3, 4) + sp.I * sqrt3 / 12
z_star = sp.Rational(1, 12) - sp.I * sqrt3 / 36
H_star = s * (1 + t_star * s + z_star * s**2)
sharp_moments = tuple(
    sp.factor(factorial_functional(H_star**j)) for j in range(1, 4)
)
require(sharp_moments[:2] == (0, 0), "sharp witness first two moments failed")
require(
    sp.simplify(sharp_moments[2] - (-18 - 4 * sp.I * sqrt3)) == 0,
    "sharp witness third moment failed",
)
gaussian_m6 = sp.binomial(6, 3) * sharp_moments[2]
require(
    sp.simplify(gaussian_m6 - (-360 - 80 * sp.I * sqrt3)) == 0,
    "Gaussian lift sixth moment failed",
)

print("THM-2812 CONSECUTIVE THREE-SLOT FACTORIAL DETECTION")
print("symbolic_moments=PASS")
print(f"Q_discriminant={disc_Q}")
print(f"cubic_remainder_A={A}")
print(f"cubic_remainder_B={B}")
print("positive_integer_control=n=1..200")
print(f"sharp_L1_L2_L3={sharp_moments}")
print(f"sharp_gaussian_M6={gaussian_m6}")
print(
    "scope=all n>=0 factorial theorem; polynomial Gaussian lift for n>=1; "
    "arbitrary three-slot supports and HYP-8765 remain open"
)
print("ALL EXACT CHECKS PASSED")
