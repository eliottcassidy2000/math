#!/usr/bin/env python3
"""Exact controls for THM-3465's cyclic-character Keller rigidity."""

import ast
import hashlib
from math import factorial
from pathlib import Path

import sympy as sp


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


z, w = sp.symbols("z w")
A, B, C, D, E = sp.symbols("A B C D E")
a, b, c, d, e = sp.symbols("a b c d e")

g4 = A * w + B * z**2 + C * z * w**2 + D * z**3 * w + E * w**4
h4 = a * z + b * w**2 + c * w * z**2 + d * w**3 * z + e * z**4
bracket4 = sp.Poly(
    sp.expand(sp.diff(g4, z) * sp.diff(h4, w) - sp.diff(g4, w) * sp.diff(h4, z)),
    z,
    w,
)

quartic_expected = (
    ((0, 0), -A * a),
    ((1, 1), -2 * (A * c - 2 * B * b + C * a)),
    ((0, 3), -A * d + 2 * C * b - 4 * E * a),
    ((3, 0), -4 * A * e + 2 * B * c - D * a),
    ((2, 2), 3 * (2 * B * d - C * c + 2 * D * b)),
    ((1, 4), C * d - 8 * E * c),
    ((4, 1), -8 * C * e + D * c),
    ((0, 6), -4 * E * d),
    ((3, 3), 8 * (D * d - 2 * E * e)),
    ((6, 0), -4 * D * e),
)
for (iz, iw), expected in quartic_expected:
    actual = bracket4.coeff_monomial(z**iz * w**iw)
    require(sp.expand(actual - expected) == 0, ("quartic coefficient", iz, iw, actual))
require(len(bracket4.terms()) == len(quartic_expected), ("unexpected quartic support", bracket4.terms()))


def jacobian(f, h, u, v):
    return sp.expand(sp.diff(f, u) * sp.diff(h, v) - sp.diff(f, v) * sp.diff(h, u))


# Equal-degree homogeneous bracket formula on the chart z/w.
t = sp.symbols("t")
for degree in range(1, 8):
    F = sum(sp.Symbol(f"f{degree}_{j}") * t**j for j in range(degree + 1))
    H = sum(sp.Symbol(f"h{degree}_{j}") * t**j for j in range(degree + 1))
    f = sp.expand(w**degree * F.subs(t, z / w))
    h = sp.expand(w**degree * H.subs(t, z / w))
    expected = sp.expand(
        degree
        * w ** (2 * degree - 2)
        * (sp.diff(F, t) * H - F * sp.diff(H, t)).subs(t, z / w)
    )
    require(sp.expand(jacobian(f, h, z, w) - expected) == 0, ("binary-form lemma", degree))


# Opposite nonreal characters stay disjoint on every homogeneous layer for
# every checked finite rotation.  Exponents are measured relative to the
# source character xi in rho(z,w)=(xi*z,xi^(-1)*w).
finite_character_checks = 0
for order in range(3, 17):
    for character in range(1, order):
        if (2 * character) % order == 0:
            continue
        for degree in range(1, 25):
            character_basis = tuple(
                (iz, degree - iz)
                for iz in range(degree + 1)
                if (iz - (degree - iz)) % order == character
            )
            inverse_basis = tuple(
                (iz, degree - iz)
                for iz in range(degree + 1)
                if (iz - (degree - iz)) % order == (-character) % order
            )
            require(
                set(character_basis).isdisjoint(inverse_basis),
                ("finite character overlap", order, character, degree),
            )
            star_basis = tuple((iw, iz) for iz, iw in character_basis)
            require(
                set(star_basis) == set(inverse_basis),
                ("finite star character", order, character, degree),
            )
            finite_character_checks += 1


# Pin the THM-3310 convention separately.
for degree in range(1, 25):
    omega_basis = tuple(
        (iz, degree - iz)
        for iz in range(degree + 1)
        if (2 * iz + degree - iz) % 3 == 1
    )
    inverse_basis = tuple(
        (iw, degree - iw)
        for iw in range(degree + 1)
        if (2 * iw + degree - iw) % 3 == 2
    )
    require(set(omega_basis).isdisjoint(inverse_basis), ("character overlap", degree))
    star_basis = tuple((iw, iz) for iz, iw in omega_basis)
    require(set(star_basis) == set(inverse_basis), ("star character", degree, star_basis, inverse_basis))


# Exact triangle coordinates and normalized third moment.
x, y = sp.symbols("x y")
sqrt3 = sp.sqrt(3)
I = sp.I
omega = (-1 + sqrt3 * I) / 2
s1, s2, s3 = x, y, 1 - x - y
zxy = sp.expand(s1 + omega * s2 + omega**2 * s3)
wxy = sp.expand(s1 + omega**2 * s2 + omega * s3)
coordinate_jacobian = sp.simplify(jacobian(zxy, wxy, x, y))
require(coordinate_jacobian == -3 * sqrt3 * I, ("coordinate Jacobian", coordinate_jacobian))


def triangle_average(poly):
    total = 0
    for (ix, iy), coefficient in sp.Poly(sp.expand(poly), x, y).terms():
        integral = sp.Rational(factorial(ix) * factorial(iy), factorial(ix + iy + 2))
        total += 2 * coefficient * integral
    return sp.simplify(total)


w3_average = triangle_average(wxy**3)
require(w3_average == sp.Rational(1, 10), ("w^3 average", w3_average))

A0 = 2 + 3 * I
h_linear = sp.conjugate(A0) * zxy
g_linear = A0 * wxy
p_linear = sp.expand((g_linear + h_linear) / 2)
q_linear = sp.expand((g_linear - h_linear) / (2 * I))
real_jacobian = sp.simplify(jacobian(p_linear, q_linear, x, y))
require(real_jacobian == -sp.Rational(39, 2) * sqrt3, ("real Jacobian", real_jacobian))
require(
    sp.simplify(triangle_average(g_linear**3) - A0**3 / sp.Integer(10)) == 0,
    "linear third moment",
)


# Load-bearing hostiles.
formal_B = sp.symbols("formal_B")
require(jacobian(w + formal_B * z**2, z, z, w) == -1, "independent-mate hostile")
P_c2, Q_c2 = x, y + x**3
require(jacobian(P_c2, Q_c2, x, y) == 1, "C2 shear Jacobian")
require(
    sp.expand(P_c2.subs({x: -x, y: -y}) + P_c2) == 0
    and sp.expand(Q_c2.subs({x: -x, y: -y}) + Q_c2) == 0,
    "C2 odd shear",
)


def factorial_functional(poly, variable):
    ans = 0
    for (degree,), coefficient in sp.Poly(sp.expand(poly), variable).terms():
        ans += coefficient * factorial(degree)
    return sp.expand(ans)


u = sp.symbols("u", real=True)
f_hostile = u**2 + (-4 + 2 * I) * u + (2 - 2 * I)
f_star = sp.conjugate(f_hostile)
hostile_values = (
    factorial_functional(f_hostile, u),
    factorial_functional(f_hostile**2, u),
    factorial_functional(f_hostile**3, u),
    factorial_functional(f_hostile * f_star, u),
)
require(hostile_values == (0, 0, 32 + 80 * I, 8), ("twisted-kernel hostile", hostile_values))

source = Path(__file__).read_text(encoding="utf-8")
require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))), "assert node")

print("THM-3465 NONREAL CYCLIC-CHARACTER KELLER RIGIDITY EXACT CONTROL")
print("quartic_bracket_support=%d coefficient_table=PASS" % len(quartic_expected))
print("binary_form_bracket_formula=PASS degrees=1..7")
print("finite_rotation_character_star_disjoint=PASS orders=3..16 "
      "degrees=1..24 checks=%d" % finite_character_checks)
print("C3_character_convention=PASS degrees=1..24")
print("Jac_xy(z,w)=%s" % coordinate_jacobian)
print("normalized_triangle_average_w3=%s" % w3_average)
print("linear_control_A=%s real_J=%s third_moment=%s" % (
    A0, real_jacobian, triangle_average(g_linear**3)
))
print("hostiles=independent_mate_bracket_-1; C2_odd_shear_J_1")
print("twisted_kernel_hostile=%s" % (hostile_values,))
print("semantic_sha256=%s" % hashlib.sha256(
    repr((
        quartic_expected,
        finite_character_checks,
        coordinate_jacobian,
        w3_average,
        hostile_values,
    )).encode("utf-8")
).hexdigest())
print("STATUS=PASS")
