#!/usr/bin/env python3
"""Exact controls for THM-3466's factorial/Keller Stokes identities."""

import ast
import hashlib
import itertools
from math import factorial
from pathlib import Path

import sympy as sp


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def factorial_functional(poly, variables):
    poly = sp.expand(poly)
    if not variables:
        return poly
    total = 0
    for exponents, coefficient in sp.Poly(poly, *variables).terms():
        weight = 1
        for exponent in exponents:
            weight *= factorial(exponent)
        total += coefficient * weight
    return sp.expand(total)


def face_functional(poly, variables, face_indices):
    substitutions = {variables[index]: 0 for index in face_indices}
    remaining = tuple(variable for index, variable in enumerate(variables) if index not in face_indices)
    return factorial_functional(sp.expand(poly.subs(substitutions)), remaining)


X, Y, Z = sp.symbols("X Y Z", real=True)
variables3 = (X, Y, Z)
h3 = (
    7
    + 2 * X
    - 3 * Y**2
    + 5 * X * Y * Z
    - 4 * X**3 * Z**2
    + 11 * Y * Z**4
)

face_checks = 0
for index, variable in enumerate(variables3):
    lhs = factorial_functional(sp.diff(h3, variable), variables3)
    rhs = factorial_functional(h3, variables3) - face_functional(h3, variables3, (index,))
    require(sp.expand(lhs - rhs) == 0, ("codimension one", index, lhs, rhs))
    face_checks += 1

for size in range(4):
    for indices in itertools.combinations(range(3), size):
        descended = h3
        for index in indices:
            descended = sp.expand(descended - sp.diff(descended, variables3[index]))
        lhs = factorial_functional(descended, variables3)
        rhs = face_functional(h3, variables3, indices)
        require(sp.expand(lhs - rhs) == 0, ("iterated face", indices, lhs, rhs))
        face_checks += 1

f3 = 1 + X - 2 * Y + X * Y + 3 * Z**2
for power in range(1, 6):
    for index, variable in enumerate(variables3):
        boundary = face_functional(f3**power, variables3, (index,))
        bulk_minus_lowering = factorial_functional(f3**power, variables3) - power * factorial_functional(
            f3 ** (power - 1) * sp.diff(f3, variable), variables3
        )
        require(sp.expand(boundary - bulk_minus_lowering) == 0, ("power face", power, index))
        face_checks += 1


def bracket(P, H, x, y):
    return sp.expand(sp.diff(P, x) * sp.diff(H, y) - sp.diff(P, y) * sp.diff(H, x))


def edge_x(poly):
    return factorial_functional(sp.expand(poly.subs(X, 0)), (Y,))


def edge_y(poly):
    return factorial_functional(sp.expand(poly.subs(Y, 0)), (X,))


keller_flux_checks = 0
for P, Q in ((X, Y), (X, Y + X**2)):
    c = bracket(P, Q, X, Y)
    require(c == 1, ("Keller control", P, Q, c))
    for aa in range(4):
        for bb in range(4):
            H = sp.expand(P**aa * Q**bb)
            general_p = (
                factorial_functional(bracket(P, H, X, Y), (X, Y))
                - factorial_functional(H * (sp.diff(P, X) - sp.diff(P, Y)), (X, Y))
                - edge_x(H * sp.diff(P, Y))
                + edge_y(H * sp.diff(P, X))
            )
            require(sp.expand(general_p) == 0, ("P flux", P, Q, aa, bb, general_p))
            expected_p = 0 if bb == 0 else bb * c * factorial_functional(P**aa * Q ** (bb - 1), (X, Y))
            require(
                sp.expand(factorial_functional(bracket(P, H, X, Y), (X, Y)) - expected_p) == 0,
                ("P monomial recurrence", aa, bb),
            )

            general_q = (
                factorial_functional(bracket(Q, H, X, Y), (X, Y))
                - factorial_functional(H * (sp.diff(Q, X) - sp.diff(Q, Y)), (X, Y))
                - edge_x(H * sp.diff(Q, Y))
                + edge_y(H * sp.diff(Q, X))
            )
            require(sp.expand(general_q) == 0, ("Q flux", P, Q, aa, bb, general_q))
            expected_q = 0 if aa == 0 else -aa * c * factorial_functional(P ** (aa - 1) * Q**bb, (X, Y))
            require(
                sp.expand(factorial_functional(bracket(Q, H, X, Y), (X, Y)) - expected_q) == 0,
                ("Q monomial recurrence", aa, bb),
            )
            keller_flux_checks += 4


t = sp.symbols("t", real=True)
I = sp.I


def triangle_integral(poly):
    total = 0
    for (ix, iy), coefficient in sp.Poly(sp.expand(poly), X, Y).terms():
        total += coefficient * sp.Rational(factorial(ix) * factorial(iy), factorial(ix + iy + 2))
    return sp.simplify(total)


triangle_edges = (
    (t, 0),
    (1 - t, t),
    (0, 1 - t),
)


def boundary_integral(form_coefficient, differential_poly):
    total = 0
    for edge_x_value, edge_y_value in triangle_edges:
        substitutions = {X: edge_x_value, Y: edge_y_value}
        coefficient = sp.expand(form_coefficient.subs(substitutions))
        differential = sp.diff(sp.expand(differential_poly.subs(substitutions)), t)
        total += sp.integrate(sp.expand(coefficient * differential), (t, 0, 1))
    return sp.simplify(total)


sqrt3 = sp.sqrt(3)
omega = (-1 + sqrt3 * I) / 2
s1, s2, s3 = X, Y, 1 - X - Y
z_triangle = sp.expand(s1 + omega * s2 + omega**2 * s3)
w_triangle = sp.expand(s1 + omega**2 * s2 + omega * s3)

triangle_controls = (
    ("identity", X, Y),
    ("triangular", X, Y + X**2),
    ("equilateral", sp.expand((z_triangle + w_triangle) / 2), sp.expand((z_triangle - w_triangle) / (2 * I))),
)

triangle_checks = 0
control_rows = []
for name, p, q in triangle_controls:
    g = sp.expand(p + I * q)
    gbar = sp.expand(p - I * q)
    c = sp.simplify(bracket(p, q, X, Y))
    require(c != 0 and not (c.has(X) or c.has(Y)), ("constant J control", name, c))
    moments = tuple(triangle_integral(g**power) for power in range(4))
    for n in range(1, 5):
        dp_response = boundary_integral(g**n, p)
        dq_response = boundary_integral(g**n, q)
        dbar_response = boundary_integral(g**n, gbar)
        expected_dp = -n * I * c * moments[n - 1]
        expected_dq = n * c * moments[n - 1]
        expected_dbar = -2 * n * I * c * moments[n - 1]
        require(sp.simplify(dp_response - expected_dp) == 0, ("dp response", name, n))
        require(sp.simplify(dq_response - expected_dq) == 0, ("dq response", name, n))
        require(sp.simplify(dbar_response - expected_dbar) == 0, ("dbar response", name, n))
        require(sp.simplify(boundary_integral(g**n, g)) == 0, ("exact dg response", name, n))
        triangle_checks += 4
    control_rows.append((name, c, moments))

equilateral_row = next(row for row in control_rows if row[0] == "equilateral")
require(
    equilateral_row[2] == (sp.Rational(1, 2), 0, 0, sp.Rational(1, 20)),
    ("equilateral hostile moments", equilateral_row),
)
equilateral_p, equilateral_q = triangle_controls[2][1:]
equilateral_g = sp.expand(equilateral_p + I * equilateral_q)
equilateral_gbar = sp.expand(equilateral_p - I * equilateral_q)
finite_prefix = tuple(boundary_integral(equilateral_g**n, equilateral_gbar) for n in (2, 3, 4))
require(
    finite_prefix == (0, 0, -3 * sqrt3 * I / 5),
    ("finite current prefix hostile", finite_prefix),
)

source = Path(__file__).read_text(encoding="utf-8")
require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))), "assert node")

semantic = (
    face_checks,
    keller_flux_checks,
    triangle_checks,
    tuple((name, str(c), tuple(map(str, moments))) for name, c, moments in control_rows),
    tuple(map(str, finite_prefix)),
)
print("THM-3466 FACTORIAL FACE STOKES AND KELLER BOUNDARY CURRENT EXACT CONTROL")
print("orthant_face_checks=%d dimensions=3 codimensions=0..3 powers=1..5" % face_checks)
print("quadrant_keller_flux_checks=%d controls=identity,triangular" % keller_flux_checks)
print("triangle_boundary_checks=%d controls=identity,triangular,equilateral n=1..4" % triangle_checks)
for name, c, moments in control_rows:
    print("%s_J=%s moments_I0_I1_I2_I3=%s" % (name, c, moments))
print("equilateral_dbar_responses_n2_n3_n4=%s" % (finite_prefix,))
print("semantic_sha256=%s" % hashlib.sha256(repr(semantic).encode("utf-8")).hexdigest())
print("STATUS=PASS")
