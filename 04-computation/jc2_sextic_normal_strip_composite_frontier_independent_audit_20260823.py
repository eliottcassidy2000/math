#!/usr/bin/env python3
"""Independent exact audit of the composite sextic normal-strip scout.

This does not import the primary scout.  It reconstructs the depressed
packets by coefficient-space homotopy, proves the (6,4) Artin certificate by
explicit ideal/Buchberger identities, and checks the (6,5) apparent scheme in
a three-variable grevlex quotient by multiplication matrices.  The scope is
the two principal balanced faces and the primitive simple-pole jet only; it
is not a sextic normal-strip theorem and has no JC(2) consequence.
"""

from __future__ import annotations

import ast
import hashlib
import itertools
import json
import math
import sys
from pathlib import Path

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


GATES = 0


def gate(condition: object, label: str) -> None:
    """Optimization-safe truth gate."""

    global GATES
    GATES += 1
    if condition is True or condition == sp.S.true:
        return
    raise RuntimeError(f"GATE FAILED: {label}: {condition}")


def canon(expression: sp.Expr) -> sp.Expr:
    return sp.cancel(sp.factor(expression))


def zero(expression: sp.Expr, label: str) -> None:
    gate(canon(expression) == 0, label)


def total_derivative(
    expression: sp.Expr,
    coordinates: tuple[sp.Symbol, ...],
    differentials: tuple[sp.Symbol, ...],
) -> sp.Expr:
    return sp.expand(sum(
        sp.diff(expression, coordinate) * differential
        for coordinate, differential in zip(coordinates, differentials)
    ))


def bucket(
    first: list[sp.Expr],
    second: list[sp.Expr],
    degree: int,
    coordinates: tuple[sp.Symbol, ...],
    differentials: tuple[sp.Symbol, ...],
) -> sp.Expr:
    """Coefficient of x^degree in A_x C_s-A_s C_x."""

    answer = sp.S.Zero
    for i, first_i in enumerate(first):
        j = degree + 1 - i
        if 0 <= j < len(second):
            answer += i * first_i * total_derivative(
                second[j], coordinates, differentials
            )
            answer -= j * total_derivative(
                first_i, coordinates, differentials
            ) * second[j]
    return sp.expand(answer)


def radial_primitive(
    form: sp.Expr,
    coordinates: tuple[sp.Symbol, ...],
    differentials: tuple[sp.Symbol, ...],
    label: str,
) -> sp.Expr:
    """Recover a polynomial potential by the radial Poincare homotopy."""

    expanded = sp.expand(form)
    coefficients = tuple(expanded.coeff(differential) for differential in differentials)
    zero(
        expanded - sum(
            coefficient * differential
            for coefficient, differential in zip(coefficients, differentials)
        ),
        f"{label}: form is linear in the declared differentials",
    )
    for left in range(len(coordinates)):
        for right in range(left + 1, len(coordinates)):
            zero(
                sp.diff(coefficients[left], coordinates[right])
                - sp.diff(coefficients[right], coordinates[left]),
                f"{label}: closed coefficient pair {left},{right}",
            )
    tau = sp.Dummy("tau")
    dilation = {coordinate: tau * coordinate for coordinate in coordinates}
    integrand = sum(
        coordinate * coefficient.subs(dilation, simultaneous=True)
        for coordinate, coefficient in zip(coordinates, coefficients)
    )
    potential = sp.cancel(sp.integrate(sp.expand(integrand), (tau, 0, 1)))
    zero(
        total_derivative(potential, coordinates, differentials) - expanded,
        f"{label}: radial primitive differentiates back",
    )
    return potential


def derive_packet(
    target_degree: int,
    target: list[sp.Expr],
    coordinates: tuple[sp.Symbol, ...],
    differentials: tuple[sp.Symbol, ...],
    integration_constants: tuple[sp.Symbol, ...],
    label: str,
) -> list[sp.Expr]:
    """Solve every triangular high row without using a displayed formula."""

    first = [sp.S.Zero for _ in range(7)]
    first[6] = R
    first[5] = D
    zero(
        bucket(
            first,
            target,
            target_degree + 4,
            coordinates,
            differentials,
        ),
        f"{label}: D is a constant",
    )
    for coefficient_degree, constant in zip(
        range(4, -1, -1), integration_constants
    ):
        row_degree = coefficient_degree + target_degree - 1
        known_form = bucket(
            first, target, row_degree, coordinates, differentials
        )
        derivative = sp.cancel(known_form / (target_degree * Q))
        first[coefficient_degree] = radial_primitive(
            derivative,
            coordinates,
            differentials,
            f"{label}: solve A_{coefficient_degree}",
        ) + constant
        zero(
            bucket(first, target, row_degree, coordinates, differentials),
            f"{label}: reconstructed high bucket x^{row_degree}",
        )
    return first


def proportional(
    left: sp.Expr,
    right: sp.Expr,
    generators: tuple[sp.Symbol, ...],
    label: str,
) -> sp.Expr:
    """Return and certify the nonzero scalar relating two polynomials."""

    left_poly = sp.Poly(sp.expand(left), *generators, domain="QQ")
    right_poly = sp.Poly(sp.expand(right), *generators, domain="QQ")
    gate(not left_poly.is_zero and not right_poly.is_zero, f"{label}: nonzero inputs")
    left_terms = left_poly.as_dict()
    right_terms = right_poly.as_dict()
    exponent = next(iter(right_terms))
    ratio = sp.cancel(left_terms.get(exponent, 0) / right_terms[exponent])
    gate(ratio != 0, f"{label}: nonzero ratio")
    zero(left - ratio * right, f"{label}: exact proportionality")
    return ratio


def laurent_leader(expression: sp.Expr, uniformizer: sp.Symbol) -> tuple[int, sp.Expr]:
    """Return the least Laurent exponent and its coefficient."""

    terms = sp.expand(expression).as_ordered_terms()
    exponents = [int(term.as_powers_dict().get(uniformizer, 0)) for term in terms]
    least = min(exponents)
    coefficient = sp.expand(sum(
        term / uniformizer**least
        for term, exponent in zip(terms, exponents)
        if exponent == least
    ))
    gate(not coefficient.has(uniformizer), "Laurent leader is uniformizer-free")
    return least, coefficient


# ---------------------------------------------------------------------------
# 1. Universal top row, UFD exponent sieve, and exact shear rows.
# ---------------------------------------------------------------------------
x = sp.symbols("x")
R, Q = sp.symbols("R Q", nonzero=True)
D, E, F, G, I, H, A0 = sp.symbols("D E F G I H A0")

gcd_packets = {}
for target_degree in range(1, 6):
    divisor = math.gcd(6, target_degree)
    exponent_w = 6 // divisor
    exponent_q = target_degree // divisor
    gcd_packets[target_degree] = (exponent_w, exponent_q)
    gate(
        6 * exponent_q == target_degree * exponent_w,
        f"top valuation equation j={target_degree}",
    )
    gate(
        math.gcd(exponent_w, exponent_q) == 1,
        f"top valuation packet is primitive j={target_degree}",
    )

gate(
    gcd_packets == {1: (6, 1), 2: (3, 1), 3: (2, 1), 4: (3, 2), 5: (6, 5)},
    "complete sextic gcd packet table",
)

h, z = sp.symbols("h z", nonzero=True)
for target_degree in (1, 2, 3):
    exponent_w, exponent_q = gcd_packets[target_degree]
    top_A = R * h**exponent_w * z**6
    top_C = Q * h**exponent_q * z**target_degree
    shear = sp.cancel(R / Q ** (6 // target_degree))
    zero(
        top_A - shear * top_C ** (6 // target_degree),
        f"target shear kills sextic top at j={target_degree}",
    )

gate(
    all(6 % degree != 0 for degree in (4, 5)),
    "only residual rows j=4,5 resist a monomial target shear",
)

# A constant target has bracket A_x*b'.  Its top x^5 coefficient forces b'=0,
# after which every coefficient vanishes.  This is the exact j=0 exclusion.
w, db0 = sp.symbols("w db0", nonzero=True)
gate(6 * w * db0 != 0, "j=0 top coefficient is nonzero when b is nonconstant")


# ---------------------------------------------------------------------------
# 2. Independently reconstruct the complete depressed (6,4) packet.
# ---------------------------------------------------------------------------
B4, N4, b4 = sp.symbols("B4 N4 b4")
dB4, dN4, db4 = sp.symbols("dB4 dN4 db4")
coords64 = (B4, N4, b4)
diffs64 = (dB4, dN4, db4)
C64 = [b4, N4, B4, sp.S.Zero, Q]
A64 = derive_packet(4, C64, coords64, diffs64, (E, F, G, I, A0), "(6,4)")

K4 = 3 * R / (2 * Q)
H4 = 5 * D / (4 * Q)
expected64 = [sp.S.Zero for _ in range(7)]
expected64[6] = R
expected64[5] = D
expected64[4] = K4 * B4 + E
expected64[3] = K4 * N4 + H4 * B4 + F
expected64[2] = (
    K4 * b4 + H4 * N4 + K4 * B4**2 / (4 * Q) + E * B4 / Q + G
)
expected64[1] = (
    5 * D * b4 / (4 * Q)
    + K4 * B4 * N4 / (2 * Q)
    + E * N4 / Q
    + H4 * B4**2 / (8 * Q)
    + 3 * F * B4 / (4 * Q)
    + I
)
expected64[0] = (
    -R * B4**3
    + 5 * D * Q * B4 * N4
    + 8 * G * Q**2 * B4
    + 12 * Q * R * B4 * b4
    + 16 * E * Q**2 * b4
    + 12 * F * Q**2 * N4
    + 6 * Q * R * N4**2
) / (16 * Q**3) + A0
for degree in range(7):
    zero(A64[degree] - expected64[degree], f"(6,4) derived coefficient A_{degree}")

row64_2 = bucket(A64, C64, 2, coords64, diffs64)
row64_1 = bucket(A64, C64, 1, coords64, diffs64)
potential64_2 = radial_primitive(row64_2, coords64, diffs64, "(6,4) x2 conserved row")
potential64_1 = radial_primitive(row64_1, coords64, diffs64, "(6,4) x1 conserved row")
expected64_2 = (
    -5 * B4**3 * D
    - 12 * B4**2 * F * Q
    - 24 * B4**2 * N4 * R
    + 40 * B4 * D * Q * b4
    + 32 * B4 * I * Q**2
    + 20 * D * N4**2 * Q
    + 96 * F * Q**2 * b4
    + 64 * G * N4 * Q**2
    + 96 * N4 * Q * R * b4
) / (32 * Q**2)
expected64_1 = (
    3 * B4**4 * R
    - 15 * B4**2 * D * N4 * Q
    - 16 * B4**2 * G * Q**2
    - 24 * B4**2 * Q * R * b4
    - 24 * B4 * F * N4 * Q**2
    - 24 * B4 * N4**2 * Q * R
    + 40 * D * N4 * Q**2 * b4
    + 64 * G * Q**3 * b4
    + 32 * I * N4 * Q**3
    + 48 * Q**2 * R * b4**2
) / (32 * Q**3)
zero(potential64_2 - expected64_2, "(6,4) independently recovered first conserved polynomial")
zero(potential64_1 - expected64_1, "(6,4) independently recovered second conserved polynomial")
zero(
    bucket(A64, C64, 0, coords64, diffs64)
    - (A64[1] * db4 - total_derivative(A64[0], coords64, diffs64) * N4),
    "(6,4) constant one-form shape",
)


# ---------------------------------------------------------------------------
# 3. Independently reconstruct the complete depressed (6,5) packet.
# ---------------------------------------------------------------------------
B5, N5, V5, b5 = sp.symbols("B5 N5 V5 b5")
dB5, dN5, dV5, db5 = sp.symbols("dB5 dN5 dV5 db5")
coords65 = (B5, N5, V5, b5)
diffs65 = (dB5, dN5, dV5, db5)
C65 = [b5, V5, N5, B5, sp.S.Zero, Q]
A65 = derive_packet(5, C65, coords65, diffs65, (E, F, G, H, A0), "(6,5)")

K5 = 6 * R / (5 * Q)
expected65 = [sp.S.Zero for _ in range(7)]
expected65[6] = R
expected65[5] = D
expected65[4] = K5 * B5 + E
expected65[3] = K5 * N5 + D * B5 / Q + F
expected65[2] = (
    K5 * V5 + D * N5 / Q + K5 * B5**2 / (10 * Q)
    + 4 * E * B5 / (5 * Q) + G
)
expected65[1] = (
    K5 * b5 + D * V5 / Q + K5 * B5 * N5 / (5 * Q)
    + 3 * F * B5 / (5 * Q) + 4 * E * N5 / (5 * Q) + H
)
expected65[0] = (
    -4 * B5**3 * R
    - 10 * B5**2 * E * Q
    + 50 * B5 * G * Q**2
    + 30 * B5 * Q * R * V5
    + 125 * D * Q**2 * b5
    + 100 * E * Q**2 * V5
    + 75 * F * N5 * Q**2
    + 15 * N5**2 * Q * R
) / (125 * Q**3) + A0
for degree in range(7):
    zero(A65[degree] - expected65[degree], f"(6,5) derived coefficient A_{degree}")

row65_3 = bucket(A65, C65, 3, coords65, diffs65)
row65_2 = bucket(A65, C65, 2, coords65, diffs65)
potential65_3 = radial_primitive(row65_3, coords65, diffs65, "(6,5) x3 conserved row")
potential65_2 = radial_primitive(row65_2, coords65, diffs65, "(6,5) x2 conserved row")
expected65_3 = (
    -15 * B5**2 * F * Q
    - 12 * B5**2 * N5 * R
    - 20 * B5 * E * N5 * Q
    + 25 * B5 * H * Q**2
    + 30 * B5 * Q * R * b5
    + 100 * E * Q**2 * b5
    + 75 * F * Q**2 * V5
    + 50 * G * N5 * Q**2
    + 30 * N5 * Q * R * V5
) / (25 * Q**2)
expected65_2 = (
    9 * B5**4 * R
    + 20 * B5**3 * E * Q
    - 75 * B5**2 * G * Q**2
    - 60 * B5**2 * Q * R * V5
    - 100 * B5 * E * Q**2 * V5
    - 150 * B5 * F * N5 * Q**2
    - 60 * B5 * N5**2 * Q * R
    - 50 * E * N5**2 * Q**2
    + 375 * F * Q**3 * b5
    + 250 * G * Q**3 * V5
    + 125 * H * N5 * Q**3
    + 150 * N5 * Q**2 * R * b5
    + 75 * Q**2 * R * V5**2
) / (125 * Q**3)
zero(potential65_3 - expected65_3, "(6,5) independently recovered first conserved polynomial")
zero(potential65_2 - expected65_2, "(6,5) independently recovered second conserved polynomial")

row65_1 = bucket(A65, C65, 1, coords65, diffs65)
coeff65_B = sp.expand(row65_1).coeff(dB5)
coeff65_N = sp.expand(row65_1).coeff(dN5)
curl65_BN = canon(sp.diff(coeff65_B, N5) - sp.diff(coeff65_N, B5))
gate(curl65_BN != 0, "(6,5) x1 row is genuinely nonclosed")
zero(
    row65_1 - (
        2 * A65[2] * db5
        - 2 * N5 * total_derivative(A65[0], coords65, diffs65)
        + A65[1] * dV5
        - V5 * total_derivative(A65[1], coords65, diffs65)
    ),
    "(6,5) x1 one-form shape",
)
zero(
    bucket(A65, C65, 0, coords65, diffs65)
    - (A65[1] * db5 - V5 * total_derivative(A65[0], coords65, diffs65)),
    "(6,5) constant one-form shape",
)


# ---------------------------------------------------------------------------
# 4. Recover both principal balanced face packets from the derived formulas.
# ---------------------------------------------------------------------------
t = sp.symbols("t", nonzero=True)
X4, Y4, Z4 = sp.symbols("X4 Y4 Z4")
gshift = 1 / t
sub64 = {
    B4: X4 / t**2,
    N4: Y4 / t**3,
    b4: Z4 / t**4,
    R: 1,
    Q: 1,
    D: 0,
    E: 0,
    F: 0,
    G: 0,
    I: 0,
    A0: 0,
}
armC64 = sp.expand(sum(C64[index] * gshift**index for index in range(len(C64))).subs(sub64) * t**4)
armA64 = sp.expand(sum(A64[index] * gshift**index for index in range(7)).subs(sub64) * t**6)
face64_2 = sp.expand(potential64_2.subs(sub64) * t**7)
face64_1 = sp.expand(potential64_1.subs(sub64) * t**8)
expected_armC64 = 1 + X4 + Y4 + Z4
expected_face64_2 = Y4 * (4 * Z4 - X4**2)
expected_face64_1 = (X4**2 - 4 * Z4) ** 2 - 8 * X4 * Y4**2
expected_armA64 = 6 * Y4**2 - (X4 + 2) ** 3
zero(armC64 - expected_armC64, "(6,4) target arm recovered from packet")
armA64_reduced = sp.expand(armA64.subs(Z4, -1 - X4 - Y4))
ratio64_A = proportional(armA64_reduced, expected_armA64, (X4, Y4), "(6,4) source arm modulo target arm")
ratio64_2 = proportional(face64_2, expected_face64_2, (X4, Y4, Z4), "(6,4) first face conservation")
ratio64_1 = proportional(face64_1, expected_face64_1, (X4, Y4, Z4), "(6,4) second face conservation")

row64_0 = bucket(A64, C64, 0, coords64, diffs64)
sub64_d = dict(sub64)
sub64_d.update({
    dB4: sp.diff(X4 / t**2, t),
    dN4: sp.diff(Y4 / t**3, t),
    db4: sp.diff(Z4 / t**4, t),
})
exponent64_0, leader64_0 = laurent_leader(row64_0.subs(sub64_d), t)
expected_bracket64 = Y4 * (X4**3 - 4 * X4 * Z4 - 6 * Y4**2)
ratio64_0 = proportional(leader64_0, expected_bracket64, (X4, Y4, Z4), "(6,4) constant face bracket")

X5, Y5, Z5, W5 = sp.symbols("X5 Y5 Z5 W5")
sub65 = {
    B5: X5 / t**2,
    N5: Y5 / t**3,
    V5: Z5 / t**4,
    b5: W5 / t**5,
    R: 1,
    Q: 1,
    D: 0,
    E: 0,
    F: 0,
    G: 0,
    H: 0,
    A0: 0,
}
armC65 = sp.expand(sum(C65[index] * gshift**index for index in range(len(C65))).subs(sub65) * t**5)
armA65 = sp.expand(sum(A65[index] * gshift**index for index in range(7)).subs(sub65) * t**6)
face65_3 = sp.expand(potential65_3.subs(sub65) * t**7)
face65_2 = sp.expand(potential65_2.subs(sub65) * t**8)
expected_armC65 = 1 + X5 + Y5 + Z5 + W5
expected_face65_3 = -2 * X5**2 * Y5 + 5 * X5 * W5 + 5 * Y5 * Z5
expected_face65_2 = 3 * X5**4 - 20 * X5**2 * Z5 - 20 * X5 * Y5**2 + 50 * Y5 * W5 + 25 * Z5**2
expected_armA65 = -4 * X5**3 + 15 * X5**2 + 30 * X5 * Y5 + 30 * X5 * Z5 + 15 * Y5**2 - 25
zero(armC65 - expected_armC65, "(6,5) target arm recovered from packet")
armA65_reduced = sp.expand(armA65.subs(W5, -1 - X5 - Y5 - Z5))
ratio65_A = proportional(armA65_reduced, expected_armA65, (X5, Y5, Z5), "(6,5) source arm modulo target arm")
ratio65_3 = proportional(face65_3, expected_face65_3, (X5, Y5, Z5, W5), "(6,5) first face conservation")
ratio65_2 = proportional(face65_2, expected_face65_2, (X5, Y5, Z5, W5), "(6,5) second face conservation")

sub65_d = dict(sub65)
sub65_d.update({
    dB5: sp.diff(X5 / t**2, t),
    dN5: sp.diff(Y5 / t**3, t),
    dV5: sp.diff(Z5 / t**4, t),
    db5: sp.diff(W5 / t**5, t),
})
exponent65_1, leader65_1 = laurent_leader(row65_1.subs(sub65_d), t)
row65_0 = bucket(A65, C65, 0, coords65, diffs65)
exponent65_0, leader65_0 = laurent_leader(row65_0.subs(sub65_d), t)
expected_row65_1 = 25 * W5 * X5**2 + 225 * W5 * Z5 + 8 * X5**3 * Y5 - 65 * X5 * Y5 * Z5 - 30 * Y5**3
expected_row65_0 = 125 * W5**2 + 25 * W5 * X5 * Y5 + 4 * X5**3 * Z5 - 30 * X5 * Z5**2 - 15 * Y5**2 * Z5
ratio65_1 = proportional(leader65_1, expected_row65_1, (X5, Y5, Z5, W5), "(6,5) x1 face row")
ratio65_0 = proportional(leader65_0, expected_row65_0, (X5, Y5, Z5, W5), "(6,5) constant face row")


# ---------------------------------------------------------------------------
# 5. (6,4): explicit ideal equality, length, support, and bracket certificate.
# ---------------------------------------------------------------------------
u, v = sp.symbols("u v")
a64 = 6 * v**2 - u**3
p64 = -v * (u**2 + 4 * v)
q64 = (u**2 + 4 * v) ** 2 - 8 * (u - 2) * v**2
g64_1 = a64
g64_2 = u**2 * (2 * u + 3 * v)
g64_3 = u**4

zero(g64_2 - (-2 * a64 - 3 * p64), "(6,4) ideal inclusion g2 in face ideal")
zero(g64_3 - (-3 * q64 - 4 * u * a64 - 24 * p64), "(6,4) ideal inclusion g3 in face ideal")
zero(p64 + (2 * g64_1 + g64_2) / 3, "(6,4) reverse ideal inclusion p")
zero(
    q64 - (-g64_3 - 4 * u * g64_1 + 16 * g64_1 + 8 * g64_2) / 3,
    "(6,4) reverse ideal inclusion q",
)

# Buchberger's three S-pairs reduce by displayed identities for lex v>u.
s12 = sp.expand(u**2 * g64_1 - 2 * v * g64_2)
s13 = sp.expand(u**4 * g64_1 - 6 * v**2 * g64_3)
s23 = sp.expand(u**2 * g64_2 - 3 * v * g64_3)
zero(s12 - ((-u + sp.Rational(8, 3)) * g64_3 - sp.Rational(4, 3) * u * g64_2), "(6,4) Buchberger S12")
zero(s13 + u**3 * g64_3, "(6,4) Buchberger S13")
zero(s23 - 2 * u * g64_3, "(6,4) Buchberger S23")

leading64 = ((2, 0), (1, 2), (0, 4))  # exponents in the order (v,u)
standard64_exponents = tuple(
    exponent
    for exponent in itertools.product(range(2), range(4))
    if not any(
        exponent[0] >= leader[0] and exponent[1] >= leader[1]
        for leader in leading64
    )
)
gate(
    standard64_exponents == ((0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1)),
    "(6,4) exact six standard monomials",
)
support64 = sp.solve_poly_system([g64_1, g64_3], u, v)
gate(support64 == [(0, 0)], "(6,4) unique geometric support")
gate(len(standard64_exponents) > len(support64), "(6,4) local algebra is nonreduced")

bracket64_local = sp.expand(
    v * ((u - 2) ** 3 - 4 * (u - 2) * (1 - u - v) - 6 * v**2)
)
certificate64 = (
    (-v + 2 * u / 3 - sp.Rational(4, 3)) * g64_1
    - sp.Rational(2, 3) * g64_2
    + sp.Rational(2, 3) * g64_3
)
zero(bracket64_local - certificate64, "(6,4) explicit leading-bracket certificate")


# Primitive simple-pole jet reconstructed from A(g) and the first conserved
# polynomial, rather than from the primary scout's displayed face formulas.
u1, u2, u3, v1, v2, v3 = sp.symbols("u1 u2 u3 v1 v2 v3")
c4, c5, c6 = sp.symbols("c4 c5 c6")
Xjet = -2 + u1 * t + u2 * t**2 + u3 * t**3
Yjet = v1 * t + v2 * t**2 + v3 * t**3
Zjet = 1 - (Xjet + 2) - Yjet + c4 * t**4 + c5 * t**5 + c6 * t**6
jet_sub = {X4: Xjet, Y4: Yjet, Z4: Zjet}

# Reinsert the lower-weight constants in the already reconstructed faces.
full64 = {
    B4: Xjet / t**2,
    N4: Yjet / t**3,
    b4: Zjet / t**4,
    R: 1,
    Q: 1,
}
Ajet64 = sp.expand(sum(A64[index] * gshift**index for index in range(7)).subs(full64) * t**6)
Kjet64 = sp.expand(potential64_2.subs(full64) * t**7)
zero(Ajet64.coeff(t, 1) - 3 * D / 8, "(6,4) primitive jet t1 source arm")
zero(Kjet64.coeff(t, 1) + 5 * D / 4, "(6,4) primitive jet t1 conserved row")
zero(Ajet64.coeff(t, 2).subs(D, 0) - 3 * v1**2 / 8, "(6,4) primitive jet t2")
zero(Kjet64.coeff(t, 3).subs({D: 0, v1: 0}) - 3 * F / 2, "(6,4) primitive jet t3 conserved")
zero(Ajet64.coeff(t, 3).subs({D: 0, v1: 0, F: 0}) + u1**3 / 16, "(6,4) primitive jet t3 arm")
zero(Kjet64.coeff(t, 4).subs({D: 0, v1: 0, F: 0, u1: 0}) + 3 * v2**2, "(6,4) primitive jet t4")
zero(Ajet64.coeff(t, 5).subs({D: 0, v1: 0, F: 0, u1: 0, v2: 0}) - I, "(6,4) primitive jet t5")
zero(Kjet64.coeff(t, 6).subs({D: 0, v1: 0, F: 0, u1: 0, v2: 0, I: 0}) + 3 * v3**2, "(6,4) primitive jet t6")


# ---------------------------------------------------------------------------
# 6. (6,5): W-elimination, grevlex quotient, and multiplication-unit audit.
# ---------------------------------------------------------------------------
W_eliminated = -1 - X5 - Y5 - Z5
f3 = sp.expand(expected_face65_3.subs(W5, W_eliminated))
f2 = sp.expand(expected_face65_2.subs(W5, W_eliminated))
a_face = sp.expand(expected_armA65.subs(W5, W_eliminated))
f1 = sp.expand(expected_row65_1.subs(W5, W_eliminated))
f0 = sp.expand(expected_row65_0.subs(W5, W_eliminated))

quotient65 = sp.groebner([f3, f2, a_face], X5, Y5, Z5, order="grevlex")
leading65 = {poly.LM(order=quotient65.order).exponents for poly in quotient65.polys}
expected_leading65 = {
    (0, 2, 2), (1, 0, 3), (0, 1, 3), (0, 0, 4), (3, 0, 0),
    (2, 1, 0), (1, 2, 0), (0, 3, 0), (2, 0, 1), (1, 1, 1),
}
gate(leading65 == expected_leading65, "(6,5) independent grevlex initial ideal")

basis65_exponents = tuple(
    exponent
    for exponent in itertools.product(range(3), range(3), range(4))
    if not any(
        all(exponent[index] >= leader[index] for index in range(3))
        for leader in leading65
    )
)
gate(len(basis65_exponents) == 14, "(6,5) quotient algebra has dimension fourteen")
basis65 = tuple(
    X5**exponent[0] * Y5**exponent[1] * Z5**exponent[2]
    for exponent in basis65_exponents
)


def quotient_vector65(expression: sp.Expr) -> sp.Matrix:
    remainder = sp.expand(quotient65.reduce(sp.expand(expression))[1])
    dictionary = sp.Poly(remainder, X5, Y5, Z5, domain="QQ").as_dict()
    gate(
        set(dictionary).issubset(set(basis65_exponents)),
        "(6,5) quotient remainder lies in standard basis",
    )
    return sp.Matrix([dictionary.get(exponent, 0) for exponent in basis65_exponents])


def multiplication65(expression: sp.Expr) -> sp.Matrix:
    return sp.Matrix.hstack(*(
        quotient_vector65(expression * monomial) for monomial in basis65
    ))


MX65 = multiplication65(X5)
chi = sp.symbols("chi")
characteristic65 = sp.Poly(MX65.charpoly(chi).as_expr(), chi, domain="QQ")
denominator65, integral65 = characteristic65.clear_denoms(convert=True)
content65, primitive65 = integral65.primitive()
if primitive65.LC() < 0:
    primitive65 = -primitive65
coefficients65 = [int(coefficient) for coefficient in primitive65.all_coeffs()]
terminal65_sha256 = hashlib.sha256(
    ",".join(str(coefficient) for coefficient in coefficients65).encode("ascii")
).hexdigest()
gate(primitive65.degree() == 14, "(6,5) multiplication-X characteristic degree fourteen")
gate(sp.gcd(primitive65, primitive65.diff()).degree() == 0, "(6,5) multiplication-X characteristic is squarefree")
gate(
    terminal65_sha256 == "0e3c6807f6d25f289fd3e1116a59d4aa01749e838efd6e6a901a051df7eab4ef",
    "(6,5) independent terminal coefficient hash",
)

MF1 = multiplication65(f1)
MF0 = multiplication65(f0)
detF1 = sp.factor(MF1.det(method="domain-ge"))
detF0 = sp.factor(MF0.det(method="domain-ge"))
gate(detF1 != 0, "(6,5) F1 multiplication is invertible: unit ideal")
gate(detF0 != 0, "(6,5) F0 multiplication is invertible: nonzero at all points")
unit_determinant_sha256 = hashlib.sha256(
    (str(detF1) + "|" + str(detF0)).encode("ascii")
).hexdigest()


# ---------------------------------------------------------------------------
# 7. Exact missing-channel ledger and frozen semantic packet.
# ---------------------------------------------------------------------------
dominant64 = tuple(
    subset
    for size in range(2, 4)
    for subset in itertools.combinations(("B", "N", "b"), size)
)
dominant65 = tuple(
    subset
    for size in range(2, 5)
    for subset in itertools.combinations(("B", "N", "V", "b"), size)
)
gate(len(dominant64) == 4, "four nonprincipal dominant-tie channels in (6,4)")
gate(len(dominant65) == 11, "eleven nonprincipal dominant-tie channels in (6,5)")

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "no optimization-disabled Python assert",
)

semantic = {
    "status": "INDEPENDENT_FINITE_EXACT_AUDIT;SCOUT_SCOPE_ONLY;JC2_OPEN",
    "top": "UFD gcd packets and j=1,2,3 target shears rederived",
    "packets": "all high coefficients recovered by radial coefficient-space homotopy",
    "conserved64": "x2 and x1 closed forms independently integrated",
    "conserved65": "x3 and x2 closed forms independently integrated;x1 nonclosed",
    "face64": "explicit ideal equality+Buchberger length6 support one bracket certificate",
    "jet64": "primitive simple-pole forces D,F,I,u1,v1,v2,v3 zero only",
    "face65": "W eliminated;3-variable grevlex dimension14;squarefree MX;F1,F0 multiplication units",
    "missing64": dominant64,
    "missing65": dominant65,
    "scope": "nonprincipal fans,regular shift,higher pole,unit/infinity,zero carriers,descent remain open",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("JC2_SEXTIC_COMPOSITE_FRONTIER_INDEPENDENT_AUDIT")
print("status=FINITE_EXACT_AUDIT;PRIMARY_SCOUT_CONFIRMED;NOT_SEXTIC_THEOREM;JC2_OPEN")
print(f"top_gcd_packets={gcd_packets};shears=j1,j2,j3")
print("packets=RADIAL_HOMOTOPY_REDERIVATION_PASS;(6,4)_two_conserved;(6,5)_two_conserved+x1_nonclosed")
print(f"face64_ratios=A:{ratio64_A},K2:{ratio64_2},K1:{ratio64_1},F0:{ratio64_0};laurent_F0={exponent64_0}")
print("face64=SUPPORT(-2,0,1);ARTIN_LENGTH6;NONREDUCED;BRACKET_CERTIFICATE_PASS")
print("jet64=PRIMITIVE_SIMPLE_POLE_ONLY;D=F=I=0;u=O(t^2);Y=O(t^4)")
print(f"face65_ratios=A:{ratio65_A},F3:{ratio65_3},F2:{ratio65_2},F1:{ratio65_1},F0:{ratio65_0};laurent=({exponent65_1},{exponent65_0})")
print(f"face65=GREVLEX_QUOTIENT_DIM14;MX_SQUAREFREE;F1_UNIT;F0_UNIT;terminal_sha256={terminal65_sha256}")
print(f"face65_unit_determinant_sha256={unit_determinant_sha256}")
print(f"missing_nonprincipal_dominant_ties=(6,4):{len(dominant64)},(6,5):{len(dominant65)}")
print("scope_open=regular_shift+higher_pole+nonprincipal_fans+unit_infinity+zero_carriers+descent_integrality")
print(f"semantic_sha256={semantic_sha256}")
print(f"GATES={GATES}")
print("RESULT=PASS")
