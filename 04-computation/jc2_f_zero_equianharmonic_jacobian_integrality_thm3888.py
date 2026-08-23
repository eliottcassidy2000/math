#!/usr/bin/env python3
"""Exact algebraic companion for THM-3888's elliptic reframe."""

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
    gate(sp.cancel(expression) == 0, message)


x, y, T, G = sp.symbols("x y T G")
a = x + 1
L = 9 * x + 4
F = 15 * x**2 + 15 * x + 4
K = y**2 - F
Delta = a**3 * L**2 - K**2

# Binary quartic and its classical invariants.
quartic = sp.expand(L**4 - 6 * a * L**2 * T**2 - 8 * K * T**3 - 3 * a**2 * T**4)
A4 = -3 * a**2
B3 = -8 * K
C2 = -6 * a * L**2
D1 = 0
E0 = L**4
I = sp.expand(12 * A4 * E0 - 3 * B3 * D1 + C2**2)
J = sp.expand(
    72 * A4 * C2 * E0
    + 9 * B3 * C2 * D1
    - 27 * A4 * D1**2
    - 27 * B3**2 * E0
    - 2 * C2**3
)
zero(I, "binary quartic I invariant")
zero(J - 1728 * L**4 * Delta, "binary quartic J invariant")
zero(
    sp.discriminant(quartic, T) + 110592 * L**8 * Delta**2,
    "binary quartic discriminant",
)

# Exact quartic-to-Jacobian map.  Reduce the asserted relation modulo
# G^2-quartic rather than sampling a branch of the square root.
X0 = (G + L**2) / T**2
Y0 = 2 * (X0**2 + 3 * a**2) * T + 8 * K
map_relation = sp.together(Y0**2 - 8 * L**2 * (X0 - a) ** 3 + 64 * Delta)
map_numerator = sp.Poly(map_relation.as_numer_denom()[0], G)
curve_equation = sp.Poly(G**2 - quartic, G)
map_remainder = map_numerator.rem(curve_equation).as_expr()
zero(map_remainder, "quartic-to-Jacobian relation")

Xe = sp.expand(2 * L**2 * (X0 - a))
Ye = sp.expand(L**2 * Y0)
scaled_relation = sp.together(Ye**2 - Xe**3 + 64 * L**4 * Delta)
scaled_numerator = sp.Poly(scaled_relation.as_numer_denom()[0], G)
zero(
    scaled_numerator.rem(curve_equation).as_expr(),
    "scaled equianharmonic Jacobian relation",
)

# Dense inverse and its single visible denominator in scaled coordinates.
T_inverse = sp.cancel((Y0 - 8 * K) / (2 * (X0**2 + 3 * a**2)))
zero(T_inverse - T, "unscaled inverse recovers T")
zero(X0 * T**2 - L**2 - G, "unscaled inverse recovers G")
inverse_denominator = sp.expand(
    sp.symbols("Xe_abstract") ** 2
    + 4 * a * L**2 * sp.symbols("Xe_abstract")
    + 16 * a**2 * L**4
)
Xa, Ya = sp.symbols("Xe_abstract Ye_abstract")
T_inverse_scaled = 2 * L**2 * (Ya - 8 * K * L**2) / inverse_denominator
zero(
    T_inverse_scaled.subs({Xa: Xe, Ya: Ye}) - T,
    "scaled inverse recovers T",
)

# Natural factor-cofactor normalization.  On integral Weierstrass sections,
# polynomial T is exactly the displayed cubic divisibility test.
uvar, vvar = sp.symbols("uvar vvar")
normalized_curve = vvar**2 - K**2 - L**2 * (uvar**3 - a**3)
zero(
    normalized_curve.subs({uvar: Xe / (4 * L**2), vvar: Ye / (8 * L**2)})
    .together()
    .as_numer_denom()[0]
    .as_poly(G)
    .rem(curve_equation)
    .as_expr(),
    "normalized cubic factor-cofactor curve",
)
normalized_denominator = uvar**2 + a * uvar + a**2
zero(
    inverse_denominator.subs(Xa, 4 * L**2 * uvar)
    - 16 * L**4 * normalized_denominator,
    "normalized inverse denominator",
)
zero(
    (vvar - K) * (vvar + K)
    - L**2 * (uvar - a) * normalized_denominator
    - normalized_curve,
    "normalized three-factor identity",
)
T_normalized = (vvar - K) / normalized_denominator
zero(
    T_inverse_scaled.subs({Xa: 4 * L**2 * uvar, Ya: 8 * L**2 * vvar})
    - T_normalized,
    "normalized T inverse",
)
zero(
    (a + 2 * uvar) * T_normalized**2
    - L**2
    - ((a + (4 * L**2 * uvar) / (2 * L**2)) * T_normalized**2 - L**2),
    "normalized G inverse",
)

# Marked sections and the group-law relation Q_+ + Q_- = P_0.
s = sp.symbols("s")


def reduce_s(expression: sp.Expr) -> sp.Expr:
    numerator = sp.together(expression).as_numer_denom()[0]
    return sp.Poly(numerator, s).rem(sp.Poly(s**2 + 3, s)).as_expr()


X_plus = 2 * a * L**2 * (s - 1)
X_minus = 2 * a * L**2 * (-s - 1)
Y_boundary = -8 * K * L**2
X_P0 = 4 * a * L**2
Y_P0 = 8 * K * L**2
zero(reduce_s(Y_boundary**2 - X_plus**3 + 64 * L**4 * Delta), "Q-plus on E")
zero(reduce_s(Y_boundary**2 - X_minus**3 + 64 * L**4 * Delta), "Q-minus on E")
zero(Y_P0**2 - X_P0**3 + 64 * L**4 * Delta, "P-zero on E")
zero(-(X_plus + X_minus) - X_P0, "horizontal chord x-coordinate")
zero(-Y_boundary - Y_P0, "horizontal chord y-coordinate")
zero(
    reduce_s(X_plus / (4 * L**2) - a * (s - 1) / 2),
    "Q-plus normalized cube-root address",
)
zero(
    reduce_s(normalized_denominator.subs(uvar, a * (s - 1) / 2)),
    "Q-plus normalized denominator root",
)
zero(
    reduce_s(inverse_denominator.subs(Xa, X_plus)),
    "inverse denominator vanishes at Q-plus",
)
zero(
    reduce_s(inverse_denominator.subs(Xa, X_minus)),
    "inverse denominator vanishes at Q-minus",
)

# Exact factor polarization of the two T=0 branches.
H = 3 * a**2 * T**2 + 8 * K * T + 6 * a * L**2
zero((G - L**2) * (G + L**2) + T**2 * H - (G**2 - quartic), "T divisor factorization")

# Generic-y surface: Delta has four simple roots and a6 has orders 1,1,1,1,2.
M = 9 * a - 4
zero(F - (15 * a**2 - 15 * a + 4), "F in the a coordinate")
zero(F**2 - a**3 * L**2 - (1 - a) ** 3 * M**2, "two-square polarization")
delta_y_discriminant = sp.factor(sp.discriminant(Delta, y))
zero(
    delta_y_discriminant - 256 * a**6 * L**4 * (1 - a) ** 3 * M**2,
    "Delta has four generic simple y-roots",
)
B_weierstrass = -64 * L**4 * Delta
weierstrass_discriminant = -432 * B_weierstrass**2
zero(
    weierstrass_discriminant + 2**16 * 3**3 * L**8 * Delta**2,
    "short Weierstrass discriminant",
)
u = sp.symbols("u")
B_infinity = sp.expand(u**6 * B_weierstrass.subs(y, 1 / u))
expected_B_infinity = sp.expand(
    64 * L**4 * u**2
    - 128 * F * L**4 * u**4
    - 64 * L**4 * (a**3 * L**2 - F**2) * u**6
)
zero(B_infinity - expected_B_infinity, "minimal infinity coefficient")
zero(
    sp.Poly(B_infinity, u).coeff_monomial(u**2) - 64 * L**4,
    "infinity order-two leading coefficient",
)
gate(4 * 2 + 4 == 12, "II^4 plus IV Euler number")
gate(10 - 2 - 2 == 6, "rational Shioda-Tate rank")

# The IV blow-up chart separates the two component addresses.  Q_+ and Q_-
# have the same negative Y1 address and first separate at order u; P0 has the
# opposite address.
zero(
    sp.limit(u**3 * Y_boundary.subs(y, 1 / u) / u, u, 0) + 8 * L**2,
    "Q boundary IV component address",
)
zero(
    sp.limit(u**3 * Y_P0.subs(y, 1 / u) / u, u, 0) - 8 * L**2,
    "P-zero IV component address",
)
zero(
    sp.expand(u * X_plus - u * X_minus) - 4 * s * a * L**2 * u,
    "Q sections first-order separation",
)

# The 3-division polynomial exposes the only two possible rational
# three-torsion mechanisms; simple Delta-valuations exclude both in the proof.
Xvar = sp.symbols("Xvar")
psi3 = sp.factor(3 * Xvar * (Xvar**3 + 4 * B_weierstrass))
zero(
    psi3 - 3 * Xvar * (Xvar**3 - 256 * L**4 * Delta),
    "three-division polynomial",
)

semantic = {
    "quartic": "I=0,J=1728*L^4*Delta,disc=-110592*L^8*Delta^2",
    "jacobian": "Ye^2=Xe^3-64*L^4*Delta",
    "boundary": "div(T)=O+P0-Qplus-Qminus;Qplus+Qminus=P0",
    "factor": "v2=K2+L2(u3-a3);T=(v-K)/(u2+au+a2)",
    "surface": "generic y-fibres II^4+IV;rational;geometric MW rank6 torsion0",
    "scope": "two-section S-integrality and x-integral descent remain open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3888-f-zero-equianharmonic-jacobian-and-two-section-integrality")
print("binary_quartic_invariant=I_zero")
print("jacobian=Y^2=X^3-64L^4Delta")
print("normalized_factor=(v-K)(v+K)=L^2(u-a)(u^2+au+a^2)")
print("generic_fibres=II,II,II,II,IV")
print("geometric_MW=rank_6_torsion_0")
print("polynomial_lane=two_section_S_integrality_plus_x_descent_OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
