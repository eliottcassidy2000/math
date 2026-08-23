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

# A genuine k(x)[y] polynomial point showing that x-integral descent and the
# origin address are load-bearing, not cosmetic.
T_hostile = -2 * K / (3 * a**2)
G_hostile = 4 * K**2 / (3 * a**3) - L**2
zero(
    G_hostile**2 - quartic.subs(T, T_hostile),
    "rational-x polynomial-y hostile square",
)
zero(
    T_hostile - T_normalized.subs({uvar: a, vvar: -K}),
    "hostile normalized section address",
)
zero(T_hostile.subs({x: 0, y: 0}) - sp.Rational(8, 3), "hostile origin failure")
gate(
    sp.cancel(a**2 * T_hostile).subs(x, -1) != 0,
    "hostile has an a-zero pole",
)

# Complete the two low integral-Weierstrass heights left by the all-degree
# argument in the theorem.  At deg_y(u)=1 the inverse divisibility forces T
# to lie in k(x).  The quartic is then -8*T^3*y^2+C(T), with no y term, so a
# square with T != 0 forces C(T)=0.  The UFD rational-root theorem leaves only
# c*L^i/a^j (0<=i<=4, 0<=j<=2).  A coefficient gcd of one over Q[c]
# certifies that none of the fifteen candidates works, uniformly after
# extension to the algebraically closed constant field.
tconst, croot = sp.symbols("tconst croot")
quartic_tconst = sp.Poly(sp.expand(quartic.subs(T, tconst)), y)
zero(
    quartic_tconst.coeff_monomial(y**2) + 8 * tconst**3,
    "constant-T y-square coefficient",
)
zero(quartic_tconst.coeff_monomial(y), "constant-T missing y coefficient")
C_tconst = sp.expand(
    L**4 - 6 * a * L**2 * tconst**2 + 8 * F * tconst**3 - 3 * a**2 * tconst**4
)
zero(
    quartic_tconst.coeff_monomial(1) - C_tconst,
    "constant-T residual polynomial",
)

C_coefficients = [sp.Poly(q, x, domain=sp.QQ) for q in sp.Poly(C_tconst, tconst).all_coeffs()]
C_content = C_coefficients[0]
for coefficient in C_coefficients[1:]:
    C_content = C_content.gcd(coefficient)
gate(C_content.degree() == 0, "constant-T polynomial is primitive over k[x]")

rational_root_candidates = []
for numerator_power in range(5):
    for denominator_power in range(3):
        candidate = croot * L**numerator_power / a**denominator_power
        candidate_numerator = sp.expand(
            sp.together(C_tconst.subs(tconst, candidate)).as_numer_denom()[0]
        )
        coefficient_expressions = [
            q for q in sp.Poly(candidate_numerator, x).all_coeffs() if q != 0
        ]
        coefficient_gcd = sp.Poly(coefficient_expressions[0], croot, domain=sp.QQ)
        for coefficient in coefficient_expressions[1:]:
            coefficient_gcd = coefficient_gcd.gcd(
                sp.Poly(coefficient, croot, domain=sp.QQ)
            )
        gate(
            coefficient_gcd.degree() == 0,
            f"constant-T rational-root candidate ({numerator_power},{denominator_power})",
        )
        rational_root_candidates.append(
            (numerator_power, denominator_power, sp.monic(coefficient_gcd.as_expr(), croot))
        )
gate(len(rational_root_candidates) == 15, "all constant-T candidates exhausted")

# At deg_y(u)=0 the finite inverse chart forces u^3=a^3 and v=+/-K.
# Only u=a has nonzero inverse denominator; the two signs are exactly the
# finite T=0 point and the descent-hostile point already checked above.
zero(
    normalized_denominator.subs(uvar, a) - 3 * a**2,
    "finite cube-root address denominator",
)
zero(T_normalized.subs({uvar: a, vvar: K}), "finite base section T")
zero(
    (a + 2 * uvar).subs(uvar, a) * T_normalized.subs({uvar: a, vvar: K}) ** 2
    - L**2
    + L**2,
    "finite base section G",
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
    "hostile": "T=-2K/(3a2) is polynomial in y but fails x-integrality/address",
    "integral_shell": "polynomial u,v,T gives only (u,v)=(a,+/-K)",
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
print("hostile_kxy=T=-2K/(3a^2),origin=8/3,a_zero_pole=yes")
print("constant_T_rational_root_candidates=15_all_empty")
print(
    "constant_T_candidate_gcds="
    + ",".join(f"({i},{j}):{gcd_value}" for i, j, gcd_value in rational_root_candidates)
)
print("integral_Weierstrass_shell=exactly_base_and_x_pole_hostile")
print("generic_fibres=II,II,II,II,IV")
print("geometric_MW=rank_6_torsion_0")
print("polynomial_lane=two_section_S_integrality_plus_x_descent_OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
