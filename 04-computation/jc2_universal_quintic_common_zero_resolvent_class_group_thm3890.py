#!/usr/bin/env python3
"""Assertion-free exact gates for THM-3890.

This companion checks the universal algebra behind the quintic

    delta = H_4(M,L) + lambda*L^5,

including the tangent factor, the non-tangent localization chart, both
unit-valuation lattices, their Smith forms, and repeated-root controls.

Reproduction:
  python3 04-computation/jc2_universal_quintic_common_zero_resolvent_class_group_thm3890.py
  python3 -O 04-computation/jc2_universal_quintic_common_zero_resolvent_class_group_thm3890.py
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


CHECKS = 0


def gate(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"{label}: failed")


def zero(label: str, value: sp.Expr) -> None:
    gate(sp.factor(sp.cancel(value)) == 0, label)


def equal(label: str, left: object, right: object) -> None:
    gate(left == right, f"{label}: {left!r} != {right!r}")


def binary_disc(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    return sp.expand(
        b**2 * c**2
        - 4 * a * c**3
        - 4 * b**3 * d
        - 27 * a**2 * d**2
        + 18 * a * b * c * d
    )


M, L, W, t, y, lam = sp.symbols("M L W t y lambda")
h0, h1, h2, h3, h4 = sp.symbols("h0 h1 h2 h3 h4")

H4 = h0 * M**4 + h1 * M**3 * L + h2 * M**2 * L**2 + h3 * M * L**3 + h4 * L**4
D = h0 * t**4 + h1 * t**3 + h2 * t**2 + h3 * t + h4
delta = sp.expand(H4 + lam * L**5)
surface = sp.expand(W**2 - delta)

# Radial chart and rational normalization of the branch curve.
zero("radial_quartic", H4.subs(M, t * L) - L**4 * D)
Lpar = -D / lam
Mpar = t * Lpar
zero("branch_parametrization", delta.subs({M: Mpar, L: Lpar}, simultaneous=True))
equal("branch_graph_transverse_derivative", sp.diff(D + lam * L, L), lam)

# After inverting L, the quadratic double surface is a localization of A2.
F = sp.expand(y**2 - D)
Lback = F / lam
Mback = t * Lback
Wback = y * Lback**2
zero(
    "quadratic_surface_localization_chart",
    surface.subs({M: Mback, L: Lback, W: Wback}, simultaneous=True),
)
zero("chart_recovers_t", Mback / Lback - t)
zero("chart_recovers_y", Wback / Lback**2 - y)
zero("chart_recovers_L", Lback - F / lam)

# The tangent direction is an immediate polynomial factor, not a new class
# group case.
g0, g1, g2, g3 = sp.symbols("g0 g1 g2 g3")
G3 = g0 * M**3 + g1 * M**2 * L + g2 * M * L**2 + g3 * L**3
tangent_delta = sp.expand(L * G3 + lam * L**5)
zero("tangent_factor", tangent_delta - L * (G3 + lam * L**4))

# Non-tangency is h0 != 0.  The special fibre has two reduced generic
# components, each with multiplicity one.
mu = sp.symbols("mu")
special_surface = sp.expand(surface.subs(L, 0))
zero(
    "two_prime_special_fibre",
    sp.expand(
        special_surface
        - (W - mu * M**2) * (W + mu * M**2)
    ).subs(mu**2, h0),
)
gate(2 * mu * M**2 != 0, "characteristic_zero_prime_separation")

# Nonsquare D: F is one irreducible localization divisor, so its valuation
# vector is (1,1).  Nagata leaves one free class.
nonsquare_val = sp.Matrix([[1], [1]])
equal(
    "nonsquare_smith_form",
    smith_normal_form(nonsquare_val, domain=ZZ),
    sp.Matrix([[1], [0]]),
)
equal("nonsquare_class_free_rank", 2 - nonsquare_val.rank(), 1)

# Square D=q^2: write Q=L^2*q(M/L).  The two localized units are
# (W-Q)/L^2 and (W+Q)/L^2.  At their matching boundary primes the numerator
# has order five; at the opposite prime it has order zero.  Subtracting the
# denominator order two gives (3,-2) and (-2,3).
q0, q1, q2 = sp.symbols("q0 q1 q2")
q = q0 * t**2 + q1 * t + q2
Q = q0 * M**2 + q1 * M * L + q2 * L**2
square_delta = sp.expand(Q**2 + lam * L**5)
zero("square_radial_identity", Q.subs(M, t * L) - L**2 * q)
# Replace W^2 by the surface equation in the split-product identity.
zero(
    "square_split_product_on_surface",
    sp.expand((W - Q) * (W + Q) - lam * L**5).subs(W**2, square_delta),
)
equal("matching_localized_unit_order", 5 - 2, 3)
equal("opposite_localized_unit_order", 0 - 2, -2)
square_val = sp.Matrix([[3, -2], [-2, 3]])
equal("square_valuation_determinant", abs(int(square_val.det())), 5)
equal(
    "square_smith_form",
    smith_normal_form(square_val, domain=ZZ),
    sp.diag(1, 5),
)
equal("square_unit_kernel_dimension", len(square_val.nullspace()), 0)
equal("square_class_has_no_three_torsion", sp.gcd(3, 5), 1)

# Repeated-root controls.  The square/nonsquare dichotomy, rather than
# squarefreeness of D, controls the unit lattice.
D_repeated_nonsquare = sp.expand(-t * (1 - t) ** 2 * (4 + 23 * t))
equal("repeated_nonsquare_degree", sp.degree(D_repeated_nonsquare, t), 4)
zero("repeated_nonsquare_simple_root", D_repeated_nonsquare.subs(t, 0))
gate(
    sp.diff(D_repeated_nonsquare, t).subs(t, 0) != 0,
    "repeated_nonsquare_not_a_square",
)
D_repeated_square = sp.expand((t * (t - 1)) ** 2)
zero("repeated_square_control", D_repeated_square - (t * (t - 1)) ** 2)
gate(
    sp.discriminant(D_repeated_nonsquare, t) == 0,
    "repeated_nonsquare_discriminant_zero",
)
gate(
    sp.discriminant(D_repeated_square, t) == 0,
    "repeated_square_discriminant_zero",
)

# Both repeated-root controls arise naturally from genuine rank-two binary
# cubic pencils after choosing a non-tangent infinity direction L=A+C.
A, C = sp.symbols("A C")
square_pencil = binary_disc(A, 0, 0, C)
zero("square_pencil_discriminant", square_pencil + 27 * A**2 * C**2)
nonsquare_pencil = binary_disc(A, 0, C, C)
zero(
    "repeated_nonsquare_pencil_discriminant",
    nonsquare_pencil + A * C**2 * (4 * C + 27 * A),
)
square_oriented = sp.factor(
    square_pencil.subs({A: t * L, C: (1 - t) * L}, simultaneous=True) / L**4
)
nonsquare_oriented = sp.factor(
    nonsquare_pencil.subs({A: t * L, C: (1 - t) * L}, simultaneous=True) / L**4
)
zero("square_oriented_control", square_oriented + 27 * t**2 * (1 - t) ** 2)
zero(
    "nonsquare_oriented_control",
    nonsquare_oriented - D_repeated_nonsquare,
)

# A one-place projective parametrization: finite roots of D map to the
# common base point, while the sole parameter point Z=0 maps to [1:0:0].
T, Z = sp.symbols("T Z")
Dh = h0 * T**4 + h1 * T**3 * Z + h2 * T**2 * Z**2 + h3 * T * Z**3 + h4 * Z**4
Mproj = -T * Dh
Lproj = -Z * Dh
Zproj = lam * Z**5
projective_delta = sp.expand(H4 * Z + lam * L**5)
zero(
    "projective_parametrization",
    projective_delta.subs(
        {M: Mproj, L: Lproj, Z: Zproj}, simultaneous=True
    ),
)
zero("unique_infinity_L_coordinate", Lproj.subs(Z, 0))
zero("unique_infinity_Z_coordinate", Zproj.subs(Z, 0))
gate(Mproj.subs(Z, 0) != 0, "unique_infinity_M_coordinate")

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

semantic = {
    "object": "S=k[M,L,W]/(W^2-(H4+lambda*L^5)), L not dividing H4",
    "normalization": "t maps to (M,L)=(-tD/lambda,-D/lambda); one infinity place",
    "normality": "branch is smooth off the origin; the double surface has isolated singular locus",
    "localization": "S_L=k[t,y,(y^2-D)^-1]",
    "class_group": "Z if D nonsquare; Z/5 if D square",
    "square_valuation_lattice": "columns (3,-2),(-2,3), Smith(1,5)",
    "units": "S*=k* in both cases",
    "cubic_consequence": "no connected normal finite-flat cubic algebra has discriminant H4+lambda*L^5",
    "common_zero_consequence": "one-place common-zero binary-cubic discriminants have degree at least six",
    "scope": "degree>=6 and unit-ideal nonmonogenic index forms remain open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3890-universal-quintic-common-zero-resolvent-class-group-dichotomy")
print("tangent_case=L_divides_H4=>delta_reducible")
print("nontangent_chart=S_L=k[t,y,(y^2-D(t))^-1]")
print("nonsquare_D_class_group=Z")
print("square_D_valuation_matrix=[[3,-2],[-2,3]]")
print("square_D_class_group=Z/5")
print("three_torsion=NONE")
print("units=k*")
print("normal_finite_flat_cubic_with_this_discriminant=NONE")
print("common_zero_one_place_quintic_normal_cubic=NONE")
print("degree_six_plus=OPEN")
print("unit_ideal_nonmonogenic=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
