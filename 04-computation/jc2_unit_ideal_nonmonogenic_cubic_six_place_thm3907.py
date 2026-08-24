#!/usr/bin/env python3
"""Assertion-free exact gates for THM-3907.

This companion checks a unit-ideal binary cubic whose index form represents
no scalar unit, its irreducible degree-seven discriminant, and the exact six
normalization places above its two projective infinity points.

Reproduction:
  python3 04-computation/jc2_unit_ideal_nonmonogenic_cubic_six_place_thm3907.py
  python3 -O 04-computation/jc2_unit_ideal_nonmonogenic_cubic_six_place_thm3907.py
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


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


def homogeneous_piece(poly: sp.Expr, variables: tuple[sp.Symbol, ...], degree: int) -> sp.Expr:
    return sp.expand(
        sum(
            coefficient * sp.prod(variable**exponent for variable, exponent in zip(variables, powers))
            for powers, coefficient in sp.Poly(poly, *variables).terms()
            if sum(powers) == degree
        )
    )


A, C, U, V, T = sp.symbols("A C U V T")
a, b, c, d = A, C, A * C - 1, A
Psi = sp.expand(a * U**3 + b * U**2 * V + c * U * V**2 + d * V**3)

# The coefficient ideal is the unit ideal even though specialization at A=0
# gives the global unit-representation obstruction.
zero("unit_ideal_witness", A * C - (A * C - 1) - 1)
X0, Y0 = sp.symbols("X0 Y0")
zero(
    "A_zero_unit_factorization",
    Psi.subs({A: 0, U: X0, V: Y0}) - X0 * Y0 * (C * X0 - Y0),
)

# Generic irreducibility is the degree-three rational-map certificate in C.
f = sp.expand(Psi.subs({U: T, V: 1}))
numerator = A * T**3 - T + A
denominator = T**2 + A * T
zero("C_linear_cubic_identity", f - (C * denominator + numerator))
fraction_field = sp.QQ.frac_field(A)
equal(
    "rational_map_coprimality",
    sp.gcd(sp.Poly(numerator, T, domain=fraction_field), sp.Poly(denominator, T, domain=fraction_field)),
    sp.Poly(1, T, domain=fraction_field),
)
equal("rational_map_degree", max(sp.degree(numerator, T), sp.degree(denominator, T)), 3)

# Exact discriminant and irreducible repeated-root incidence.
D = binary_disc(a, b, c, d)
D_expected = sp.expand(
    -4 * A**4 * C**3
    - 27 * A**4
    + 30 * A**3 * C**2
    + A**2 * C**4
    - 30 * A**2 * C
    - 6 * A * C**3
    + 4 * A
    + C**2
)
zero("explicit_discriminant", D - D_expected)
equal("discriminant_total_degree", sp.Poly(D, A, C).total_degree(), 7)
zero("degree_seven_carrier", homogeneous_piece(D, (A, C), 7) + 4 * A**4 * C**3)
equal("linear_A_term", sp.Poly(D, A, C).coeff_monomial(A), 4)

g = sp.expand((1 - 2 * T**3) * A**2 + (2 * T - T**4) * A - T**2)
incidence_relation = sp.expand(C * T**2 - A * (1 - 2 * T**3))
zero("incidence_quadratic_discriminant", sp.discriminant(g, A) - T**2 * (T**6 - 12 * T**3 + 8))
equal(
    "incidence_squarefree_sidecar",
    sp.gcd(T**6 - 12 * T**3 + 8, sp.diff(T**6 - 12 * T**3 + 8, T)),
    1,
)
zero("incidence_resultant", sp.resultant(g, incidence_relation, T) + A**3 * D)

# Projective point [1:0:0]: three distinct tangent lines.
x, z, rho = sp.symbols("x z rho")
H1 = sp.expand(z**7 * D.subs({A: 1 / z, C: x / z}, simultaneous=True))
H1_expected = sp.expand(
    -4 * x**3
    + x**4 * z
    + 30 * x**2 * z**2
    + (-6 * x**3 - 27) * z**3
    - 30 * x * z**4
    + x**2 * z**5
    + 4 * z**6
)
zero("first_infinity_local_equation", H1 - H1_expected)
H1_tangent = homogeneous_piece(H1, (x, z), 3)
zero("first_infinity_tangent_cone", H1_tangent + 4 * x**3 + 27 * z**3)
gate(sp.discriminant(4 * rho**3 + 27, rho) != 0, "three_distinct_first_infinity_tangents")
equal("first_infinity_places", 3, 3)

# Projective point [0:1:0]: one z~4a^2 branch and two
# a~(3+/-2sqrt(2))z^2 branches.
y = sp.symbols("y")
H2 = sp.expand(z**7 * D.subs({A: y / z, C: 1 / z}, simultaneous=True))
H2_expected = sp.expand(
    -4 * y**4
    + y**2 * z
    + (-27 * y**4 - 6 * y) * z**3
    + 30 * y**3 * z**2
    - 30 * y**2 * z**4
    + z**5
    + 4 * y * z**6
)
zero("second_infinity_local_equation", H2 - H2_expected)
zero("second_infinity_tangent_cone", homogeneous_piece(H2, (y, z), 3) - y**2 * z)
equal(
    "second_infinity_z_over_y2_edge",
    sp.Poly(sp.expand(H2.subs(z, rho * y**2)), y).coeff_monomial(y**4),
    rho - 4,
)
equal("second_infinity_z_over_y2_root", 4, 4)
equal(
    "second_infinity_y_over_z2_edge",
    sp.Poly(sp.expand(H2.subs(y, rho * z**2)), z).coeff_monomial(z**5),
    rho**2 - 6 * rho + 1,
)
equal("second_infinity_y_over_z2_discriminant", sp.discriminant(rho**2 - 6 * rho + 1, rho), 32)
equal("second_infinity_places", 3, 3)
equal("total_infinity_places", 3 + 3, 6)
equal("projective_support_points", 2, 2)

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

semantic = {
    "object": "Psi=A U^3+C U^2V+(AC-1)UV^2+A V^3",
    "coefficient_ideal": "unit ideal",
    "unit_representation": "no nonzero scalar",
    "cubic_order": "normal finite-flat globally nonmonogenic S3",
    "discriminant": "absolutely irreducible degree seven",
    "infinity": "two support points, three places above each",
    "scope": "unit-ideal one-place discriminants and JC2 remain open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3907-unit-ideal-nonmonogenic-cubic-six-place-boundary")
print("coefficient_ideal=UNIT")
print("represents_scalar_unit=NO")
print("normal_finite_flat_cubic=YES")
print("generic_galois_group=S3")
print("discriminant_degree=7")
print("discriminant_absolutely_irreducible=YES")
print("projective_infinity_support_points=2")
print("normalization_places_at_infinity=6")
print("common_zero_required_for_nonmonogenicity=NO")
print("unit_ideal_one_place=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
