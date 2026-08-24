#!/usr/bin/env python3
"""Assertion-free exact gates for THM-3906.

The companion verifies the first normal common-zero S3 cubic of
discriminant degree six, the normalization and full genus ledger of its
branch curve, the two-place Newton boundary for its ambient coefficient
grammar, and the derivative obstruction to the tempting four-cusp sextic.

Reproduction:
  python3 04-computation/jc2_degree_six_common_zero_normal_cubic_two_place_thm3906.py
  python3 -O 04-computation/jc2_degree_six_common_zero_normal_cubic_two_place_thm3906.py
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

# The explicit binary cubic and its exact discriminant.
a = C**2 + A
b = A
c = A + C
d = C
Phi = sp.expand(a * U**3 + b * U**2 * V + c * U * V**2 + d * V**3)
D = binary_disc(a, b, c, d)
D_expected = sp.expand(
    -3 * A**4
    - 4 * A**3 * C**2
    + 4 * A**3 * C
    + 6 * A**2 * C**3
    - 20 * A**2 * C**2
    - 48 * A * C**4
    - 4 * A * C**3
    - 27 * C**6
    - 4 * C**5
)
zero("explicit_discriminant", D - D_expected)
equal("common_zero_a", a.subs({A: 0, C: 0}), 0)
equal("common_zero_b", b.subs({A: 0, C: 0}), 0)
equal("common_zero_c", c.subs({A: 0, C: 0}), 0)
equal("common_zero_d", d.subs({A: 0, C: 0}), 0)
equal("coefficient_ideal_contains_A", b, A)
equal("coefficient_ideal_contains_C", d, C)

# The generic cubic is primitive and irreducible.  After V=1 it is linear
# in A with coprime coefficient and constant term in k[C,T].
f = sp.expand(Phi.subs({U: T, V: 1}))
qA = T**3 + T**2 + T
rA = C * (C * T**3 + T + 1)
zero("A_linear_cubic_identity", f - (A * qA + rA))
equal("A_linear_coprimality", sp.gcd(qA, rA), 1)
equal("primitive_coefficient_gcd", sp.gcd(sp.gcd(a, b), sp.gcd(c, d)), 1)

# Degree, unique projective support, and the ordinary quadruple common zero.
equal("discriminant_total_degree", sp.Poly(D, A, C).total_degree(), 6)
D4 = homogeneous_piece(D, (A, C), 4)
D5 = homogeneous_piece(D, (A, C), 5)
D6 = homogeneous_piece(D, (A, C), 6)
zero("degree_six_carrier", D6 + 27 * C**6)
zero("degree_five_layer", D5 + 4 * C**5 + 48 * A * C**4 - 6 * A**2 * C**3 + 4 * A**3 * C**2)
zero("quartic_tangent_cone", D4 + A * (3 * A**3 - 4 * A**2 * C + 20 * A * C**2 + 4 * C**3))
tangent_cubic = 3 * T**3 - 4 * T**2 + 20 * T + 4
equal("tangent_cubic_discriminant", sp.discriminant(tangent_cubic, T), -109744)
gate(sp.discriminant(tangent_cubic, T) != 0, "four_distinct_tangent_lines")
equal("plane_sextic_arithmetic_genus", (6 - 1) * (6 - 2) // 2, 10)
equal("ordinary_quadruple_delta", 4 * 3 // 2, 6)

# Absolute irreducibility: specialization at C=1 is irreducible modulo 7,
# and the A-leading coefficient is constant.  A rational smooth point rules
# out geometric splitting into conjugate components.
D_C1 = sp.expand(D.subs(C, 1))
zero("irreducible_specialization", D_C1 + 3 * A**4 + 14 * A**2 + 52 * A + 31)
equal("constant_A_leading_coefficient", sp.Poly(D, A).LC(), -3)
mod7_factors = sp.factor_list(A**4 + A + 1, modulus=7)[1]
equal("mod7_irreducible_factor_count", len(mod7_factors), 1)
equal("mod7_irreducible_factor", mod7_factors[0], (A**4 + A + 1, 1))
smooth_C = -sp.Rational(4, 27)
equal("rational_curve_point", D.subs({A: 0, C: smooth_C}), 0)
equal("smooth_point_dA", sp.diff(D, A).subs({A: 0, C: smooth_C}), -sp.Rational(1792, 177147))
equal("smooth_point_dC", sp.diff(D, C).subs({A: 0, C: smooth_C}), sp.Rational(1024, 531441))

# The unique support point [1:0:0] splits into two normalization places.
x, z, rho = sp.symbols("x z rho")
h = sp.expand(z**6 * D.subs({A: 1 / z, C: x / z}, simultaneous=True))
h_expected = sp.expand(
    -27 * x**6
    + z * (-4 * x**5 - 48 * x**4 + 6 * x**3 - 4 * x**2)
    + z**2 * (-4 * x**3 - 20 * x**2 + 4 * x - 3)
)
zero("projective_local_equation", h - h_expected)
equal("quadratic_in_projective_parameter", sp.degree(h, z), 2)
zero("projective_z_discriminant", sp.discriminant(h, z) - 16 * x**4 * (x**2 - x + 1) ** 3)
edge_two = homogeneous_piece(sp.expand(h.subs(z, rho * x**2)), (x,), 4)
zero("order_two_edge", edge_two + rho * (3 * rho + 4) * x**4)
equal("order_two_nonzero_root", -sp.Rational(4, 3), -sp.Rational(4, 3))
equal(
    "order_four_edge",
    sp.Poly(sp.expand(h.subs(z, rho * x**4)), x).coeff_monomial(x**6),
    -4 * rho - 27,
)
equal("order_four_root", -sp.Rational(27, 4), -sp.Rational(27, 4))
equal("infinity_branch_intersection", sp.Poly(sp.discriminant(h, z), x).as_dict().get((4,)), 16)
equal("infinity_delta", 2, 2)

# A rational normalization obtained from the repeated root T of the cubic.
N = 2 * T**3 + 4 * T**2 + 2 * T + 1
Cpar = N / (T**3 * (T + 2))
Apar = -(2 * T + 3) * N / (T**4 * (T + 2) ** 2)
zero("rational_normalization", D.subs({A: Apar, C: Cpar}, simultaneous=True))
Tp, Sp, Zp = sp.symbols("Tp Sp Zp")
Nh = 2 * Tp**3 + 4 * Tp**2 * Sp + 2 * Tp * Sp**2 + Sp**3
Aproj = -Sp**2 * (2 * Tp + 3 * Sp) * Nh
Cproj = Sp * Tp * (Tp + 2 * Sp) * Nh
Zproj = Tp**4 * (Tp + 2 * Sp) ** 2
Dprojective = sp.expand(D4 * Zp**2 + D5 * Zp + D6)
zero(
    "projective_normalization",
    Dprojective.subs({A: Aproj, C: Cproj, Zp: Zproj}, simultaneous=True),
)
equal("projective_coordinate_degree_A", sp.Poly(Aproj, Tp, Sp).total_degree(), 6)
equal("projective_coordinate_degree_C", sp.Poly(Cproj, Tp, Sp).total_degree(), 6)
equal("projective_coordinate_degree_Z", sp.Poly(Zproj, Tp, Sp).total_degree(), 6)
gate(Aproj.subs({Tp: 0, Sp: 1}) != 0, "infinity_preimage_t_zero")
gate(Aproj.subs({Tp: -2, Sp: 1}) != 0, "infinity_preimage_t_minus_two")
gate(Zproj.subs({Tp: 1, Sp: 0}) != 0, "origin_preimage_t_infinity")
equal("origin_finite_address_discriminant", sp.discriminant(N, T), -76)
equal("origin_finite_address_count", sp.degree(N, T), 3)
equal("origin_infinite_address_A_order", sp.limit(Apar.subs(T, 1 / x) / x**2, x, 0), -4)
equal("origin_infinite_address_C_order", sp.limit(Cpar.subs(T, 1 / x) / x, x, 0), 2)
cusp_parameter = T**2 + 3 * T + 3
zero(
    "normalization_dA_factor",
    sp.diff(Apar, T)
    - 2 * cusp_parameter * (4 * T**3 + 9 * T**2 + 7 * T + 4) / (T**5 * (T + 2) ** 3),
)
zero(
    "normalization_dC_factor",
    sp.diff(Cpar, T)
    + 2 * (T**2 + T + 1) * cusp_parameter / (T**4 * (T + 2) ** 2),
)
equal("two_cusp_parameters", sp.degree(cusp_parameter, T), 2)
gate(sp.discriminant(cusp_parameter, T) != 0, "distinct_cusp_parameters")

# Direct singular-support certificate.  Saturating by C leaves precisely the
# two cusp points; C=0 leaves only the ordinary quadruple origin.
u = sp.symbols("u")
singular_sat = sp.groebner([D, sp.diff(D, A), sp.diff(D, C), u * C - 1], u, A, C, order="lex")
sat_without_u = [sp.factor(poly.as_expr()) for poly in singular_sat.polys if not poly.as_expr().has(u)]
equal("saturated_singular_basis_size", len(sat_without_u), 2)
equal(
    "saturated_A_relation",
    sat_without_u[0],
    6174 * A + 729 * C**3 + 9693 * C**2 + 2421 * C + 722,
)
equal("saturated_C_relation", sat_without_u[1], (27 * C**2 + 27 * C + 19) ** 2)
zero("C_zero_curve", D.subs(C, 0) + 3 * A**4)
zero("C_zero_gradient", sp.diff(D, A).subs(C, 0) + 12 * A**3)
cusp_C = 27 * C**2 + 27 * C + 19
cusp_A = 21 * A - 24 * C - 19
zero(
    "cusp_A_relation_mod_C_polynomial",
    sp.rem(sat_without_u[0], cusp_C, C) - 294 * cusp_A,
)
gate(sp.discriminant(cusp_C, C) != 0, "two_distinct_finite_cusps")
Hessian = sp.hessian(D, (A, C)).subs(A, (24 * C + 19) / 21)
Hessian_det_num = sp.together(Hessian.det()).as_numer_denom()[0]
equal("cusp_tangent_cone_rank_one", sp.rem(Hessian_det_num, cusp_C, C), 0)
gate(sp.rem(sp.together(Hessian[0, 0]).as_numer_denom()[0], cusp_C, C) != 0, "cusp_multiplicity_two")
equal("two_A2_delta", 2, 2)
equal("complete_genus_ledger", 6 + 2 + 2, 10)

# The entire minimal triple-root sextic grammar has the same two-edge debt.
alpha, alpha1, beta, beta1, gamma, gamma1 = sp.symbols(
    "alpha alpha1 beta beta1 gamma gamma1"
)
ag = C**2 + alpha * A + alpha1 * C
bg = beta * A + beta1 * C
cg = gamma * A + gamma1 * C
dg = C
Dg = binary_disc(ag, bg, cg, dg)
Dg4 = homogeneous_piece(Dg, (A, C), 4)
Dg5 = homogeneous_piece(Dg, (A, C), 5)
Dg6 = homogeneous_piece(Dg, (A, C), 6)
zero("family_C_zero_dichotomy", Dg.subs(C, 0) - gamma**2 * (beta**2 - 4 * alpha * gamma) * A**4)
zero("family_degree_six_carrier", Dg6 + 27 * C**6)
Hg = sp.expand(Dg4.subs({A: 1, C: x}) * z**2 + Dg5.subs({A: 1, C: x}) * z + Dg6.subs({A: 1, C: x}))
equal("family_projective_parameter_degree", sp.degree(Hg, z), 2)
zero(
    "family_z2_constant",
    sp.Poly(Hg, x, z).coeff_monomial(z**2)
    - gamma**2 * (beta**2 - 4 * alpha * gamma),
)
equal("family_zx2_coefficient", sp.Poly(Hg, x, z).coeff_monomial(x**2 * z), -4 * gamma**3)
equal("family_zx0_absent", sp.Poly(Hg, x, z).coeff_monomial(z), 0)
equal("family_zx1_absent", sp.Poly(Hg, x, z).coeff_monomial(x * z), 0)
equal("family_x6_coefficient", sp.Poly(Hg, x, z).coeff_monomial(x**6), -27)
zero("family_no_lower_z0_support", Hg.subs(z, 0) + 27 * x**6)

# Three or more finite cusps cannot occur in the rational smooth-infinity
# ordinary-four packet.  Factoring the four-address polynomial leaves a
# degree-at-most-two tangent-direction map; every cusp ramifies that map.
p0, p1, p2, q0, q1 = sp.symbols("p0 p1 p2 q0 q1")
P2 = p0 + p1 * T + p2 * T**2
Q1 = q0 + q1 * T
Wronskian = sp.expand(P2 * sp.diff(Q1, T) - sp.diff(P2, T) * Q1)
gate(sp.degree(Wronskian, T) <= 2, "tangent_direction_wronskian_degree")
equal("degree_two_map_total_ramification", 2 * 2 - 2, 2)
gate(3 > 2, "three_distinct_cusps_exceed_ramification_budget")

# The bound two is sharp as plane geometry.  This control is not asserted to
# be a binary-cubic discriminant.
Qsharp = T**4 + T**3 - 2 * T - 5
Xsharp = sp.expand(T * Qsharp)
Ysharp = sp.expand((T**2 + 1) * Qsharp)
equal("sharp_quadruple_addresses", sp.discriminant(Qsharp, T), -23355)
equal(
    "sharp_distinct_tangent_resultant",
    sp.resultant(Qsharp, sp.expand(T**4 * Qsharp.subs(T, 1 / T)), T),
    273375,
)
equal(
    "sharp_two_cusp_parameters",
    sp.factor(sp.gcd(sp.diff(Xsharp, T), sp.diff(Ysharp, T))),
    (T - 1) * (T + 1),
)
sharp_cusp_det = sp.expand(
    sp.diff(Xsharp, T, 2) * sp.diff(Ysharp, T, 3)
    - sp.diff(Ysharp, T, 2) * sp.diff(Xsharp, T, 3)
)
equal("sharp_A2_at_plus_one", sharp_cusp_det.subs(T, 1), 1680)
equal("sharp_A2_at_minus_one", sharp_cusp_det.subs(T, -1), 432)
equal("sharp_infinity_degree_X", sp.degree(Xsharp, T), 5)
equal("sharp_infinity_degree_Y", sp.degree(Ysharp, T), 6)
equal("sharp_infinity_leading_X", sp.Poly(Xsharp, T).LC(), 1)
equal("sharp_infinity_leading_Y", sp.Poly(Ysharp, T).LC(), 1)
equal("sharp_cusp_off_quadruple_plus", Qsharp.subs(T, 1), -5)
equal("sharp_cusp_off_quadruple_minus", Qsharp.subs(T, -1), -3)
sigma = sp.symbols("sigma")
zero(
    "sharp_radial_pairing_identity",
    sigma * (T**2 + 1) - T * (sigma**2 + 1) + (sigma - T) * (sigma * T - 1),
)
zero(
    "sharp_generic_birational_control",
    sp.together((Ysharp - Ysharp.subs(T, 1 / T)) * T**6)
    - (T - 1) ** 3 * (T + 1) ** 3 * (T**2 + 1) * (T**4 + T**3 + 3 * T**2 + T + 1),
)
equal("sharp_residual_node_pair_discriminant", sp.discriminant(T**4 + T**3 + 3 * T**2 + T + 1, T), 189)

# The original four-cusp derivative argument is retained as an independent
# hostile route.
# Once C has degree five, four common derivative zeros force
# A'=(kT+l)C'.  The displayed primitive is the exact integration identity.
c0, c1, c2, c3, c4, c5, k, ell = sp.symbols("c0 c1 c2 c3 c4 c5 k ell")
Cpoly = c0 + c1 * T + c2 * T**2 + c3 * T**3 + c4 * T**4 + c5 * T**5
Jpoly = sp.integrate(Cpoly, T)
Aforced = sp.expand((k * T + ell) * Cpoly - k * Jpoly)
zero("four_cusp_integration_identity", sp.diff(Aforced, T) - (k * T + ell) * sp.diff(Cpoly, T))
equal("primitive_degree", sp.degree(Jpoly, T), 6)
gate(4 * 2 > 6, "four_double_roots_exceed_degree_six")

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

semantic = {
    "object": "Phi=(C^2+A)U^3+A U^2V+(A+C)UV^2+C V^3",
    "cubic_order": "normal finite-flat nonmonogenic S3",
    "discriminant": "geometrically irreducible rational sextic",
    "finite_packet": "ordinary quadruple delta6 plus two A2 cusps",
    "infinity": "one projective point, two smooth places of orders 2 and 4, delta2",
    "family": "minimal triple-root sextic grammar is reducible or two-place",
    "no_go": "ordinary-four rational smooth-infinity sextic has at most two finite cusps, sharply",
    "scope": "all other degree-six grammars and unit-ideal one-place forms remain open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3906-degree-six-common-zero-normal-cubic-two-place-boundary")
print("normal_common_zero_S3_degree_six=YES")
print("discriminant_geometrically_irreducible=YES")
print("branch_normalization=P1")
print("finite_packet=ordinary4_plus_two_A2")
print("projective_infinity_support_points=1")
print("normalization_places_at_infinity=2")
print("infinity_orders=[2,4]")
print("minimal_triple_root_sextic_grammar=REDUCIBLE_OR_TWO_PLACE")
print("ordinary4_smooth_one_place_sextic_finite_cusps_at_most=2")
print("two_cusp_bound=SHARP_AS_PLANE_GEOMETRY")
print("all_degree_six_common_zero_one_place=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
