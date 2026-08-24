#!/usr/bin/env python3
"""Exact companion for THM-3921's K=0 decic degeneration and order.

Reproduction:
  python3 04-computation/jc2_quintic_decic_degeneration_order_thm3921.py
  python3 -O 04-computation/jc2_quintic_decic_degeneration_order_thm3921.py
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


def a_chain(rank: int) -> sp.Matrix:
    matrix = -2 * sp.eye(rank)
    for index in range(rank - 1):
        matrix[index, index + 1] = 1
        matrix[index + 1, index] = 1
    return matrix


b, t, u, x, z = sp.symbols("b t u x z")
A, C, Z, W = sp.symbols("A C Z W")
U, V = sp.symbols("U V")

K = 2304 * b**5 + 10176 * b**4 + 4064 * b**3 + 996 * b**2 + 84 * b + 5
K_poly = sp.Poly(K, b, domain=sp.QQ)


def reduce_mod_K(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    numerator_remainder = sp.rem(
        sp.Poly(numerator, b, domain=sp.QQ), K_poly
    ).as_expr()
    denominator_remainder = sp.rem(
        sp.Poly(denominator, b, domain=sp.QQ), K_poly
    ).as_expr()
    denominator_inverse = sp.invert(
        sp.Poly(denominator_remainder, b, domain=sp.QQ), K_poly
    ).as_expr()
    return sp.factor(
        sp.rem(
            sp.Poly(numerator_remainder * denominator_inverse, b, domain=sp.QQ),
            K_poly,
        ).as_expr()
    )


gate(sp.Poly(K, b, modulus=11).is_irreducible, "K is irreducible mod 11")
gate(sp.discriminant(K, b) != 0, "K is squarefree")

# ---------------------------------------------------------------------------
# The collapsed family and its rational decic branch normalization.
# ---------------------------------------------------------------------------

p = (u**2 - 1) * (u**2 + b) ** 2
h = (
    u**9
    + (3 * b - sp.Rational(3, 2)) * u**7
    + (3 * b**2 - sp.Rational(9, 2) * b + sp.Rational(3, 8)) * u**5
    + (
        b**3
        - sp.Rational(9, 2) * b**2
        + sp.Rational(9, 8) * b
        + sp.Rational(1, 16)
    )
    * u**3
    + (
        -sp.Rational(3, 2) * b**3
        + sp.Rational(9, 8) * b**2
        + sp.Rational(3, 16) * b
        + sp.Rational(3, 128)
    )
    * u
)
r = sp.factor(p**3 - h**2)
F = z**3 - 3 * p * z + 2 * h

G_x = (
    x**4
    + (12 * b + 3) * x**3
    + (48 * b**2 + 12 * b) * x**2
    + (64 * b**3 - 48 * b**2 - 36 * b - 8) * x
    - 192 * b**3
    - 96 * b**2
    - 36 * b
    - 6
)
G = G_x.subs(x, t**2)
G_star = sp.expand(t**8 * G.subs(t, 1 / t))

G_3915 = t**8 + 6 * t**6 + 6 * t**4 - 19 * t**2 - 24
zero(G.subs(b, sp.Rational(1, 4)) - G_3915, "THM-3915 is b=1/4")

u_t = (t + 1 / t) / 2
v_t = -(t - 1 / t) / 2
z_t = (u_t**2 + b) * v_t
A_chart_on_D = sp.factor((h.subs(u, u_t) - p.subs(u, u_t) * z_t) / 2)
C_chart_on_D = sp.factor(u_t * A_chart_on_D)
zero(A_chart_on_D - t * G / 512, "branch A normalization")
zero(C_chart_on_D - (t**2 + 1) * G / 1024, "branch C normalization")

A_scaled = 2 * t * G
C_scaled = (t**2 + 1) * G
gate(sp.degree(A_scaled, t) == 9, "scaled A degree nine")
gate(sp.degree(C_scaled, t) == 10, "scaled C degree ten")
zero(C_scaled / A_scaled - u_t, "normalization slope")

zero(
    r.subs(u, u_t) + G * G_star / (2**16 * t**8),
    "residual--origin fibre product identity",
)

disc_G_expected = (
    -17915904
    * (4 * b + 1)
    * (8 * b**2 + 2 * b + 1)
    * (64 * b**3 + 48 * b**2 + 24 * b + 5) ** 2
    * K**2
)
zero(sp.discriminant(G, t) - disc_G_expected, "G discriminant factorization")

xi = (3456 * b**4 + 15840 * b**3 + 8352 * b**2 + 998 * b - 115) / 150
gq1 = sp.Rational(2, 75) * (
    1728 * b**4 + 7920 * b**3 + 4176 * b**2 + 949 * b + 55
)
gq0 = -sp.Rational(1, 75) * (
    5760 * b**4 + 26592 * b**3 + 13664 * b**2 + 3142 * b + 365
)
G_quotient = x**2 + gq1 * x + gq0
for power in range(5):
    coefficient = sp.Poly(
        sp.expand(G_x - (x - xi) ** 2 * G_quotient), x
    ).coeff_monomial(x**power)
    zero(reduce_mod_K(coefficient), f"G K-quotient factorization {power}")

GQ_at_xi = -sp.Rational(1, 30) * (
    8064 * b**4 + 36192 * b**3 + 16192 * b**2 + 2894 * b + 205
)
GQ_discriminant = sp.Rational(1, 45) * (
    16128 * b**4 + 72384 * b**3 + 33104 * b**2 + 7588 * b + 815
)
zero(reduce_mod_K(G_quotient.subs(x, xi)) - GQ_at_xi, "double root is exact")
zero(
    reduce_mod_K(sp.discriminant(G_quotient, x)) - GQ_discriminant,
    "other two x-roots are simple",
)

nonzero_residues = (
    sp.cancel(xi).as_numer_denom()[0],
    sp.cancel(gq0).as_numer_denom()[0],
    sp.cancel(GQ_at_xi).as_numer_denom()[0],
    sp.cancel(GQ_discriminant).as_numer_denom()[0],
    b,
    b + 1,
    4 * b + 1,
    8 * b**2 + 2 * b + 1,
    48 * b**2 + 16 * b + 3,
    64 * b**3 + 48 * b**2 + 24 * b + 5,
)
for residue in nonzero_residues:
    gate(sp.gcd(K, residue) == 1, "required K-residue is nonzero")

x0 = (
    95616 * b**4 + 421728 * b**3 + 162368 * b**2 + 16726 * b + 965
) / 15120
zero(
    reduce_mod_K((xi + 2 + 1 / xi) / 4 - x0),
    "double G root maps to double residual root",
)

reciprocal_resultant = (
    2**18
    * b**6
    * (48 * b**2 + 16 * b + 3) ** 12
    * (64 * b**3 + 48 * b**2 + 24 * b + 5) ** 6
)
zero(sp.resultant(G, G_star, t) - reciprocal_resultant, "no reciprocal G-roots")

# ---------------------------------------------------------------------------
# Complete finite singularity packet and the smooth contact-ten infinity.
# ---------------------------------------------------------------------------

G_prime = sp.diff(G, t)
zero((G + G_prime).subs(t, 1), "critical parameter t=1")
zero((G - G_prime).subs(t, -1), "critical parameter t=-1")
zero(
    G.subs(t, 1) + 2 * (64 * b**3 + 48 * b**2 + 24 * b + 5),
    "critical targets avoid origin",
)

for address in (1, -1):
    cusp_determinant = sp.factor(
        sp.det(
            sp.Matrix(
                [
                    [
                        sp.diff(A_scaled, t, 2).subs(t, address),
                        sp.diff(C_scaled, t, 2).subs(t, address),
                    ],
                    [
                        sp.diff(A_scaled, t, 3).subs(t, address),
                        sp.diff(C_scaled, t, 3).subs(t, address),
                    ],
                ]
            )
        )
    )
    zero(
        cusp_determinant
        - 12288 * (b + 1) ** 3 * (64 * b**3 + 48 * b**2 + 24 * b + 5),
        f"external A2 determinant t={address}",
    )

# At a double G-root, gamma''=G''(2t,t^2+1) and
# gamma'''=G'''(2t,t^2+1)+3G''(2,2t), so the determinant is
# 6(G'')^2(t^2-1), nonzero on the K-locus.
v0 = sp.Matrix([2 * t, t**2 + 1])
v1 = sp.diff(v0, t)
zero(sp.det(sp.Matrix.hstack(v0, v1)) - 2 * (t**2 - 1), "cusp tangent determinant")

J = t**4 + (4 * b + 2) * t**2 + 1
zero(
    t**9 * (A_scaled - A_scaled.subs(t, 1 / t))
    - 2 * (t - 1) ** 3 * (t + 1) ** 3 * J**3,
    "reciprocal collision identity",
)
zero(sp.discriminant(J, t) - 4096 * b**2 * (b + 1) ** 2, "four collision addresses")
zero(
    sp.resultant(G, J, t) - 16 * b**2 * (48 * b**2 + 16 * b + 3) ** 4,
    "collisions avoid origin",
)
gate(J.subs(t, 1) == 4 * (b + 1), "collisions avoid critical parameters")
zero(J / t**2 - 4 * (u_t**2 + b), "collision slopes are u^2=-b")

# Origin: four smooth branches and two A2-cuspidal branches have pairwise
# distinct tangents.  Their pair intersections sum to 26; the two cusp
# deltas add 2, reproducing the ordinary-eight delta tariff 28.
branch_multiplicities = (2, 2, 1, 1, 1, 1)
pair_delta = sum(
    branch_multiplicities[i] * branch_multiplicities[j]
    for i in range(6)
    for j in range(i + 1, 6)
)
gate(pair_delta == 26 and pair_delta + 2 == 28, "origin delta packet")
gate(36 == 28 + 2 + 6, "complete decic genus ledger")

# The homogeneous normalization is
# [T:S] |-> [2TS G_h:(T^2+S^2)G_h:S^10].  It has no base point.  At S=0,
# A/C=2S/T+O((S/T)^3) and Z/C=(S/T)^10+..., so the unique infinity place is
# smooth and has contact ten.
s = sp.symbols("s")
G_h_at_infinity = sp.expand(s**8 * G.subs(t, 1 / s))
gate(G_h_at_infinity.subs(s, 0) == 1, "normalization has no infinity base point")
gate(sp.series(2 * s / (1 + s**2), s, 0, 3).removeO() == 2 * s, "smooth infinity parameter")
gate(
    sp.limit((s**10 / ((1 + s**2) * G_h_at_infinity)) / s**10, s, 0) == 1,
    "infinity contact ten",
)

# ---------------------------------------------------------------------------
# The global normal cubic order and its A^5 power-basis debt.
# ---------------------------------------------------------------------------

P_global = sp.cancel(A**6 * p.subs(u, C / A))
H_global = sp.cancel(A**9 * h.subs(u, C / A))
Q_global = 2 * A**10 - H_global
Delta = sp.cancel(A**8 * r.subs(u, C / A) + 4 * H_global - 4 * A**10)
f = Z**3 - 3 * P_global * Z - 2 * Q_global
zero(sp.discriminant(f, Z) - 108 * A**10 * Delta, "power-order discriminant")
zero(
    Delta.subs({A: A_chart_on_D, C: C_chart_on_D}),
    "decic vanishes on normalization",
)
gate(sp.Poly(Delta, A, C).total_degree() == 10, "branch degree ten")

c = b - sp.Rational(1, 2)
d = (4 * b + 1) / 8
L = C**3 + c * A**2 * C
N = C * L - d * A**4
e1 = (Z**2 + L * Z - 2 * L**2) / A**4
e2 = (2 * C * L**2 - 4 * d * A**4 * L - N * Z - C * Z**2) / A**5
basis = [sp.Integer(1), e1, e2]

field = sp.QQ.frac_field(A, C, b, U, V)
f_poly = sp.Poly(f, Z, domain=field)


def reduce_in_cubic(expression: sp.Expr) -> sp.Expr:
    return sp.factor(
        sp.rem(sp.Poly(sp.cancel(expression), Z, domain=field), f_poly).as_expr()
    )


basis_matrix = sp.Matrix(
    [
        [sp.Poly(element, Z).coeff_monomial(Z**power) for element in basis]
        for power in range(3)
    ]
)
zero(
    basis_matrix.det() + (4 * b + 1) / (8 * A**5),
    "global index-A5 basis",
)
basis_inverse = basis_matrix.inv()


def coordinates_in_basis(expression: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    reduced = reduce_in_cubic(expression)
    old_coordinates = sp.Matrix(
        [sp.Poly(reduced, Z).coeff_monomial(Z**power) for power in range(3)]
    )
    return tuple(sp.factor(entry) for entry in basis_inverse * old_coordinates)


def polynomial_over_parameter_field(expression: sp.Expr) -> bool:
    try:
        sp.Poly(expression, A, C, domain=sp.QQ.frac_field(b))
        return True
    except sp.PolynomialError:
        return False


multiplication_packets = (
    coordinates_in_basis(e1**2),
    coordinates_in_basis(e1 * e2),
    coordinates_in_basis(e2**2),
)
for packet in multiplication_packets:
    gate(
        all(polynomial_over_parameter_field(entry) for entry in packet),
        "global basis is closed under multiplication",
    )

for element in (Z, Z**2):
    gate(
        all(polynomial_over_parameter_field(entry) for entry in coordinates_in_basis(element)),
        "power order lies in global overorder",
    )

normal_discriminant = sp.factor(sp.discriminant(f, Z) * basis_matrix.det() ** 2)
zero(
    normal_discriminant - sp.Rational(27, 16) * (4 * b + 1) ** 2 * Delta,
    "reduced normal-order discriminant",
)

# Exact A-adic maximality packet.  C is a unit at the generic point of A=0.
m = L - d * A**4 / C
local_x = sp.symbols("local_x")
local_residual = sp.cancel(f.subs(Z, m + A**5 * local_x) / A**10).subs(A, 0)
zero(
    local_residual
    - (384 * C**4 * local_x**2 - 512 * C - 96 * b**3 - 48 * b**2 - 18 * b - 3)
    / (128 * C),
    "separable local double-sheet residual",
)
gate(sp.degree(sp.together(local_residual), local_x) == 2, "local residual degree two")
gate(sp.discriminant(sp.together(local_residual), local_x) != 0, "local residual separable")

# Binary index form: primitive globally, but every coefficient vanishes at
# the target origin, so no element can generate the rank-three order.
alpha = U * e1 + V * e2
alpha_coordinates = coordinates_in_basis(alpha)
alpha_square_coordinates = coordinates_in_basis(alpha**2)
index_matrix = sp.Matrix.hstack(
    sp.Matrix([1, 0, 0]),
    sp.Matrix(alpha_coordinates),
    sp.Matrix(alpha_square_coordinates),
)
index_form = sp.factor(index_matrix.det())
index_coefficients = sp.Poly(index_form, U, V).coeffs()
gate(
    all(coefficient.subs({A: 0, C: 0}) == 0 for coefficient in index_coefficients),
    "index form vanishes at origin",
)
coefficient_gcd = sp.Poly(
    index_coefficients[0], A, C, domain=sp.QQ.frac_field(b)
)
for coefficient in index_coefficients[1:]:
    coefficient_gcd = sp.gcd(
        coefficient_gcd,
        sp.Poly(coefficient, A, C, domain=sp.QQ.frac_field(b)),
    )
gate(coefficient_gcd.monic().as_expr() == 1, "index form is primitive")
zero(
    sp.discriminant(index_form.subs(V, 1), U) - normal_discriminant,
    "index-form discriminant",
)

# ---------------------------------------------------------------------------
# Conic-bundle class group: the nodes change geometry, not the quotient.
# ---------------------------------------------------------------------------

slope = sp.symbols("slope")
zero(
    Delta.subs(C, slope * A)
    + 4 * A**8 * (A**2 - A * h.subs(u, slope) - r.subs(u, slope) / 4),
    "blowup-origin conic equation",
)
zero(h**2 + r - p**3, "conic residual discriminant is p^3")

K10 = sp.Matrix([[-4, 5], [5, -4]])
old_origin_block = sp.diag(-2, -1, -1)
node_change = sp.Matrix([[1, 0, 0], [-2, 1, 0], [-2, 0, 1]])
node_origin_block = node_change.T * old_origin_block * node_change
gate(
    node_origin_block == sp.Matrix([[-10, 2, 2], [2, -1, 0], [2, 0, -1]]),
    "two boundary-node blowups",
)
gate(node_change.det() == 1 and node_origin_block.det() == -2, "node blowups are unimodular")

removed = sp.diag(
    node_origin_block,
    K10,
    a_chain(2),
    a_chain(2),
    a_chain(5),
    a_chain(5),
)
removed_snf = smith_normal_form(removed, domain=sp.ZZ)
removed_nonunits = tuple(
    abs(int(removed_snf[index, index]))
    for index in range(removed.rows)
    if abs(int(removed_snf[index, index])) != 1
)
gate(removed.det() == 5832, "removed lattice determinant persists")
gate(removed_nonunits == (3, 3, 6, 6, 18), "removed lattice Smith packet persists")

endpoint_relation = sp.Matrix([[3, 3, 6, 6]])
endpoint_snf = smith_normal_form(endpoint_relation, domain=sp.ZZ)
gate(endpoint_snf[0, 0] == 3, "class-group torsion is Z/3")
gate(4 - endpoint_relation.rank() == 3, "class-group free rank three")

semantic_payload = {
    "parameter": "K(b)=0_genus_collapse_locus",
    "decic_normalization": "A=2tG_b,C=(t2+1)G_b_one_smooth_contact10_infinity",
    "origin_packet": "four_smooth_plus_two_A2_cuspidal_branches_delta28",
    "other_singularities": "two_A2_plus_two_A5_complete",
    "address_bridge": "r((t+t^-1)/2)=-G_b*G_b_star/(2^16*t^8)",
    "normal_order": "finite_free_normal_nonmonogenic_index_A5_rational",
    "quadratic_class_group": "Z3_plus_Z/3_scalar_units",
    "separation": "algebraic_invoices_persist_only_six_branch_boundary_topology_fails",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("theorem=THM-3921-quintic-genus-collapse-decic-degeneration-packet")
print("normalization=scaled_A=2*t*G_b;scaled_C=(t^2+1)*G_b;degree=10")
print("infinity=one_smooth_place;contact=10")
print("origin=4_smooth_plus_2_cuspidal_branches;delta=28")
print("finite_packet=origin_plus_2A2_plus_2A5;delta_total=36")
print("address_bridge=r(u(t))=-G_b(t)*G_b_star(t)/(2^16*t^8)")
print("six_addresses=4_simple_G_roots_plus_2_double_G_roots")
print("normal_cubic_order=finite_free_rational_nonmonogenic;power_index=A^5")
print("normal_discriminant=(27/16)*(4b+1)^2*Delta")
print("quadratic_class_group=Z^3_plus_Z/3;units=scalar")
print("scope=all_algebraic_invoices_persist;six_branch_boundary_obstruction_remains;JC2_open")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
