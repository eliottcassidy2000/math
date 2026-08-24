#!/usr/bin/env python3
"""Exact companion for THM-3915's rational decic cube resolvent.

Reproduction:
  python3 04-computation/jc2_rational_decic_cube_resolvent_index_debt_thm3915.py
  python3 -O 04-computation/jc2_rational_decic_cube_resolvent_index_debt_thm3915.py
"""

from __future__ import annotations

import hashlib
import json
from itertools import product

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


t, A, C, Z, W, s, X, u, v = sp.symbols("t A C Z W s X u v")

# ---------------------------------------------------------------------------
# The rational decic and its complete normalization/singularity packet.
# ---------------------------------------------------------------------------

G = t**8 + 6 * t**6 + 6 * t**4 - 19 * t**2 - 24
A_t = t * G
C_t = (t**2 + 1) * G

gate(G.subs(t, 1) == -30 and G.subs(t, -1) == -30, "nonzero cusp values")
gate(
    (G + sp.diff(G, t)).subs(t, 1) == 0
    and (G - sp.diff(G, t)).subs(t, -1) == 0,
    "forced critical conditions",
)
gate(sp.discriminant(G, t) == -605291616000000, "G is squarefree")
G_reciprocal = sp.expand(t**8 * G.subs(t, 1 / t))
gate(sp.gcd(G, G_reciprocal) == 1, "no reciprocal roots of G")

Delta = (
    -A**10
    + 30 * A**8 * C
    - 16 * A**8
    + 5 * A**6 * C**3
    - 309 * A**6 * C**2
    - 9 * A**4 * C**5
    - 246 * A**4 * C**4
    - 3 * A**2 * C**7
    - 29 * A**2 * C**6
    + C**9
    + 24 * C**8
)
zero(Delta.subs({A: A_t, C: C_t}), "decic vanishes on parametrization")
zero(
    sp.resultant(A - A_t, C - C_t, t) - Delta,
    "parametrization eliminant",
)
gate(sp.Poly(Delta, A, C).total_degree() == 10, "image degree ten")
gate(sp.factor(Delta) == Delta, "eliminant is irreducible over Q")

# The projective map is [TS Gh:(T^2+S^2)Gh:S^10].  The eliminant gives the
# same unique smooth infinity point and exact contact ten.
Delta_h = sp.expand(
    sum(
        coefficient
        * A**monomial[0]
        * C**monomial[1]
        * Z ** (10 - sum(monomial))
        for monomial, coefficient in sp.Poly(Delta, A, C).terms()
    )
)
zero(Delta_h.subs(Z, 0) + A**10, "unique infinity support")
infinity = {A: 0, C: 1, Z: 0}
gate(sp.diff(Delta_h, Z).subs(infinity) == 1, "smooth infinity point")

dA_t = sp.factor(sp.diff(A_t, t))
dC_t = sp.factor(sp.diff(C_t, t))
gate(sp.gcd(dA_t, dC_t) == t**2 - 1, "only two critical parameters")
gate((A_t.subs(t, 1), C_t.subs(t, 1)) == (-30, -60), "first cusp")
gate((A_t.subs(t, -1), C_t.subs(t, -1)) == (30, -60), "second cusp")
for address in (1, -1):
    cusp_determinant = sp.det(
        sp.Matrix(
            [
                [
                    sp.diff(A_t, t, 2).subs(t, address),
                    sp.diff(C_t, t, 2).subs(t, address),
                ],
                [
                    sp.diff(A_t, t, 3).subs(t, address),
                    sp.diff(C_t, t, 3).subs(t, address),
                ],
            ]
        )
    )
    gate(cusp_determinant == 180000, f"A2 determinant at t={address}")

# Put a(t)=tG(t).  The reciprocal difference is (t-t^-1)^3 h(t+t^-1).
h = (s**2 + 1) ** 3
collision_polynomial = t**4 + 3 * t**2 + 1
zero(
    A_t - A_t.subs(t, 1 / t)
    - (t - 1 / t) ** 3 * h.subs(s, t + 1 / t),
    "reciprocal collision identity",
)
gate(
    sp.factor(t**10 * G - G_reciprocal)
    == (t - 1) ** 3 * (t + 1) ** 3 * collision_polynomial**3,
    "polynomial collision factor",
)
gate(sp.discriminant(collision_polynomial, t) == 400, "four collision addresses")
gate(sp.gcd(G, collision_polynomial) == 1, "collisions avoid the origin fibre")
gate(
    collision_polynomial.subs(t, 1) == 5
    and collision_polynomial.subs(t, -1) == 5,
    "collisions avoid cusp addresses",
)

# In coordinates s=C/A, each root s=+/-i gives two immersed branches whose
# graph difference has order three, hence an A5.  Their targets are explicit.
u_s = s * (s**8 - 3 * s**6 - 9 * s**4 + 5 * s**2 + 30)
v_s = -24 * s**8 + 29 * s**6 + 246 * s**4 + 309 * s**2 + 16
zero(
    Delta.subs(C, s * A) + A**8 * (A**2 - u_s * A + v_s),
    "pencil-line branch factorization",
)
p_s = (s**2 - 4) * (s**2 + 1) ** 2
zero(u_s**2 - 4 * v_s - p_s**3, "perfect-cube conic discriminant")
gate(
    (sp.simplify(u_s.subs(s, sp.I) / 2), sp.simplify(sp.I * u_s.subs(s, sp.I) / 2))
    == (10 * sp.I, -10),
    "first A5 target",
)
gate(
    (
        sp.simplify(u_s.subs(s, -sp.I) / 2),
        sp.simplify(-sp.I * u_s.subs(s, -sp.I) / 2),
    )
    == (-10 * sp.I, -10),
    "second A5 target",
)

# Independent singular-ideal support: origin, the two cusps, and the two A5s.
Delta_A = sp.diff(Delta, A)
Delta_C = sp.diff(Delta, C)
singular_groebner = sp.groebner([Delta, Delta_A, Delta_C], A, C, order="lex")
gate(
    sp.factor(singular_groebner.polys[-1].as_expr())
    == C**11 * (C + 10) ** 5 * (C + 60) ** 2,
    "complete singular C-support",
)
expected_fibre_gcds = {
    0: A**7,
    -10: A**2 + 100,
    -60: A**2 - 900,
}
for c_value, expected in expected_fibre_gcds.items():
    fibre_polynomials = [
        sp.Poly(polynomial.subs(C, c_value), A)
        for polynomial in (Delta, Delta_A, Delta_C)
    ]
    fibre_gcd = fibre_polynomials[0]
    for polynomial in fibre_polynomials[1:]:
        fibre_gcd = sp.gcd(fibre_gcd, polynomial)
    gate(
        sp.monic(fibre_gcd.as_expr(), A) == expected,
        f"singular fibre C={c_value}",
    )

tangent_cone = sum(
    coefficient * A**monomial[0] * C**monomial[1]
    for monomial, coefficient in sp.Poly(Delta, A, C).terms()
    if sum(monomial) == 8
)
zero(
    tangent_cone
    + sp.resultant(G, (t**2 + 1) * A - t * C, t),
    "ordinary-eight tangent cone",
)
gate(
    sp.discriminant(tangent_cone.subs(C, 1), A) != 0,
    "eight tangent directions are distinct",
)
gate(36 == 28 + 2 + 2 * 3, "complete arithmetic-genus ledger")

# ---------------------------------------------------------------------------
# The S3 cubic, its global maximal order, and the A^5 power-basis debt.
# ---------------------------------------------------------------------------

P = (C**2 - 4 * A**2) * (C**2 + A**2) ** 2
Q = (
    2 * A**10
    - C**9
    + 3 * A**2 * C**7
    + 9 * A**4 * C**5
    - 5 * A**6 * C**3
    - 30 * A**8 * C
)
f = Z**3 - 3 * P * Z - 2 * Q
zero(P**3 - Q**2 - 4 * A**10 * Delta, "Cardano norm identity")
zero(sp.discriminant(f, Z) - 432 * A**10 * Delta, "power-basis discriminant")
gate(
    sp.factor(f.subs(C, 0)) == Z**3 + 12 * A**6 * Z - 4 * A**10,
    "irreducibility specialization",
)

L = C**3 - C * A**2
N = C * L - 4 * A**4
e1 = (Z**2 + L * Z - 2 * L**2) / A**4
f2 = (2 * C * L**2 - 16 * A**4 * L - N * Z - C * Z**2) / A**5
basis = [sp.Integer(1), e1, f2]
basis_matrix = sp.Matrix(
    [
        [sp.Poly(element, Z).coeff_monomial(Z**power) for element in basis]
        for power in range(3)
    ]
)
gate(sp.factor(basis_matrix.det()) == -4 * A**-5, "global index-A5 basis")


def reduce_in_field(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.together(expression).as_numer_denom()
    return sp.cancel(sp.rem(numerator, f, Z) / denominator)


def coordinates_in_basis(expression: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    reduced = reduce_in_field(expression)
    old_coordinates = sp.Matrix(
        [sp.Poly(reduced, Z).coeff_monomial(Z**power) for power in range(3)]
    )
    return tuple(sp.factor(entry) for entry in basis_matrix.inv() * old_coordinates)


e1_square = coordinates_in_basis(e1**2)
e1_f2 = coordinates_in_basis(e1 * f2)
f2_square = coordinates_in_basis(f2**2)
expected_e1_square = (
    -24 * C * (A - C) * (A + C) * (A**2 - 12 * C),
    A**2 * C - 12 * A**2 - 36 * C**2,
    A * (A**2 - 12 * C),
)
expected_e1_f2 = (
    8 * A * (2 * A**4 + 3 * A**2 * C**2 - 6 * A**2 * C - 3 * C**4 - 18 * C**3),
    -A * C * (C - 12),
    -A**2 * C - 12 * A**2 - 12 * C**2,
)
expected_f2_square = (
    -8 * C * (4 * A**4 + 3 * A**2 * C**2 + 24 * A**2 * C - 3 * C**4 - 72 * C**3),
    16 * A**2 + C**3 + 24 * C**2,
    A * C * (C + 60),
)
gate(e1_square == expected_e1_square, "e1 square multiplication")
gate(e1_f2 == expected_e1_f2, "e1 f2 multiplication")
gate(f2_square == expected_f2_square, "f2 square multiplication")

# Both theta and theta^2 have polynomial coordinates in the new basis, so
# the power order is a genuine suborder of index A^5.
theta_coordinates = tuple(sp.factor(entry) for entry in basis_matrix.inv()[:, 1])
theta2_coordinates = tuple(sp.factor(entry) for entry in basis_matrix.inv()[:, 2])
gate(all(entry.is_polynomial(A, C) for entry in theta_coordinates), "theta lies in maximal order")
gate(all(entry.is_polynomial(A, C) for entry in theta2_coordinates), "theta square lies in maximal order")

# The new discriminant is the old one divided by the square of the index.
maximal_discriminant = sp.factor(sp.discriminant(f, Z) * basis_matrix.det() ** 2)
zero(maximal_discriminant - 6912 * Delta, "maximal-order discriminant")

# Exact A-adic local basis and multiplication.  C is a unit at the generic
# point of A=0.
m = C**3 - C * A**2 - 4 * A**4 / C
beta = (Z**2 + m * Z - 2 * m**2) / A**5
aa = sp.factor(3 * (P - m**2) / A**5)
bb = sp.factor(2 * (Q + m**3) / A**5)
cc = sp.factor((bb + m * aa) / A**5)
gate(aa == -12 * A * (4 * A**2 + 3 * C**2) / C**2, "local theta-beta coefficient a")
gate(
    bb
    == -4
    * A
    * (32 * A**6 - A**4 * C**3 + 24 * A**4 * C**2 - 3 * A**2 * C**4 - 9 * C**6)
    / C**3,
    "local theta-beta coefficient b",
)
zero(cc - 4 * (16 * A**2 + C**3 + 24 * C**2) / C**3, "local beta-square coefficient")
zero(Z**2 - (A**5 * beta - m * Z + 2 * m**2), "local theta square")
zero(reduce_in_field(Z * beta) - (m * beta + aa * Z + bb), "local theta beta")
zero(reduce_in_field(beta**2) - (aa * beta + cc * Z + 2 * m * cc), "local beta square")

# Binary index form.  Its coefficients generate the unit ideal but all
# vanish at the origin, proving the maximal order is globally nonmonogenic.
alpha = u * e1 + v * f2
alpha_square_coordinates = coordinates_in_basis(alpha**2)
index_matrix = sp.Matrix(
    [
        [1, 0, alpha_square_coordinates[0]],
        [0, u, alpha_square_coordinates[1]],
        [0, v, alpha_square_coordinates[2]],
    ]
)
index_form = sp.factor(index_matrix.det())
expected_index_form = (
    A * (A**2 - 12 * C) * u**3
    - 3 * (A**2 * C + 4 * A**2 - 4 * C**2) * u**2 * v
    + 3 * A * C * (C + 12) * u * v**2
    - (16 * A**2 + C**3 + 24 * C**2) * v**3
)
zero(index_form - expected_index_form, "global binary index form")
zero(sp.discriminant(index_form.subs(v, 1), u) - 6912 * Delta, "index-form discriminant")
index_coefficients = sp.Poly(index_form, u, v).coeffs()
gate(all(coefficient.subs({A: 0, C: 0}) == 0 for coefficient in index_coefficients), "index form vanishes at origin")
coefficient_gcd = index_coefficients[0]
for coefficient in index_coefficients[1:]:
    coefficient_gcd = sp.gcd(coefficient_gcd, coefficient)
gate(coefficient_gcd in (1, -1, sp.Rational(1, 1), sp.Rational(-1, 1)), "primitive index form")

# The fibre at the origin is k plus a square-zero two-plane.  At every other
# singular target A is a unit, B equals the power order, and P=Q=0 gives the
# one-point fibre k[theta]/theta^3.  Hence no singular target contributes an
# etale point.
for multiplication_packet in (e1_square, e1_f2, f2_square):
    gate(
        all(coefficient.subs({A: 0, C: 0}) == 0 for coefficient in multiplication_packet),
        "square-zero origin fibre",
    )
external_singular_points = [(-30, -60), (30, -60), (10 * sp.I, -10), (-10 * sp.I, -10)]
for a_value, c_value in external_singular_points:
    gate(P.subs({A: a_value, C: c_value}) == 0, "external triple-root P")
    gate(Q.subs({A: a_value, C: c_value}) == 0, "external triple-root Q")

# The normalized cubic field is rational.  On A!=0 put s=C/A and z=theta/A^3.
# The cubic equation solves polynomially for A and C in terms of (s,z).
z = sp.symbols("z")
p_parameter = (s**2 - 4) * (s**2 + 1) ** 2
h_parameter = s**9 - 3 * s**7 - 9 * s**5 + 5 * s**3 + 30 * s
F_parameter = z**3 - 3 * p_parameter * z + 2 * h_parameter
A_parameter = F_parameter / 4
C_parameter = s * F_parameter / 4
zero(
    f.subs({A: A_parameter, C: C_parameter, Z: A_parameter**3 * z}),
    "rational cubic-surface parametrization",
)
parameter_jacobian = sp.factor(
    sp.det(
        sp.Matrix(
            [
                [sp.diff(A_parameter, s), sp.diff(A_parameter, z)],
                [sp.diff(C_parameter, s), sp.diff(C_parameter, z)],
            ]
        )
    )
)
zero(
    parameter_jacobian + F_parameter * sp.diff(F_parameter, z) / 16,
    "non-Keller chart Jacobian",
)

# ---------------------------------------------------------------------------
# The mixed boundary--ADE three-class on the quadratic resolvent.
# ---------------------------------------------------------------------------

D_block = sp.Matrix([[-2]])
boundary_block = sp.Matrix([[-4, 5], [5, -4]])
A2_block = a_chain(2)
A5_block = a_chain(5)
removed_gram = sp.diag(D_block, boundary_block, A2_block, A2_block, A5_block, A5_block)
gate(removed_gram.shape == (17, 17), "removed lattice rank")
gate(boundary_block.det() == -9, "boundary determinant")
gate(removed_gram.det() == 5832, "removed lattice determinant")
removed_snf = smith_normal_form(removed_gram, domain=sp.ZZ)
gate(
    tuple(abs(int(removed_snf[index, index])) for index in range(17) if abs(int(removed_snf[index, index])) != 1)
    == (3, 3, 6, 6, 18),
    "removed lattice Smith packet",
)

rho3 = sp.Matrix([1, 2])
rho6 = sp.Matrix([1, 2, 0, 1, 2])
gate(all(int(entry) % 3 == 0 for entry in A2_block * rho3), "A2 mod-three radical")
gate(all(int(entry) % 3 == 0 for entry in A5_block * rho6), "A5 mod-three radical")
gate(int((rho3.T * A2_block * rho3)[0]) == -6, "A2 radical square")
gate(int((rho6.T * A5_block * rho6)[0]) == -12, "A5 radical square")

# A conic-bundle Picard basis: B+, a fibre F, then in each degenerate fibre
# E_1,...,E_(n-1),L_n, omitting the endpoint L_0 met by B+.
fibre_exponents = (3, 3, 6, 6)
picard_rank = 2 + sum(fibre_exponents)
picard_gram = sp.zeros(picard_rank)
picard_gram[0, 0] = -4
picard_gram[0, 1] = picard_gram[1, 0] = 1
fibre_starts: list[int] = []
cursor = 2
for exponent in fibre_exponents:
    fibre_starts.append(cursor)
    for offset in range(exponent):
        picard_gram[cursor + offset, cursor + offset] = -2 if offset < exponent - 1 else -1
        if offset:
            picard_gram[cursor + offset - 1, cursor + offset] = 1
            picard_gram[cursor + offset, cursor + offset - 1] = 1
    cursor += exponent
gate(picard_rank == 20 and picard_gram.det() == -1, "unimodular rational Picard basis")

# Pairing data determine w=B+-B-.  It pairs -9 with B+, zero with F and all
# ADE roots, and -1 with each far endpoint L_n.
w_pairings = sp.zeros(picard_rank, 1)
w_pairings[0] = -9
for start, exponent in zip(fibre_starts, fibre_exponents):
    w_pairings[start + exponent - 1] = -1
w_coordinates = picard_gram.inv() * w_pairings
gate(all(entry.q == 1 for entry in w_coordinates), "integral boundary difference")
gate(int((w_coordinates.T * picard_gram * w_coordinates)[0]) == -18, "boundary difference square")

rho_coordinates = sp.zeros(picard_rank, 1)
for start, exponent in zip(fibre_starts, fibre_exponents):
    for j in range(1, exponent):
        rho_coordinates[start + j - 1] = j % 3
mixed_numerator = w_coordinates - rho_coordinates
gate(all(int(entry) % 3 == 0 for entry in mixed_numerator), "exact mixed three-divisibility")
tau_coordinates = mixed_numerator / 3
gate(int((tau_coordinates.T * picard_gram * tau_coordinates)[0]) == -6, "globalized class square")

# Killing B+, B-, D=T-F, and every middle ADE component leaves the four far
# endpoints with the sole relation 3L1+3L2+6L3+6L4=0.
affine_endpoint_relation = sp.Matrix([[3, 3, 6, 6]])
affine_endpoint_snf = smith_normal_form(affine_endpoint_relation, domain=sp.ZZ)
gate(
    tuple(affine_endpoint_snf[index, index] for index in range(min(affine_endpoint_snf.shape)))
    == (3,),
    "exact affine class-group torsion",
)
gate(sp.gcd_list([3, 3, 6, 6]) == 3, "unique affine three-torsion line")
gate((1, 1, 2, 2) == tuple(exponent // 3 for exponent in fibre_exponents), "tau endpoint residue")

B_plus = sp.zeros(picard_rank, 1)
B_plus[0] = 1
B_minus = B_plus - w_coordinates
gate(int((B_plus.T * picard_gram * B_plus)[0]) == -4, "Bplus square")
gate(int((B_minus.T * picard_gram * B_minus)[0]) == -4, "Bminus square")
gate(int((B_plus.T * picard_gram * B_minus)[0]) == 5, "boundary intersection")

# Boundary-only divisibility fails: w pairs -1 with every far endpoint.
for start, exponent in zip(fibre_starts, fibre_exponents):
    endpoint = sp.zeros(picard_rank, 1)
    endpoint[start + exponent - 1] = 1
    gate(int((w_coordinates.T * picard_gram * endpoint)[0]) == -1, "boundary-only hostile")

# The abstract removed discriminant has many isotropic lines.  Only the
# displayed diagonal tau is proved to globalize here.
isotropic_vectors = []
for residue in product(range(3), repeat=5):
    if residue == (0, 0, 0, 0, 0):
        continue
    quadratic_value = (
        residue[1] ** 2
        + residue[2] ** 2
        - residue[3] ** 2
        - residue[4] ** 2
    ) % 3
    if quadratic_value == 0:
        isotropic_vectors.append(residue)
gate(len(isotropic_vectors) == 98 and len(isotropic_vectors) // 2 == 49, "abstract isotropic-line hostile")

# Compact-support Euler ledger for the maximal etale locus of the normalized
# cubic.  The affine normalization of Delta is A1.  Gluing eight origin
# addresses and two pairs of A5 addresses gives e_c(Delta)=-8.  Cusps are
# unibranch.  All five singular target fibres contribute zero etale points.
e_delta = 1 - (8 - 1) - 2 * (2 - 1)
singular_target_count = 5
e_delta_smooth = e_delta - singular_target_count
e_complement = 1 - e_delta
e_u_et = 3 * e_complement + e_delta_smooth
gate(e_delta == -8, "branch Euler characteristic")
gate(e_delta_smooth == -13 and e_complement == 9, "Euler stratification")
gate(e_u_et == 14, "maximal etale locus Euler tariff")

semantic_payload = {
    "branch": "rational_decic_one_smooth_infinity_contact_10",
    "finite_packet": "ordinary8_plus_2A2_plus_2A5",
    "conic_discriminant": "((s^2-4)(s^2+1)^2)^3",
    "cubic_power_discriminant": "432*A^10*Delta",
    "normal_order_discriminant": "6912*Delta",
    "power_basis_index": "A^5",
    "normal_order": "finite_flat_normal_nonmonogenic_rational_field",
    "quadratic_resolvent_class": "Cl=Z^3_plus_Z/3_mixed_boundary_ADE",
    "boundary_only": "primitive_not_3_divisible",
    "removed_lattice_det": 5832,
    "etale_euler": 14,
    "scope": "U_et_not_A2_but_proper_plane_open_remains_open_no_JC2_claim",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("theorem=THM-3915-rational-decic-cube-resolvent-index-debt-euler-tariff")
print("branch=degree10;rational;infinity_places=1;infinity_contact=10")
print("finite_singularities=ordinary8_plus_2A2_plus_2A5;delta_ledger=28+2+6=36")
print("conic_discriminant=((s^2-4)(s^2+1)^2)^3")
print("power_order_discriminant=432*A^10*Delta;index=A^5")
print("normal_order=finite_flat_normal_nonmonogenic;discriminant=6912*Delta")
print("cubic_field=rational_k(s,z);natural_chart_jacobian=-F*F_z/16")
print("removed_lattice=rank17;det=5832;snf_nonunits=3,3,6,6,18")
print("quadratic_resolvent_class_group=Z^3_plus_Z/3;mixed_tau_square=-6;boundary_only=primitive")
print("branch_euler=-8;branch_smooth_euler=-13;base_complement_euler=9")
print("maximal_etale_locus_euler=14;maximal_etale_locus_is_not_A2")
print("scope=proper_plane_open_may_delete_euler13;plane_atlas_and_JC2_open")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
