#!/usr/bin/env python3
"""Independent hostile audit for THM-3847.

This checker deliberately starts from the polynomial map

    (C,Z) |-> (A(C,Z), C)

obtained by solving the displayed monic cubic for A.  It does not import the
primary companion or use its multiplication table.  All gates remain active
under ``python -O``.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


A, C, B, Z, T, X, Y, W, E, RHO, H, V = sp.symbols(
    "A C B Z T X Y W E RHO H V"
)
GATES = 0


def fail(label: str, detail: object) -> None:
    raise RuntimeError(f"{label}: {detail}")


def gate_zero(label: str, value: sp.Expr) -> None:
    global GATES
    GATES += 1
    value = sp.factor(sp.cancel(value))
    if value != 0:
        fail(label, value)


def gate_equal(label: str, left: sp.Expr, right: sp.Expr) -> None:
    gate_zero(label, left - right)


def gate_true(label: str, value: object) -> None:
    global GATES
    GATES += 1
    if value is not True and value != sp.S.true:
        fail(label, value)


def reduce_quadratic_extensions(value: sp.Expr) -> sp.Expr:
    """Reduce modulo W^2=L, E^2=-8, RHO^2=-2."""
    numerator, denominator = sp.together(value).as_numer_denom()
    domain = sp.QQ.frac_field(A, C, B, V)
    poly = sp.Poly(numerator, W, E, RHO, domain=domain)
    divisors = [
        sp.Poly(W**2 + 2 * C**2 - 9 * B, W, E, RHO, domain=domain),
        sp.Poly(E**2 + 8, W, E, RHO, domain=domain),
        sp.Poly(RHO**2 + 2, W, E, RHO, domain=domain),
    ]
    remainder = sp.reduced(poly, divisors)[1].as_expr()
    return sp.factor(remainder / denominator)


def reduce_rho(value: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.together(value).as_numer_denom()
    domain = sp.QQ.frac_field(A, C, B, V)
    remainder = sp.rem(
        sp.Poly(numerator, RHO, domain=domain),
        sp.Poly(RHO**2 + 2, RHO, domain=domain),
    ).as_expr()
    return sp.factor(remainder / denominator)


# -------------------------------------------------------------------------
# 1. Recover the surface as A^2 and the finite cubic map independently.
# -------------------------------------------------------------------------

G = (
    Z**3
    - 4 * C * Z**2
    + (4 * C**2 + 6 * B) * Z
    - 4 * A * B**2
    - 12 * B * C
)
A_map = (
    Z**3
    - 4 * C * Z**2
    + (4 * C**2 + 6 * B) * Z
    - 12 * B * C
) / (4 * B**2)

gate_equal("G_is_linear_in_A_with_unit_coefficient", sp.diff(G, A), -4 * B**2)
gate_equal("solve_G_for_A", G.subs(A, A_map), 0)
gate_equal("monic_cubic_degree", sp.Poly(G, Z).degree(), 3)
gate_equal("monic_cubic_lead", sp.Poly(G, Z).LC(), 1)

# For fixed B in k*, the quotient by G is literally k[C,Z].  In particular
# it is smooth and normal; G_A is a nowhere-zero scalar.  The following
# specialization is a separate algorithmic irreducibility control for the
# proof's C=0 route.
G_C0 = sp.Poly(
    G.subs(C, 0), Z, domain=sp.QQ.frac_field(A, B)
)
gate_true("C0_specialization_irreducible", G_C0.is_irreducible)
gate_equal("C0_specialization", G.subs(C, 0), Z**3 + 6 * B * Z - 4 * A * B**2)

# The beta=0 boundary is hostile and must fail both irreducibility and
# smoothness exactly as the theorem's scope says.
gate_equal("hostile_beta_zero_factorization", G.subs(B, 0), Z * (Z - 2 * C) ** 2)


# -------------------------------------------------------------------------
# 2. Depressed generator, discriminants, critical and companion sheets.
# -------------------------------------------------------------------------

p = sp.Rational(3, 2) + A * C
u = 1 + A * C + A**2 * B
F = T**3 - 2 * p * T + 2 * u
t_map = (2 * C * Z - Z**2 - 4 * B) / (2 * B)
Delta = (
    -27 * A**2 * B**2
    + 8 * A * C**3
    - 54 * A * C * B
    + 9 * C**2
    - 54 * B
)
gate_equal("hostile_beta_zero_branch", Delta.subs(B, 0), C**2 * (8 * A * C + 9))

gate_equal("depressed_relation_on_A2", F.subs({T: t_map, A: A_map}), 0)
gate_equal(
    "field_generator_inverse_relation",
    (t_map**2 + t_map - 2).subs(A, A_map),
    A_map * Z,
)
gate_equal("depressed_discriminant", sp.discriminant(F, T), 4 * A**2 * Delta)
gate_equal("power_discriminant", sp.discriminant(G, Z), 16 * B**2 * Delta)

different = sp.diff(G, Z)
jacobian = sp.diff(A_map, Z)
gate_equal("map_jacobian_is_different", jacobian, different.subs(A, A_map) / (4 * B**2))
gate_equal(
    "different_norm_is_branch",
    sp.resultant(G, different, Z),
    -16 * B**2 * Delta,
)

# Pulling the branch equation to the source displays the ramified square and
# the visible unramified companion without using the primary multiplication
# table.
companion = 3 * Z**2 - 8 * C * Z + 24 * B
gate_equal(
    "ramified_square_and_companion",
    Delta.subs(A, A_map),
    -(different**2) * companion / (16 * B**2),
)
gate_equal("different_quadratic_discriminant", sp.discriminant(different, Z), -8 * (-2 * C**2 + 9 * B))
gate_equal("companion_quadratic_discriminant", sp.discriminant(companion, Z), -32 * (-2 * C**2 + 9 * B))
gate_true("different_is_nonconstant_in_source", sp.diff(different, Z) != 0)


# -------------------------------------------------------------------------
# 3. Branch irreducibility and exact normalization/puncture packet.
# -------------------------------------------------------------------------

L = -2 * C**2 + 9 * B
branch_poly = sp.Poly(Delta, A)
gate_equal("branch_A_degree", branch_poly.degree(), 2)
gate_equal("branch_A_leading_coefficient", branch_poly.LC(), -27 * B**2)
gate_equal("branch_A_discriminant", sp.discriminant(Delta, A), -8 * L**3)
gate_equal("L_squarefree_certificate", sp.discriminant(L, C), 72 * B)

Bquad = 8 * C**3 - 54 * B * C
A_norm = (Bquad - E * L * W) / (54 * B**2)
gate_zero("normalization_map", reduce_quadratic_extensions(Delta.subs(A, A_norm)))

# The critical curve in the smooth source is the same conic normalization.
Z_crit = (8 * C - E * W) / 6
gate_zero(
    "critical_curve_parameter",
    reduce_quadratic_extensions(different.subs(Z, Z_crit)),
)
gate_zero(
    "critical_curve_maps_to_normalization_formula",
    reduce_quadratic_extensions(A_map.subs(Z, Z_crit) - A_norm),
)

# On the ramification normalization, the repeated depressed-cubic root is a
# Laurent fold through the single standard cusp.  With E=2*RHO and
# U=W+RHO*C it is invariant under U |-> -U; both affine cusp addresses W=0
# map to repeated root zero.  This is the same mechanism exposed by the
# newly incoming THM-3844 one-place design, despite the different passport.
root_on_ramification = t_map.subs(Z, Z_crit)
gate_zero(
    "repeated_root_formula",
    reduce_quadratic_extensions(
        root_on_ramification - W * (C * E - 2 * W) / (18 * B)
    ),
)
gate_zero(
    "depressed_cusp_p_on_ramification",
    reduce_quadratic_extensions(
        p.subs(A, A_norm) - sp.Rational(3, 2) * root_on_ramification**2
    ),
)
gate_zero(
    "depressed_cusp_u_on_ramification",
    reduce_quadratic_extensions(
        u.subs(A, A_norm) - root_on_ramification**3
    ),
)
gate_zero(
    "Gm_to_depressed_cusp_fold",
    reduce_quadratic_extensions(
        root_on_ramification.subs(E, 2 * RHO)
        + sp.Rational(1, 2)
        + 9 * B / (2 * (W + RHO * C) ** 2)
    ),
)

# The only singular branch images lie over W=0 (equivalently L=0).  There
# are two such normalization points, and the map is pointwise injective:
# away from L=0, A recovers W; at L=0 the conic forces W=0.  Thus these are
# cuspidal failures, not the finite self-identifications forced in THM-3843.
gate_zero(
    "branch_gradient_A_on_normalization",
    reduce_quadratic_extensions(sp.diff(Delta, A).subs(A, A_norm) - E * W * L),
)
gate_zero(
    "branch_gradient_C_has_L_factor",
    reduce_quadratic_extensions(
        sp.diff(Delta, C).subs(A, A_norm)
        - L
        * (72 * B * C + 9 * B * E * W - 16 * C**3 - 4 * C**2 * E * W)
        / (9 * B**2)
    ),
)

# Smoothness has the explicit Jacobian-ideal certificate
#   conic - W/2*d_W(conic) - C/2*d_C(conic) = -9B,
# a unit for constant B != 0.
conic = W**2 + 2 * C**2 - 9 * B
smooth_certificate = conic - W * sp.diff(conic, W) / 2 - C * sp.diff(conic, C) / 2
gate_equal("affine_conic_smooth_certificate", smooth_certificate, -9 * B)

projective_conic = W**2 + 2 * C**2 - 9 * B * H**2
gate_zero(
    "infinity_factorization_extension",
    reduce_quadratic_extensions(
        projective_conic.subs(H, 0) - (W - RHO * C) * (W + RHO * C)
    ),
)

# An explicit Laurent coordinate proves that the affine normalization is
# G_m and shows exactly which two points were deleted.
Uplus = W + RHO * C
Uminus = W - RHO * C
gate_zero(
    "Gm_unit_product",
    reduce_quadratic_extensions(Uplus * Uminus - 9 * B),
)
W_from_U = (Uplus + 9 * B / Uplus) / 2
C_from_U = (Uplus - 9 * B / Uplus) / (2 * RHO)
gate_zero("Gm_inverse_W", reduce_quadratic_extensions(W_from_U - W))
gate_zero("Gm_inverse_C", reduce_quadratic_extensions(C_from_U - C))


# -------------------------------------------------------------------------
# 4. Binary index form, derived in the power basis rather than from (10).
# -------------------------------------------------------------------------

def reduce_mod_G(value: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.together(value).as_numer_denom()
    remainder = sp.rem(numerator, G, Z)
    return sp.cancel(remainder / denominator)


def power_vector(value: sp.Expr) -> sp.Matrix:
    reduced = sp.Poly(reduce_mod_G(value), Z)
    return sp.Matrix([
        reduced.coeff_monomial(1),
        reduced.coeff_monomial(Z),
        reduced.coeff_monomial(Z**2),
    ])


basis_change = sp.Matrix.hstack(
    power_vector(1), power_vector(t_map), power_vector(Z)
)
gate_equal("basis_change_determinant", basis_change.det(), 1 / (2 * B))

theta = X * t_map + Y * Z
power_index_det = sp.Matrix.hstack(
    power_vector(1), power_vector(theta), power_vector(theta**2)
).det()
index_form = sp.factor(power_index_det / basis_change.det())
index_expected = A * X**3 + 3 * X**2 * Y - 2 * C * X * Y**2 + 2 * B * Y**3
gate_equal("independent_binary_index", index_form, index_expected)
gate_equal("z_index_is_scalar_unit", index_form.subs({X: 0, Y: 1}), 2 * B)

# Discriminant transformation through the independently computed basis matrix.
gate_equal(
    "trace_basis_discriminant",
    basis_change.det() ** 2 * sp.discriminant(G, Z),
    4 * Delta,
)


# -------------------------------------------------------------------------
# 5. Hostile profile and scope gates.
# -------------------------------------------------------------------------

b = C * (2 * C + 1) / 9
Delta_b = sp.factor(Delta.subs(B, b))
gate_equal("profile_linearizes_radicand", L.subs(B, b), C)
gate_equal(
    "profile_branch_factorization",
    Delta_b,
    -C
    * (
        4 * A**2 * C**3
        + 4 * A**2 * C**2
        + A**2 * C
        + 12 * A * C**2
        + 18 * A * C
        + 9 * C
        + 18
    )
    / 3,
)

# The naive monogenic order for the nonconstant profile is singular all along
# C=Z=0: this is why its normalization, not the same power-basis order, is the
# only legitimate place to search for a new completion.
G_b = sp.factor(G.subs(B, b))
profile_singular_subs = {C: 0, Z: 0}
gate_equal("profile_singular_G", G_b.subs(profile_singular_subs), 0)
gate_equal("profile_singular_dA", sp.diff(G_b, A).subs(profile_singular_subs), 0)
gate_equal("profile_singular_dC", sp.diff(G_b, C).subs(profile_singular_subs), 0)
gate_equal("profile_singular_dZ", sp.diff(G_b, Z).subs(profile_singular_subs), 0)

# Sharpen that hostile: the residual component in (25) still has the exact
# three-puncture passport prohibited by THM-3841.  Its A-discriminant is
# -72*C.  With RHO^2=-2 and C=V^2, its normalization map simplifies to the
# displayed rational A(V); the apparent pole at V=RHO/2 cancels, while the
# poles V=0 and V=-RHO/2 and parameter infinity are genuinely missing.
Q_profile = (
    4 * A**2 * C**3
    + 4 * A**2 * C**2
    + A**2 * C
    + 12 * A * C**2
    + 18 * A * C
    + 9 * C
    + 18
)
gate_equal("profile_residual_A_discriminant", sp.discriminant(Q_profile, A), -72 * C)
A_profile_uncancelled = (
    -6 * V**2 * (2 * V**2 + 3) + 6 * RHO * V
) / (2 * V**2 * (2 * V**2 + 1) ** 2)
A_profile = -3 * (V + RHO) / (2 * V * (V + RHO / 2) ** 2)
gate_zero("profile_normalization_cancellation", reduce_rho(A_profile_uncancelled - A_profile))
gate_zero(
    "profile_residual_normalization_map",
    reduce_rho(Q_profile.subs({A: A_profile, C: V**2})),
)
gate_zero(
    "profile_finite_cancelled_point",
    reduce_rho(A_profile.subs(V, RHO / 2) - sp.Rational(9, 4)),
)
gate_true(
    "profile_pole_at_zero",
    reduce_rho(sp.cancel(V * A_profile).subs(V, 0)) != 0,
)
gate_true(
    "profile_pole_at_minus_rho_over_two",
    reduce_rho(
        sp.cancel((V + RHO / 2) ** 2 * A_profile).subs(V, -RHO / 2)
    )
    != 0,
)


# Pin the primary artifact only after all independent derivations have passed.
repo = Path(__file__).resolve().parents[1]
primary_script = repo / "04-computation" / "jc2_two_place_cubic_monogenic_unit_debt_thm3847.py"
primary_output = repo / "05-knowledge" / "results" / "jc2_two_place_cubic_monogenic_unit_debt_thm3847.out"
script_hash = hashlib.sha256(primary_script.read_bytes()).hexdigest()
output_hash = hashlib.sha256(primary_output.read_bytes()).hexdigest()
gate_true(
    "primary_script_hash",
    script_hash == "158de2197bc9e0c9d68923864b42792549cccb1f4357129f72c2567d46c5ae02",
)
gate_true(
    "primary_output_hash",
    output_hash == "1477e231edf1a8e5bd60e1b647b98adada8d27252ff7fe8afaefe6d4778d5f3e",
)


semantic_packet = (
    "completion is A2_CZ because G is linear in A with scalar coefficient",
    "finite flat cubic map has irreducible Gm ramification normalization",
    "branch pullback is different squared times one companion",
    "binary index represents one and the different is a nonconstant etale unit",
    "two punctures remain a Jelonek obstruction; the constant profile fails twice",
    "the advertised cheap nonconstant profile has a three-puncture residual branch",
    "next cubic move must normalize a nonconstant degeneration and become nonmonogenic",
)

print("THM3847_AUDIT", "PASS")
print("COMPLETION", "Spec(S_beta)=A2_(C,Z), hence smooth and normal")
print("IRREDUCIBILITY", "linear-in-A unit certificate plus irreducible C=0 control")
print("FINITE_FLAT_RANK", 3)
print("MAP_A", sp.factor(A_map))
print("MAP_JACOBIAN", sp.factor(jacobian))
print("BRANCH", Delta)
print("BRANCH_NORMALIZATION", "W^2+2*C^2=9*beta = G_m")
print("BRANCH_INFINITY_PLACES", 2)
print("RAMIFICATION", different)
print("COMPANION", companion)
print("BRANCH_PULLBACK", "Delta=-(different^2*companion)/(16*beta^2)")
print("INDEX", index_form)
print("INDEX_REPRESENTS_ONE", "(X,Y)=(0,(2*beta)^(-1/3))")
print("ETALE_UNIT", different)
print("HOSTILE_BETA_ZERO", sp.factor(G.subs(B, 0)))
print("HOSTILE_PROFILE_SINGULAR_LINE", "C=Z=0")
print("HOSTILE_PROFILE_RESIDUAL_NORMALIZATION", "P1 minus {0,-rho/2,infinity}")
print("THM3841_CONNECTION", "G_m has no dominant polynomial A1 curve, so the branch passport is still nonplanar")
print("THM3843_CONNECTION", "normalization map preserves C and is injective on points; it cannot be the forced self-identifying Russell arm")
print("THM3844_CONNECTION", "the ramification root is -1/2-9*beta/(2*U^2): two cusp addresses fold to one depressed cusp")
print("THM3846_CONNECTION", "next move is a higher-normal algebraizing seed, not another constant monogenic resummation")
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("PRIMARY_SCRIPT_SHA256", script_hash)
print("PRIMARY_OUTPUT_SHA256", output_hash)
print("GATES", GATES)
