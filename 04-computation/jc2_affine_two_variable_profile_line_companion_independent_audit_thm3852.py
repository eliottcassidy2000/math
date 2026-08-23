#!/usr/bin/env python3
"""Independent hostile audit of THM-3852.

This checker deliberately uses different elimination orders and proof
coordinates from the canonical companion.  It audits both charts of the
dual affine-line space, the normalized residual, every exceptional
parameter, and the vertical companion without relying on assertions (so
``python -O`` is an equally strong replay).
"""

from __future__ import annotations

import hashlib
import sys
from itertools import product

import sympy as sp

sys.stdout.reconfigure(newline="\n")


A, C, Z = sp.symbols("A C Z")
b0, bA, bC = sp.symbols("b0 bA bC")
d, e, m, n = sp.symbols("d e m n")
tau = sp.symbols("tau")
GATES = 0


def zero(label: str, value: sp.Expr) -> None:
    global GATES
    GATES += 1
    value = sp.factor(value)
    if value != 0:
        raise RuntimeError(f"{label}: expected zero, got {value}")


def nonzero(label: str, value: sp.Expr) -> None:
    global GATES
    GATES += 1
    value = sp.factor(value)
    if value == 0:
        raise RuntimeError(f"{label}: expected nonzero")


def equal(label: str, left: sp.Expr, right: sp.Expr) -> None:
    zero(label, left - right)


def same_ideal(
    label: str,
    left: list[sp.Expr],
    right: list[sp.Expr],
    *generators: sp.Symbol,
) -> None:
    """Check ideal equality by two independent Groebner reductions."""
    global GATES
    left_basis = sp.groebner(left, *generators, order="lex")
    right_basis = sp.groebner(right, *generators, order="lex")
    for item in left:
        remainder = sp.factor(right_basis.reduce(item)[1])
        if remainder != 0:
            raise RuntimeError(f"{label}: left remainder {remainder}")
    for item in right:
        remainder = sp.factor(left_basis.reduce(item)[1])
        if remainder != 0:
            raise RuntimeError(f"{label}: right remainder {remainder}")
    GATES += 1


def delta(profile: sp.Expr) -> sp.Expr:
    return sp.expand(
        -27 * A**2 * profile**2
        + 8 * A * C**3
        - 54 * A * C * profile
        + 9 * C**2
        - 54 * profile
    )


def homogenize(poly: sp.Expr, degree: int) -> sp.Expr:
    answer = 0
    for (a_degree, c_degree), coefficient in sp.Poly(poly, A, C).terms():
        answer += coefficient * A**a_degree * C**c_degree * Z ** (
            degree - a_degree - c_degree
        )
    return sp.expand(answer)


# -------------------------------------------------------------------------
# 1. Both affine-line charts: nonvertical A=mC+n and vertical C=c.
# -------------------------------------------------------------------------

profile = b0 + bA * A + bC * C
restricted_profile = e + d * C
nonvertical_restriction = sp.expand(
    delta(restricted_profile).subs(A, m * C + n)
)
nonvertical_coefficients = sp.Poly(nonvertical_restriction, C).all_coeffs()

# A lexicographic order different from the primary script collapses the five
# coefficient equations to this three-row triangular ideal.
alternative_triangular_basis = [
    d - e**2 * m * n - e * n,
    -18 * e * m**2 - 12 * m + n**2,
    (3 * e * m + 1) ** 2,
]
same_ideal(
    "nonvertical_coefficient_ideal_alternative_order",
    nonvertical_coefficients,
    alternative_triangular_basis,
    d,
    n,
    m,
    e,
)

# On k-points, the square in the third row forces em=-1/3.  The other two
# rows then force m=n^2/6, e=-2/n^2, and d=-4/(3n); in particular n is a
# unit.  Check the implication after clearing only powers known to be units.
e_solution = -2 / n**2
m_solution = n**2 / 6
d_solution = -sp.Rational(4, 3) / n
for index, coefficient in enumerate(nonvertical_coefficients):
    equal(
        f"nonvertical_solution_{index}",
        coefficient.subs({e: e_solution, m: m_solution, d: d_solution}),
        0,
    )
equal(
    "triangular_row_two_after_em",
    alternative_triangular_basis[1].subs(e, -1 / (3 * m)),
    n**2 - 6 * m,
)
equal(
    "triangular_row_one_after_solution",
    alternative_triangular_basis[0].subs(
        {e: e_solution, m: m_solution, d: d_solution}
    ),
    0,
)
nonzero("n_must_be_nonzero", 6)

lam = sp.symbols("lam")
line_n = A - n**2 * C / 6 - n
classified_profile = lam * line_n - 2 / n**2 - sp.Rational(4, 3) * C / n
equal("classified_line_restriction", classified_profile.subs(A, m_solution * C + n), e_solution + d_solution * C)
equal("classified_profile_A_coefficient", sp.diff(classified_profile, A), lam)
equal("classified_line_converse", delta(classified_profile).subs(A, m_solution * C + n), 0)

# The complementary chart is C=c.  Vanishing of the A^4, A^2, and A rows
# successively gives bA=0, b0+bC*c=0, c=0; hence b0=0.
c = sp.symbols("c")
vertical_polynomial = sp.Poly(delta(profile).subs(C, c), A)
equal("vertical_A4", vertical_polynomial.coeff_monomial(A**4), -27 * bA**2)
after_bA = sp.Poly(vertical_polynomial.as_expr().subs(bA, 0), A)
vertical_value = b0 + bC * c
equal("vertical_A2", after_bA.coeff_monomial(A**2), -27 * vertical_value**2)
equal(
    "vertical_A_after_value_zero",
    after_bA.coeff_monomial(A).subs(b0, -bC * c),
    8 * c**3,
)
equal("vertical_converse", delta(bC * C).subs(C, 0), 0)

# A small hostile census in the rational box catches chart, sign, and
# intercept mistakes independently of the symbolic ideal comparison.
nonvertical_hits = []
for b0_value, bA_value, bC_value, m_value, n_value in product(
    range(-2, 3), range(-2, 3), range(-2, 3), range(-3, 4), range(-3, 4)
):
    candidate = delta(b0_value + bA_value * A + bC_value * C).subs(
        A, m_value * C + n_value
    )
    if sp.Poly(candidate, C).is_zero:
        nonvertical_hits.append(
            (b0_value, bA_value, bC_value, m_value, n_value)
        )
        nonzero("finite_census_n_nonzero", n_value)
        equal("finite_census_m", sp.Rational(m_value), sp.Rational(n_value**2, 6))
        equal(
            "finite_census_b0",
            sp.Rational(b0_value),
            -sp.Rational(bA_value * n_value) - sp.Rational(2, n_value**2),
        )
        equal(
            "finite_census_bC",
            sp.Rational(bC_value),
            -sp.Rational(bA_value * n_value**2, 6)
            - sp.Rational(4, 3 * n_value),
        )
equal("finite_census_hit_count", sp.Integer(len(nonvertical_hits)), 0)

vertical_hits = []
for b0_value, bA_value, bC_value, c_value in product(
    range(-2, 3), range(-2, 3), range(-2, 3), range(-3, 4)
):
    candidate = delta(b0_value + bA_value * A + bC_value * C).subs(C, c_value)
    if sp.Poly(candidate, A).is_zero:
        vertical_hits.append((b0_value, bA_value, bC_value, c_value))
        equal("finite_vertical_b0", sp.Integer(b0_value), 0)
        equal("finite_vertical_bA", sp.Integer(bA_value), 0)
        equal("finite_vertical_c", sp.Integer(c_value), 0)
equal("finite_vertical_hit_count", sp.Integer(len(vertical_hits)), 5)


# -------------------------------------------------------------------------
# 2. Normalize the nonvertical chart and classify every residual factor.
# -------------------------------------------------------------------------

Ap, Cp = sp.symbols("Ap Cp")
scale = {A: n * Ap / 6, C: 6 * Cp / n, lam: 216 * tau / n**3}
equal("symplectic_product", scale[A] * scale[C], Ap * Cp)
equal("line_scaling", line_n.subs(scale), n * (Ap - 6 * Cp - 6) / 6)
normalized_profile = tau * (A - 6 * C - 6) - sp.Rational(1, 18) - sp.Rational(2, 9) * C
equal(
    "profile_scaling",
    n**2 * classified_profile.subs(scale) / 36,
    normalized_profile.subs({A: Ap, C: Cp}),
)
equal(
    "branch_covariance",
    delta(classified_profile).subs(scale),
    36 * delta(normalized_profile).subs({A: Ap, C: Cp}) / n**2,
)

line = A - 6 * C - 6
normalized_delta = sp.factor(delta(normalized_profile))
R = sp.factor(-12 * normalized_delta / line)
R_expected = (
    324 * tau**2 * A**3
    - 1944 * tau**2 * A**2 * C
    - 144 * tau * A**2 * C
    - 1944 * tau**2 * A**2
    - 36 * tau * A**2
    + 16 * A * C**2
    + 648 * tau * A * C
    + 8 * A * C
    + A
    + 18 * C
    + 648 * tau
    + 6
)
equal("normalized_factorization", normalized_delta, -line * R_expected / 12)
equal("residual_formula", R, R_expected)

# tau=0 is a primitive linear polynomial f(C)A+g(C).  An explicit Bezout
# identity makes f a unit in its quotient, proving that the ring is exactly
# k[C,1/(4C+1)] and not a hidden reducible or nonnormal model.
R_zero = sp.factor(R.subs(tau, 0))
f = (4 * C + 1) ** 2
g = 6 * (3 * C + 1)
equal("zero_seam", R_zero, f * A + g)
equal("zero_seam_gcd", sp.gcd(f, g), 1)
bezout_u = 9
bezout_v = -8 * C - sp.Rational(4, 3)
equal("zero_seam_bezout", bezout_u * f + bezout_v * g, 1)
equal(
    "zero_seam_inverse_in_quotient",
    sp.rem(f * (bezout_u - bezout_v * A) - 1, R_zero, A),
    0,
)

# For tau != 0, R has exact total degree three.  Any proper factorization
# must contain an affine line.  The vertical chart is impossible.  In the
# nonvertical chart, a second elimination order gives a triangular basis in
# N rather than the primary script's triangular basis in tau.
equal("residual_A3", sp.Poly(R, A).coeff_monomial(A**3), 324 * tau**2)
M, N = sp.symbols("M N")
line_factor_coefficients = sp.Poly(sp.expand(R.subs(A, M * C + N)), C).all_coeffs()
alternative_second_line_basis = [
    -2 * N**2 + 33 * N + 1944 * tau - 162,
    9 * M - N**2 - 6 * N + 18,
    (N - 6) ** 3,
]
same_ideal(
    "residual_line_ideal_alternative_order",
    line_factor_coefficients,
    alternative_second_line_basis,
    tau,
    M,
    N,
)
collision = {tau: sp.Rational(1, 54), M: 6, N: 6}
for index, coefficient in enumerate(line_factor_coefficients):
    equal(f"collision_line_{index}", coefficient.subs(collision), 0)
equal(
    "collision_tau_from_triangular_row",
    alternative_second_line_basis[0].subs({N: 6, tau: sp.Rational(1, 54)}),
    0,
)
equal(
    "collision_slope_from_triangular_row",
    alternative_second_line_basis[1].subs({N: 6, M: 6}),
    0,
)

# Homogenize independently, rather than copying only the top-degree terms.
# The support P0=[0:1:0] comes from A=0; q(0,1)=16 keeps every q-root away
# from P0.  Algebraic closure supplies a q-root even at the double-root row.
R_projective = homogenize(R, 3)
R_infinity = sp.factor(R_projective.subs(Z, 0))
q = sp.factor(R_infinity / A)
equal("infinity_full_factor", R_infinity, A * q)
equal("P0_on_residual", R_projective.subs({A: 0, C: 1, Z: 0}), 0)
equal("q_avoids_P0", q.subs({A: 0, C: 1}), 16)
equal("q_degree", sp.Poly(q, A, C).total_degree(), 2)
equal(
    "q_discriminant",
    sp.discriminant(q.subs(C, 1), A),
    139968 * tau**3 * (27 * tau + 4),
)
equal(
    "hostile_double_root",
    q.subs(tau, -sp.Rational(4, 27)),
    sp.Rational(16, 9) * (2 * A - 3 * C) ** 2,
)
nonzero("hostile_double_root_distinct_from_P0", (2 * A - 3 * C).subs({A: 0, C: 1}))

# At the unique collision the residual is L*K/9, so the full branch has a
# doubled L but its reduced companion is a genuinely irreducible smooth
# conic.  Its infinity packet is two distinct supports.
K = A**2 - 24 * A * C - 6 * A - 27
equal("collision_residual", R.subs(tau, sp.Rational(1, 54)), line * K / 9)
equal(
    "collision_full_branch",
    normalized_delta.subs(tau, sp.Rational(1, 54)),
    -line**2 * K / 108,
)
K_discriminant = sp.factor(sp.discriminant(K, A))
equal("collision_conic_discriminant", K_discriminant, 144 * (4 * C**2 + 2 * C + 1))
equal("collision_conic_squarefree", sp.gcd(4 * C**2 + 2 * C + 1, 8 * C + 2), 1)
nonzero("collision_conic_nonsquare", sp.discriminant(4 * C**2 + 2 * C + 1, C))
K_projective = homogenize(K, 2)
same_ideal(
    "collision_conic_smooth",
    [
        K_projective,
        sp.diff(K_projective, A),
        sp.diff(K_projective, C),
        sp.diff(K_projective, Z),
    ],
    [A, C, Z],
    A,
    C,
    Z,
)
equal("collision_conic_infinity", K_projective.subs(Z, 0), A * (A - 24 * C))
nonzero("collision_conic_infinity_distinct", 24)
equal(
    "collision_whole_infinity_packet",
    R_infinity.subs(tau, sp.Rational(1, 54)),
    sp.Rational(1, 9) * A * (A - 6 * C) * (A - 24 * C),
)


# -------------------------------------------------------------------------
# 3. Rebuild the vertical companion, including its endpoint, from scratch.
# -------------------------------------------------------------------------

kappa, V = sp.symbols("kappa V")
kinv = sp.symbols("kinv")
vertical_delta = sp.factor(delta(kappa * C))
H = sp.factor(vertical_delta / C)
equal("vertical_factorization", vertical_delta, C * H)
equal("vertical_residual_over_C0", H.subs(C, 0), -54 * kappa)
H_coefficients = sp.Poly(H, A).all_coeffs()
equal("vertical_residual_primitive", sp.gcd_list(H_coefficients), 1)
equal(
    "vertical_residual_discriminant",
    sp.factor(sp.discriminant(H, A)),
    -8 * C * (9 * kappa - 2 * C) ** 3,
)

# The normalization model V^2=C(9*kappa-2C) is a smooth conic when kappa is
# nonzero.  Its C=0 point is absent from H because H(A,0)=-54*kappa; its two
# projective infinity points are distinct.
vertical_conic = V**2 - C * (9 * kappa * Z - 2 * C)
same_ideal(
    "vertical_conic_smooth_after_kappa_unit",
    [
        vertical_conic,
        sp.diff(vertical_conic, V),
        sp.diff(vertical_conic, C),
        sp.diff(vertical_conic, Z),
        kappa * kinv - 1,
    ],
    [V, C, Z, kappa * kinv - 1],
    V,
    C,
    Z,
    kappa,
    kinv,
)
equal("vertical_conic_C0", vertical_conic.subs({C: 0, Z: 1}), V**2)
equal("vertical_conic_infinity", vertical_conic.subs(Z, 0), V**2 + 2 * C**2)
nonzero("vertical_conic_two_infinity_points", sp.discriminant(V**2 + 2, V))

# At kappa=0 the residual hyperbola is itself a smooth conic with the two
# supports [1:0:0] and [0:1:0] removed.
endpoint_delta = sp.factor(delta(0))
endpoint_residual = 8 * A * C + 9
equal("vertical_endpoint_packet", endpoint_delta, C**2 * endpoint_residual)
endpoint_projective = homogenize(endpoint_residual, 2)
same_ideal(
    "endpoint_hyperbola_projective_smooth",
    [
        endpoint_projective,
        sp.diff(endpoint_projective, A),
        sp.diff(endpoint_projective, C),
        sp.diff(endpoint_projective, Z),
    ],
    [A, C, Z],
    A,
    C,
    Z,
)
equal("endpoint_hyperbola_infinity", endpoint_projective.subs(Z, 0), 8 * A * C)


semantic_packet = (
    "THM-3852 independent hostile audit",
    "dual affine-line charts exhausted",
    "nonvertical ideal triangulated in d,n,m,e order",
    "residual line ideal triangulated in tau,M,N order",
    "tau zero localization has explicit Bezout inverse",
    "generic residual has P0 plus a disjoint quadratic infinity support",
    "tau=-4/27 retains two supports after q coalescence",
    "tau=1/54 is the unique doubled-line row with a smooth two-place conic",
    "vertical nonzero profiles have a three-place conic companion",
    "vertical zero endpoint has a two-place hyperbola companion",
)

print("THM3852_HOSTILE_AUDIT", "PASS")
print("LINE_CHARTS", "A=mC+n and C=c complete")
print("NONVERTICAL_SOLUTION", "n!=0,m=n^2/6,e=-2/n^2,d=-4/(3n)")
print("FINITE_NONVERTICAL_HITS", len(nonvertical_hits))
print("FINITE_VERTICAL_HITS", len(vertical_hits))
print("NORMALIZED_TRANSVERSE", "tau=n^3*lambda/216")
print("ZERO_SEAM_RING", "k[C,1/(4C+1)]")
print("RESIDUAL_REDUCIBILITY", "tau=1/54 only; repeated original line")
print("HOSTILE_INFINITY", "tau=-4/27: P0 plus q=(16/9)(2A-3C)^2")
print("COLLISION_COMPANION", "smooth irreducible conic; two infinity supports")
print("VERTICAL_COMPANION", "kappa!=0: three punctures; kappa=0: two")
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("GATES", GATES)
