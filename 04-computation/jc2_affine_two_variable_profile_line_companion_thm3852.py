#!/usr/bin/env python3
"""Assertion-free exact gates for the THM-3852 proof candidate.

Reproduction:
  python3 04-computation/jc2_affine_two_variable_profile_line_companion_thm3852.py
  python3 -O 04-computation/jc2_affine_two_variable_profile_line_companion_thm3852.py
"""

from __future__ import annotations

import hashlib

import sympy as sp


A, C, Z = sp.symbols("A C Z")
b0, bA, bC, m, n = sp.symbols("b0 bA bC m n")
CHECKS = 0


def zero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def nonzero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value == 0:
        raise AssertionError(f"{label}: unexpectedly zero")


def equal(label: str, left: sp.Expr, right: sp.Expr) -> None:
    zero(label, left - right)


def same_ideal(label: str, left: list[sp.Expr], right: list[sp.Expr], *gens: sp.Symbol) -> None:
    global CHECKS
    G_left = sp.groebner(left, *gens, order="lex")
    G_right = sp.groebner(right, *gens, order="lex")
    for poly in left:
        if sp.factor(G_right.reduce(poly)[1]) != 0:
            raise AssertionError(f"{label}: left generator does not reduce")
    for poly in right:
        if sp.factor(G_left.reduce(poly)[1]) != 0:
            raise AssertionError(f"{label}: right generator does not reduce")
    CHECKS += 1


def branch(profile: sp.Expr) -> sp.Expr:
    return sp.expand(
        -27 * A**2 * profile**2
        + 8 * A * C**3
        - 54 * A * C * profile
        + 9 * C**2
        - 54 * profile
    )


# General affine two-variable profile and a nonvertical line A=mC+n.
b = b0 + bA * A + bC * C
Delta = branch(b)
e, d = sp.symbols("e d")
line_restriction = e + d * C
Delta_line = sp.expand(branch(line_restriction).subs(A, m * C + n))
line_coeffs = sp.Poly(Delta_line, C).all_coeffs()
equal(
    "affine_profile_line_restriction",
    b.subs(A, m * C + n),
    (b0 + bA * n) + (bA * m + bC) * C,
)

line_groebner = [
    27 * d**3 * n + 54 * d**2 + 16 * e,
    (3 * d * n + 4) ** 2,
    18 * d * m + 3 * d * n**2 + 8 * n,
    d * n**3 - 4 * m + 2 * n**2,
    (6 * m - n**2) ** 2,
]
same_ideal("nonvertical_line_coefficient_ideal", line_coeffs, line_groebner, e, d, m, n)

m_solution = n**2 / 6
d_solution = -sp.Rational(4, 3) / n
e_solution = -2 / n**2
for index, coefficient in enumerate(line_coeffs):
    equal(
        f"nonvertical_solution_{index}",
        coefficient.subs({m: m_solution, d: d_solution, e: e_solution}),
        0,
    )
equal("nonvertical_solution_dn", d_solution * n, -sp.Rational(4, 3))
nonzero("nonvertical_intercept_forced_by_3dn_plus_4", 4)

# A vertical line C=c is possible only in the inherited one-variable boundary
# b=bC*C, C=0.  The A^4, A^2, and A coefficients force this successively.
c = sp.symbols("c")
Delta_vertical = sp.Poly(Delta.subs(C, c), A)
equal("vertical_A4", Delta_vertical.coeff_monomial(A**4), -27 * bA**2)
vertical_after_bA = sp.Poly(Delta_vertical.as_expr().subs(bA, 0), A)
e_vertical = b0 + bC * c
equal("vertical_A2", vertical_after_bA.coeff_monomial(A**2), -27 * e_vertical**2)
equal(
    "vertical_A_after_profile_zero",
    vertical_after_bA.coeff_monomial(A).subs(b0, -bC * c),
    8 * c**3,
)
equal("vertical_boundary", branch(bC * C).subs(C, 0), 0)
equal("constant_zero_vertical_packet", branch(0), C**2 * (8 * A * C + 9))
equal(
    "constant_zero_residual_infinity",
    (8 * A * C + 9 * Z**2).subs(Z, 0),
    8 * A * C,
)

# The diagonal symplectic rescaling reduces every n!=0 family to n=6.
tau = sp.symbols("tau")
line = A - 6 * C - 6
b_normal = tau * line - sp.Rational(1, 18) - sp.Rational(2, 9) * C
Delta_normal = sp.factor(branch(b_normal))

# If t is the free coefficient b_A before normalization, then
# A=(n/6)A', C=(6/n)C' and b'=(n^2/36)b take t to
# tau=t*n^3/216.  This checks both the profile and the branch scaling.
Ap, Cp, tfree = sp.symbols("Ap Cp tfree")
b_family = (
    tfree * (A - n**2 * C / 6 - n)
    - 2 / n**2
    - sp.Rational(4, 3) * C / n
)
scale_subs = {A: n * Ap / 6, C: 6 * Cp / n, tfree: 216 * tau / n**3}
b_scaled = sp.factor(n**2 * b_family.subs(scale_subs) / 36)
equal(
    "symplectic_profile_normalization",
    b_scaled,
    tau * (Ap - 6 * Cp - 6) - sp.Rational(1, 18) - sp.Rational(2, 9) * Cp,
)
equal(
    "branch_covariance_under_normalization",
    branch(b_family).subs(scale_subs),
    sp.Rational(36, 1)
    / n**2
    * branch(b_normal).subs({A: Ap, C: Cp}),
)
R = sp.factor(-12 * Delta_normal / line)
R_expected = (
    324 * A**3 * tau**2
    - 1944 * A**2 * C * tau**2
    - 144 * A**2 * C * tau
    - 1944 * A**2 * tau**2
    - 36 * A**2 * tau
    + 16 * A * C**2
    + 648 * A * C * tau
    + 8 * A * C
    + A
    + 18 * C
    + 648 * tau
    + 6
)
equal("normalized_branch_factorization", R, R_expected)
equal("normalized_line_factor", Delta_normal.subs(A, 6 * C + 6), 0)

# t=0 is inherited from THM-3850 and has an irreducible rational residual
# with one finite denominator puncture and infinity.
R0 = sp.factor(R.subs(tau, 0))
equal("zero_transverse_residual", R0, A * (4 * C + 1) ** 2 + 6 * (3 * C + 1))
equal("zero_transverse_denominator_root", (4 * C + 1).subs(C, -sp.Rational(1, 4)), 0)
nonzero("zero_transverse_numerator_at_root", (3 * C + 1).subs(C, -sp.Rational(1, 4)))

# For tau!=0, the residual is cubic.  Its infinity form has A=0 and at least
# one distinct root of the quadratic packet because q(0,C)=16C^2.
R_infinity = sp.expand(
    324 * tau**2 * A**3
    - (1944 * tau**2 + 144 * tau) * A**2 * C
    + 16 * A * C**2
)
q_infinity = sp.factor(R_infinity / A)
equal("infinity_factorization", R_infinity, A * q_infinity)
equal("infinity_quadratic_avoids_A0", q_infinity.subs(A, 0), 16 * C**2)
equal("infinity_quadratic_degree", sp.degree(q_infinity, A), 2)
equal(
    "infinity_quadratic_discriminant",
    sp.discriminant(q_infinity, A),
    139968 * C**2 * tau**3 * (27 * tau + 4),
)
equal("double_infinity_hostile", sp.discriminant(q_infinity, A).subs(tau, -sp.Rational(4, 27)), 0)
nonzero("double_infinity_still_avoids_A0", q_infinity.subs({tau: -sp.Rational(4, 27), A: 0, C: 1}))

# A reducible cubic residual would have an affine line factor.  Substituting
# A=M C+N gives the exact coefficient ideal below; the only solution is the
# original line at tau=1/54.  Vertical lines are excluded by the A^3 term.
Mline, Nline = sp.symbols("Mline Nline")
second_coeffs = sp.Poly(sp.expand(R.subs(A, Mline * C + Nline)), C).all_coeffs()
second_groebner = [
    9 * Mline - 233280 * tau**2 + 12528 * tau - 206,
    9 * Nline - 93312 * tau**2 + 5400 * tau - 122,
    (54 * tau - 1) ** 3,
]
same_ideal(
    "second_line_coefficient_ideal",
    second_coeffs,
    second_groebner,
    Mline,
    Nline,
    tau,
)
equal("second_line_tau", (54 * tau - 1).subs(tau, sp.Rational(1, 54)), 0)
equal(
    "second_line_slope",
    second_groebner[0].subs({tau: sp.Rational(1, 54), Mline: 6}),
    0,
)
equal(
    "second_line_intercept",
    second_groebner[1].subs({tau: sp.Rational(1, 54), Nline: 6}),
    0,
)
equal("residual_vertical_A3", sp.Poly(R.subs(C, c), A).coeff_monomial(A**3), 324 * tau**2)

# At the collision parameter the reduced packet is one line plus an
# irreducible conic with two distinct infinity points.
tau_collision = sp.Rational(1, 54)
K = A**2 - 24 * A * C - 6 * A - 27
equal(
    "collision_factorization",
    Delta_normal.subs(tau, tau_collision),
    -line**2 * K / 108,
)
equal("collision_conic_discriminant", sp.discriminant(K, A), 144 * (4 * C**2 + 2 * C + 1))
nonzero("collision_conic_nonsquare", sp.discriminant(4 * C**2 + 2 * C + 1, C))
K_infinity = A**2 - 24 * A * C
equal("collision_conic_infinity", K_infinity, A * (A - 24 * C))
nonzero("collision_conic_two_distinct_infinities", 24)
K_projective = A**2 - 24 * A * C - 6 * A * Z - 27 * Z**2
same_ideal(
    "collision_conic_projective_smoothness",
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

print("THM3852_LINE_CLASSIFICATION", "m=n^2/6,d=-4/(3n),e=-2/n^2,n!=0")
print("THM3852_VERTICAL_BOUNDARY", "b=bC*C,C=0; inherited THM3850")
print("THM3852_NORMAL_PROFILE", b_normal)
print("THM3852_BRANCH", "Delta=(A-6C-6)*(-R_tau/12)")
print("THM3852_ZERO_TRANSVERSE", R0)
print("THM3852_INFINITY_FORM", sp.factor(R_infinity))
print("THM3852_SECOND_LINE", "tau=1/54,M=N=6 only")
print("THM3852_COLLISION_REDUCED", "A-6C-6 plus A^2-24AC-6A-27")
print("THM3852_SCOPE", "affine b(A,C) with an affine-line branch; no arbitrary-profile or JC(2) claim")
semantic_packet = (
    "affine two-variable profile",
    "complete vertical/nonvertical affine-line classification",
    "normalized one-parameter transverse family",
    "zero transverse inherits two-place b(C) wall",
    "nonzero residual cubic irreducible except line-square collision",
    "two distinct projective infinity supports force two normalization places",
    "collision reduced packet is A1 plus two-place conic",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
