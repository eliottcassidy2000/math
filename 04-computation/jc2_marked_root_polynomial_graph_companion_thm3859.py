#!/usr/bin/env python3
"""Assertion-free exact gates for the THM-3859 proof candidate.

Reproduction:
  python3 04-computation/jc2_marked_root_polynomial_graph_companion_thm3859.py
  python3 -O 04-computation/jc2_marked_root_polynomial_graph_companion_thm3859.py
"""

from __future__ import annotations

import hashlib

import sympy as sp


A, C, S, Q, V, Baxis = sp.symbols("A C S Q V Baxis")
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


def branch(profile: sp.Expr) -> sp.Expr:
    return sp.expand(
        -27 * A**2 * profile**2
        + 8 * A * C**3
        - 54 * A * C * profile
        + 9 * C**2
        - 54 * profile
    )


# Cusp normalization on a marked branch.  S and Q stand for arbitrary
# values s(A),q(A); no derivative or degree assumption enters the identities.
F = C - 6 * S * (1 + A * S)
b0 = 2 * S**2 * (3 + 4 * A * S)
b = b0 + Q * F
Delta = branch(b)
H = sp.cancel(Delta / F)
zero("companion_is_polynomial", sp.denom(H) - 1)
H = sp.expand(H)
equal("marked_root_factorization", Delta, F * H)

z = 1 + 2 * A * S
P = 1 + sp.Rational(2, 3) * A * C
U = 1 + A * C + A**2 * b
equal("cusp_square_on_graph", P.subs(C, 6 * S * (1 + A * S)), z**2)
equal("cusp_cube_on_graph", U.subs(C, 6 * S * (1 + A * S)), z**3)
equal("profile_on_graph", b.subs(C, 6 * S * (1 + A * S)), b0)

H_expected = (
    8 * A * C**2
    + C
    * (
        -27 * A**2 * Q**2
        + 48 * A**2 * S**2
        - 54 * A * Q
        + 48 * A * S
        + 9
    )
    + 162 * A**3 * Q**2 * S**2
    - 432 * A**3 * Q * S**3
    + 288 * A**3 * S**4
    + 162 * A**2 * Q**2 * S
    - 324 * A**2 * Q * S**2
    + 144 * A**2 * S**3
    + 18 * A * S**2
    - 54 * Q
    + 54 * S
)
equal("explicit_companion", H, H_expected)

# The residual quadratic has two canonical discriminant flanks.
Hpoly = sp.Poly(H, C)
equal("quadratic_leading_coefficient", Hpoly.coeff_monomial(C**2), 8 * A)
B = sp.factor(Hpoly.coeff_monomial(C))
E = sp.factor(Hpoly.coeff_monomial(1))
D = sp.factor(sp.discriminant(H, C))
Lplus = A * Q + 4 * A * S + 3
Lminus = 3 * A * Q - 4 * A * S + 1
equal("companion_discriminant", D, 27 * Lplus * Lminus**3)
equal("special_fibre_companion", H.subs(A, 0), 9 * C - 54 * Q + 54 * S)
equal("special_fibre_B", B.subs(A, 0), 9)
equal("special_fibre_discriminant", D.subs(A, 0), 81)
nonzero("two_unramified_special_places", D.subs({A: 0, S: 2, Q: 5}))

# Complete reducible case.  If D=W^2, choose the sign W(0)=9 and write
# W=9+A*V.  One linear factor is divisible by A; the other has constant 18.
B1 = sp.cancel((B - 9) / A)
zero("B_minus_9_divisible_by_A", sp.denom(B1) - 1)
W = 9 + A * V
small_factor = sp.expand(16 * A * C + B - W)
pole_factor = sp.expand(16 * A * C + B + W)
equal("small_factor_A_divisibility", small_factor, A * (16 * C + B1 - V))
equal("pole_factor_special_unit", pole_factor.subs(A, 0), 18)
equal(
    "quadratic_difference_of_squares",
    (16 * A * C + B) ** 2 - D,
    32 * A * H,
)
equal(
    "formal_square_root_factorization_error",
    small_factor * pole_factor - 32 * A * H,
    D - W**2,
)
nonzero("pole_component_missing_A0", pole_factor.subs(A, 0))

# Positive and hostile controls.
equal("zero_profile_F", F.subs({S: 0, Q: 0}), C)
equal("zero_profile_packet", Delta.subs({S: 0, Q: 0}), C**2 * (8 * A * C + 9))
equal(
    "nonconstant_graph_control",
    F.subs({S: A, Q: 0}),
    C - 6 * A * (1 + A**2),
)
equal(
    "nonconstant_graph_discriminant",
    D.subs({S: A, Q: 0}),
    27 * (4 * A**2 + 3) * (1 - 4 * A**2) ** 3,
)

# The unique THM-3852 line collision appears as a square-discriminant control.
collision = {S: -sp.Rational(1, 6), Q: -sp.Rational(1, 3)}
line = A - 6 * C - 6
conic = A**2 - 24 * A * C - 6 * A - 27
equal("collision_graph", F.subs(collision), -line / 6)
equal("collision_profile", b.subs(collision), (A - 18 * C - 9) / 54)
equal("collision_companion", H.subs(collision), line * conic / 18)
equal("collision_square_discriminant", D.subs(collision), (A - 3) ** 4)
equal("excluded_axis_boundary", branch(Baxis).subs(A, 0), 9 * C**2 - 54 * Baxis)

print("THM3859_MARKED_ROOT", "z=1+2As,C=6s(1+As),b|F=2s^2(3+4As)")
print("THM3859_FACTOR", "Delta_b=F*H, F=C-6s(1+As)")
print("THM3859_COMPANION_LC", "8A")
print("THM3859_COMPANION_DISC", "27(Aq+4As+3)(3Aq-4As+1)^3")
print("THM3859_SPECIAL_FIBRE", "H(0,C)=9C-54q(0)+54s(0), disc(0)=81")
print("THM3859_IRREDUCIBLE", "one A=0 pole plus at least one A=infinity place")
print("THM3859_REDUCIBLE", "one graph factor plus one G_m factor")
print("THM3859_SCOPE", "s,q in k[A]; C-dependent transverse quotient remains open")
semantic_packet = (
    "cusp-UFD marked-root law",
    "all polynomial graph branches C=6s(1+As)",
    "one-variable transverse quotient b=b0+qF",
    "exact residual quadratic and two-flank discriminant",
    "irreducible companion has finite pole plus infinity",
    "reducible companion contains canonical Gm factor",
    "THM3852 collision recovered exactly",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
