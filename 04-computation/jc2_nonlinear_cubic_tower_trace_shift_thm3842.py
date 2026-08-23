#!/usr/bin/env python3
"""Exact gates for THM-3842.

Reproduction:
  python3 04-computation/jc2_nonlinear_cubic_tower_trace_shift_thm3842.py
  python3 -O 04-computation/jc2_nonlinear_cubic_tower_trace_shift_thm3842.py
"""

from __future__ import annotations

import hashlib

import sympy as sp


A, C, T, tau, p, u, q = sp.symbols("A C T tau p u q")
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


# THM-3811's marked-root characteristic cubic and intrinsic discriminant.
f = T**3 - C * T**2 + 7 * A**2 * T + 3 * A**3 - A**2 * C**2
Delta = (
    A * (C + 5 * A) * (4 * C + 19 * A) * (3 * C - 17 * A)
    + C**2 * (162 * A**3 + 126 * A**2 * C - 4 * C**3)
    - 27 * A**2 * C**4
)

# The unique trace shift gives the m=2 THM-3774 trinomial.
p0 = (C**2 - 21 * A**2) / 6
u0 = (81 * A**3 - 27 * A**2 * C**2 + 63 * A**2 * C - 2 * C**3) / 54
shifted = sp.expand(f.subs(T, tau + C / 3))
equal("trace_shift", shifted, tau**3 - 2 * p0 * tau + 2 * u0)
equal("marked_trace", sp.expand(-sp.diff(f, T, 2).subs(T, 0) / 2), C)

# The square A^2 is exactly the index-form square of the marked root.
equal("index_square_discriminant", sp.discriminant(f, T), A**2 * Delta)
equal("tower_cusp_pullback", 4 * (8 * p0**3 - 27 * u0**2), A**2 * Delta)

J = sp.factor(sp.diff(p0, A) * sp.diff(u0, C) - sp.diff(p0, C) * sp.diff(u0, A))
J_expected = A * (42 * A**2 * C - 49 * A**2 - 9 * A * C + 2 * C**3) / 6
equal("coefficient_jacobian", J, J_expected)
nonzero("coefficient_jacobian_nonconstant", sp.diff(J, C))

# Over k(p), the coefficient field is the smooth conic C^2-21A^2=6p.
# On that conic u has two poles, each of order four.  The displayed identity
# also supplies a finite elimination check of total degree eight.
conic = C**2 - 21 * A**2 - 6 * p
u_conic_numerator = -567 * A**4 + 81 * A**3 - 162 * p * A**2 + (21 * A**2 - 12 * p) * C
equal(
    "conic_u_identity",
    54 * u0 - u_conic_numerator,
    (27 * A**2 + 2 * C) * (21 * A**2 - C**2 + 6 * p),
)
equal("pole_lead_both_infinities", sp.LC(sp.Poly(u_conic_numerator, A), A), -567)

F1 = 21 * A**2 - C**2 + 6 * p
F2 = -81 * A**3 + 27 * A**2 * C**2 - 63 * A**2 * C + 2 * C**3 + 54 * u
elim_A = sp.factor(sp.resultant(F1, F2, C) / -27)
equal("degree_eight_eliminant", sp.degree(elim_A, A), 8)
equal("degree_eight_lead", sp.LC(sp.Poly(elim_A, A), A), 11907)
nonzero("generic_finiteness", J)

# THM-3811's branch normalization and the tower-cusp normalization.
R = (q - 3) * (q + 1) * (q + 2)
Aq = -2 * q**2 * R / (3 * q**2 + 7) ** 2
Cq = -q * R / (3 * q**2 + 7)
equal("branch_parametrization", Delta.subs({A: Aq, C: Cq}), 0)

r = sp.factor(Aq * q - Cq / 3)
r_expected = -q * R * (3 * q**2 - 7) / (3 * (3 * q**2 + 7) ** 2)
equal("double_root_parameter", r, r_expected)
equal("cusp_p", p0.subs({A: Aq, C: Cq}), sp.Rational(3, 2) * r**2)
equal("cusp_u", u0.subs({A: Aq, C: Cq}), r**3)

r_num, r_den = map(sp.factor, sp.together(r).as_numer_denom())
equal("branch_map_degree", max(sp.degree(r_num, q), sp.degree(r_den, q)), 6)
equal("branch_map_coprime", sp.gcd(r_num, r_den), 1)
equal("zero_fibre_six_simple", sp.factor(r_num), -q * R * (3 * q**2 - 7))
equal("infinity_fibre_three_double", sp.factor(r_den), 3 * (3 * q**2 + 7) ** 2)

H7 = 9 * q**7 + 63 * q**5 + 27 * q**4 - 539 * q**3 - 378 * q**2 + 343 * q + 147
equal("branch_derivative", sp.diff(r, q), -2 * H7 / (3 * (3 * q**2 + 7) ** 3))
nonzero("critical_polynomial_squarefree", sp.resultant(H7, sp.diff(H7, q), q))
nonzero("critical_avoids_zero_fibre", sp.resultant(H7, q * R * (3 * q**2 - 7), q))
nonzero("critical_avoids_pole_fibre", sp.resultant(H7, 3 * q**2 + 7, q))
equal("riemann_hurwitz_total", 3 + sp.degree(H7, q), 2 * 6 - 2)

# The marked symplectic branch-residue passport.  Multiplying a discriminant
# representative by a square changes the residue by a square, so the odd
# part of its divisor on the branch normalization is invariant.
Delta_C = sp.diff(Delta, C)
eta_nonlinear = sp.factor(sp.diff(Aq, q) / Delta_C.subs({A: Aq, C: Cq}))
eta_nonlinear_expected = (3 * q**2 + 7) ** 4 / (
    2 * q**3 * R**3 * (3 * q**2 - 7) ** 2
)
equal("nonlinear_branch_residue", eta_nonlinear, eta_nonlinear_expected)

eta_tower_pullback = sp.factor(-sp.diff(r, q) / (18 * r**2))
J_branch = sp.factor(J.subs({A: Aq, C: Cq}))
equal("residue_pullback_multiplier", eta_tower_pullback, 4 * J_branch * eta_nonlinear / Aq**2)

# eta_tower=-dr/(18r^2) has even divisor.  Modulo a rational square the
# nonlinear coefficient is q*R; its four distinct roots give the odd packet.
square_witness = sp.factor(eta_nonlinear / (q * R))
square_root_witness = (3 * q**2 + 7) ** 2 / (
    sp.sqrt(2) * q**2 * R**2 * (3 * q**2 - 7)
)
equal("nonlinear_residue_squareclass", square_witness, square_root_witness**2)
equal("four_odd_residue_points", sp.degree(q * R, q), 4)
nonzero("four_odd_residue_points_distinct", sp.discriminant(q * R, q))

print("THM3842_TRACE_SHIFT", sp.factor(shifted))
print("THM3842_COEFFICIENT_P", sp.factor(p0))
print("THM3842_COEFFICIENT_U", sp.factor(u0))
print("THM3842_INDEX_SQUARE", "disc(f)=A^2*Delta=4*(8*p^3-27*u^2)")
print("THM3842_COEFFICIENT_JACOBIAN", J)
print("THM3842_CONIC_POLE_DIVISOR", "4*infinity_plus+4*infinity_minus")
print("THM3842_COEFFICIENT_FIELD_DEGREE", 8)
print("THM3842_BRANCH_PARAMETER", r)
print("THM3842_BRANCH_ZERO_PACKET", "q=0,3,-1,-2,+sqrt(7/3),-sqrt(7/3)")
print("THM3842_BRANCH_INFINITY_PACKET", "three points of index 2")
print("THM3842_BRANCH_OTHER_RAMIFICATION", "seven simple critical points")
print("THM3842_TOWER_RESIDUE_PARITY", "zero")
print("THM3842_NONLINEAR_RESIDUE_PARITY", sp.factor(q * R))
print("THM3842_SCOPE", "exact marked pullback; no planar Keller map and no unmarked birational classification")
semantic_packet = (
    "trace-shift m2 tower",
    "index square A^2",
    "coefficient field degree 8",
    "branch normalization degree 6",
    "zero fibre 4+1+1",
    "ramification 2+2+2 plus seven simple",
    "marked residue parity 0 versus q(q-3)(q+1)(q+2)",
    "no planar Keller or unmarked birational conclusion",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
