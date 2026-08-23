#!/usr/bin/env python3
"""Exact gates for the THM-3853 inverse-discriminant obstruction.

Reproduction:
  python3 04-computation/jc2_inverse_discriminant_quadratic_depth_thm3853.py
  python3 -O 04-computation/jc2_inverse_discriminant_quadratic_depth_thm3853.py

The two Groebner computations are over QQ.  Each saturates the exact
coefficient ideal by the requested nonzero fifth-degree coefficient.
"""

from __future__ import annotations

import hashlib

import sympy as sp


A, C, X, Y, t, lam, eta = sp.symbols("A C X Y t lambda eta")
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


def equal(label: str, left: object, right: object) -> None:
    global CHECKS
    CHECKS += 1
    if left != right:
        raise AssertionError(f"{label}: {left!r} != {right!r}")


def binary_discriminant(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    return sp.expand(
        b**2 * c**2
        - 4 * a * c**3
        - 4 * b**3 * d
        - 27 * a**2 * d**2
        + 18 * a * b * c * d
    )


# The rational homogeneous-linear packet of THM-3808.
a0, b0, c0, d0 = A, C, 7 * A, -3 * A
Delta0 = sp.factor(A * (C + 5 * A) * (4 * C + 19 * A) * (3 * C - 17 * A))
zero("base_discriminant", binary_discriminant(a0, b0, c0, d0) - Delta0)

# Two genuine one-place orientations.  In coordinates (M,L), write
# Delta0=L^4 D(M/L); Delta0+lambda L^5 is normalized by
# L=-D(t)/lambda, M=tL.
D_C = sp.factor(Delta0.subs({A: t, C: 1}))
C_C = -D_C / lam
A_C = t * C_C
delta_C = sp.expand(Delta0 + lam * C**5)
zero("C_orientation_parametrization", delta_C.subs({A: A_C, C: C_C}))
equal("C_orientation_degree", sp.degree(D_C, t), 4)
nonzero("C_orientation_four_distinct_origins", sp.discriminant(D_C, t))
nonzero("C_orientation_not_base_factor", Delta0.subs(C, 1))

L_sum = A + C
D_sum = sp.factor(Delta0.subs({A: t, C: 1 - t}))
L_param = -D_sum / lam
A_sum = t * L_param
C_sum = (1 - t) * L_param
delta_sum = sp.expand(Delta0 + lam * L_sum**5)
zero("sum_orientation_parametrization", delta_sum.subs({A: A_sum, C: C_sum}))
equal("sum_orientation_degree", sp.degree(D_sum, t), 4)
nonzero("sum_orientation_four_distinct_origins", sp.discriminant(D_sum, t))
nonzero("sum_orientation_not_base_factor", Delta0.subs({A: t, C: 1 - t}))

# General homogeneous quadratic perturbation of all four coefficients.
u = sp.symbols("u0:12")
quadratic_monomials = (A**2, A * C, C**2)
qa, qb, qc, qd = (
    sum(u[3 * row + j] * quadratic_monomials[j] for j in range(3))
    for row in range(4)
)
a, b, c, d = A + qa, C + qb, 7 * A + qc, -3 * A + qd
Delta = binary_discriminant(a, b, c, d)
error = sp.Poly(sp.expand(Delta - Delta0), A, C)

# The linear part is frozen, so degree four is exactly Delta0.  Every index
# coefficient remains in the maximal ideal (A,C).
for degree in range(0, 5):
    for i in range(degree + 1):
        zero(
            f"no_error_degree_{degree}_{i}",
            error.coeff_monomial(A**i * C ** (degree - i)),
        )
for label, coefficient in (("a", a), ("b", b), ("c", c), ("d", d)):
    zero(f"index_coefficient_{label}_at_origin", coefficient.subs({A: 0, C: 0}))


def target_equations(line: sp.Expr) -> tuple[list[sp.Expr], sp.Expr]:
    """Equations for Delta-Delta0=lambda*line^5, with lambda eliminated."""

    target = sp.Poly(sp.expand(line**5), A, C)
    # Both chosen lines have nonzero C^5 coefficient, so the C^5 row fixes
    # the target scalar without adding another variable.
    target_C5 = target.coeff_monomial(C**5)
    nonzero("target_C5_coefficient", target_C5)
    scalar = sp.cancel(error.coeff_monomial(C**5) / target_C5)
    equations: list[sp.Expr] = []
    for degree in range(5, 9):
        for i in range(degree + 1):
            monomial = A**i * C ** (degree - i)
            wanted = scalar * target.coeff_monomial(monomial) if degree == 5 else 0
            equation = sp.expand(error.coeff_monomial(monomial) - wanted)
            if equation != 0:
                equations.append(equation)
    return equations, sp.factor(scalar)


eq_C, scalar_C = target_equations(C)
equal("C_orientation_equation_count", len(eq_C), 29)
zero("C_orientation_scalar_formula", scalar_C + 4 * u[11])

# Lambda=0 is a necessary hostile control: the unperturbed packet satisfies
# all unsaturated coefficient equations.  The obvious q_d=C^2 perturbation
# supplies the desired -4 C^5 row but leaves active mixed/higher debt.
zero_substitution = {variable: 0 for variable in u}
for number, equation in enumerate(eq_C):
    zero(f"C_zero_control_{number}", equation.subs(zero_substitution))
qd_C2 = {variable: 0 for variable in u}
qd_C2[u[11]] = 1
equal("C_active_control_scalar", scalar_C.subs(qd_C2), -4)
nonzero(
    "C_active_control_has_residual_debt",
    sum(value**2 for value in (equation.subs(qd_C2) for equation in eq_C)),
)

G_C = sp.groebner(eq_C + [eta * scalar_C - 1], eta, *u, order="grevlex")
equal("C_orientation_saturated_basis", list(G_C), [sp.Integer(1)])

eq_sum, scalar_sum = target_equations(A + C)
equal("sum_orientation_equation_count", len(eq_sum), 29)
zero("sum_orientation_scalar_formula", scalar_sum + 4 * u[11])
for number, equation in enumerate(eq_sum):
    zero(f"sum_zero_control_{number}", equation.subs(zero_substitution))

G_sum = sp.groebner(
    eq_sum + [eta * scalar_sum - 1], eta, *u, order="grevlex"
)
equal("sum_orientation_saturated_basis", list(G_sum), [sp.Integer(1)])

print("THM3853_BASE_PACKET", "(A,C,7A,-3A)")
print("THM3853_BASE_DISCRIMINANT", Delta0)
print("THM3853_C_NORMALIZATION_D", D_C)
print("THM3853_SUM_NORMALIZATION_D", D_sum)
print("THM3853_C_TARGET_SCALAR", scalar_C)
print("THM3853_SUM_TARGET_SCALAR", scalar_sum)
print("THM3853_C_EQUATIONS", len(eq_C))
print("THM3853_SUM_EQUATIONS", len(eq_sum))
print("THM3853_C_SATURATED_GROEBNER", list(G_C))
print("THM3853_SUM_SATURATED_GROEBNER", list(G_sum))
print("THM3853_ZERO_CONTROL", "lambda=0 survives at q_a=q_b=q_c=q_d=0")
print(
    "THM3853_ACTIVE_CONTROL",
    "q_d=C^2 gives lambda=-4 but leaves mixed degree-five and degree-six debt",
)
print(
    "THM3853_SCOPE",
    "no homogeneous-quadratic perturbation realizes either named one-place target",
)
semantic_packet = (
    "base binary cubic (A,C,7A,-3A)",
    "two one-place targets Delta0+lambda*C^5 and Delta0+lambda*(A+C)^5",
    "normalization t=M/L, L=-D(t)/lambda, M=tL",
    "four distinct finite normalization addresses glue to the origin",
    "all quadratic-perturbed index coefficients remain in (A,C)",
    "29 exact coefficient equations per orientation",
    "saturation by nonzero target coefficient has Groebner basis [1]",
    "lambda zero and q_d=C^2 hostile controls",
    "first unresolved algebraic coefficient depth is at least three",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
