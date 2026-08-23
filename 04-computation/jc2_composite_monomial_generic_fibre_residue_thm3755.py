#!/usr/bin/env python3
"""Exact companion for THM-3755's composite-monomial residue theorem."""

from __future__ import annotations

import ast
import hashlib
from itertools import product
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


X, T, W, L = sp.symbols("X T W L")


def jac(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, X) * sp.diff(second, T)
        - sp.diff(first, T) * sp.diff(second, X)
    )


def total_degree_basis(bound: int) -> list[sp.Expr]:
    return [
        X**i * T**j
        for i in range(bound + 1)
        for j in range(bound + 1 - i)
    ]


# The monomial chart and the smoothness identity are uniform in m and F.
sample_profiles = (
    W,
    W**2,
    W + W**2,
    W - 2 * W**2 + W**3,
    2 * W + W**3,
)
smooth_controls = 0
factor_controls = 0
for m in range(2, 9):
    monomial = X**m * T
    gate(jac(X, monomial) == X**m, "monomial chart Jacobian")
    for profile in sample_profiles:
        component = sp.expand(X + profile.subs(W, monomial))
        qx = sp.diff(component, X)
        qt = sp.diff(component, T)
        gate(
            sp.expand(X * (qx - 1) - m * T * qt) == 0,
            "universal smoothness identity",
        )
        critical = sp.groebner((qx, qt), X, T, order="lex")
        gate(critical.contains(sp.Integer(1)), "smooth composite control")
        quotient = sp.cancel(profile / W)
        expected_factor = X * (
            1 + X ** (m - 1) * T * quotient.subs(W, monomial)
        )
        gate(sp.expand(component - expected_factor) == 0,
             "reducible zero-fibre factorization")
        smooth_controls += 1
        factor_controls += 1


# Formula for the generic-fibre residue.  If D_F=(1/F')d/dW and
# omega=dW/(L-F)^m, then at a simple root alpha of F-L,
# res_alpha(omega)=(-1)^m D_F^(m-1)(1/F')(alpha)/(m-1)!.
residue_controls = 0
for profile, target, root in ((W**2, 1, 1), (W + W**2, 2, 1), (W**3, 1, 1)):
    derivative = sp.diff(profile, W)
    for m in range(2, 8):
        invariant = 1 / derivative
        for _ in range(m - 1):
            invariant = sp.cancel(sp.diff(invariant, W) / derivative)
        predicted = sp.cancel(
            (-1) ** m * invariant.subs(W, root) / sp.factorial(m - 1)
        )
        actual = sp.residue(1 / (target - profile) ** m, W, root)
        gate(sp.cancel(actual - predicted) == 0,
             "generic-fibre residue formula")
        residue_controls += 1


# Hostile finite census for the integration lemma: no sampled nonlinear F
# makes D_F^(m-1)(1/F') vanish.  The theorem proves this at all degrees by
# integrating inside the differential field k(W), then using degrees.
invariant_controls = 0
for coefficients in product((-1, 0, 1), repeat=4):
    profile = sum(coefficients[index - 1] * W**index for index in range(1, 5))
    if sp.degree(profile, W) is sp.S.NegativeInfinity or sp.degree(profile, W) < 2:
        continue
    derivative = sp.diff(profile, W)
    for m in range(2, 7):
        invariant = 1 / derivative
        for _ in range(m - 1):
            invariant = sp.cancel(sp.diff(invariant, W) / derivative)
        gate(invariant != 0, "nonlinear residue invariant")
        invariant_controls += 1


# The linear profile has an exact rational primitive, but the primitive has
# negative X-order.  Verify its Jacobian and the polynomial positive boundary.
a, c = sp.symbols("a c", nonzero=True)
linear_boundary_controls = 0
for m in range(2, 10):
    component = X + a * X**m * T
    rational_mate = -c * X ** (1 - m) / (a * (m - 1))
    gate(sp.cancel(jac(rational_mate, component) - c) == 0,
         "linear-profile rational primitive")
    gate(jac(-c * T, X) == c, "coordinate polynomial mate")
    linear_boundary_controls += 1


# Bounded polynomial-mate systems are hostile controls for both the rationally
# obstructed nonlinear profiles and the rational-but-nonpolynomial boundary.
mate_profiles = (
    X + X**2 * T,
    X + (X**2 * T) ** 2,
    X + X**2 * T + (X**2 * T) ** 2,
    X + X**3 * T + (X**3 * T) ** 2,
    X + X**4 * T - 2 * (X**4 * T) ** 2 + (X**4 * T) ** 3,
)
mate_obstructions = 0
rank_gap_min = None
for profile_index, component in enumerate(mate_profiles):
    for bound in range(0, 11):
        basis = total_degree_basis(bound)
        coefficients = sp.symbols(
            f"p_{profile_index}_{bound}_0:{len(basis)}"
        )
        prospective = sum(
            coefficient * monomial
            for coefficient, monomial in zip(coefficients, basis)
        )
        equation = sp.Poly(jac(prospective, component) - 1, X, T)
        rows = [coefficient for _, coefficient in equation.terms()]
        matrix, rhs = sp.linear_eq_to_matrix(rows, coefficients)
        rank = matrix.rank()
        augmented_rank = matrix.row_join(rhs).rank()
        gap = augmented_rank - rank
        gate(gap == 1, "bounded mate augmented-rank obstruction")
        rank_gap_min = gap if rank_gap_min is None else min(rank_gap_min, gap)
        mate_obstructions += 1


# The advertised nonlinear multi-charge example has three distinct Euler
# charges for wt(X)=1, wt(T)=-1 and is smooth with a reducible zero fibre.
multi_charge = X + X**3 * T + X**6 * T**2
charges = tuple(i - j for (i, j), _ in sp.Poly(multi_charge, X, T).terms())
gate(charges == (4, 2, 1), "three Euler-charge sectors")
gate(
    sp.groebner(
        (sp.diff(multi_charge, X), sp.diff(multi_charge, T)),
        X,
        T,
        order="lex",
    ).contains(sp.Integer(1)),
    "multi-charge smoothness",
)
gate(
    sp.expand(multi_charge - X * (1 + X**2 * T + X**5 * T**2)) == 0,
    "multi-charge reducible fibre",
)


semantic_rows = (
    "family:Q=X+F(W);W=X^m*T;m>=2;F(0)=0",
    "smooth:Q_T=X^m*Fprime;X*(Q_X-1)=m*T*Q_T",
    "noncoordinate:F_nonzero_implies_reducible_zero_fibre",
    "fibre:dP/dW_at_L=-c*(L-F)^(-m)",
    "residue:(-1)^(m+1)c*D_F^(m-1)(1/Fprime)/(m-1)!",
    "nonlinear:no_rational_mate;linear:rational_not_polynomial;zero:coordinate",
    "control:X+X^3*T+X^6*T^2;charges=1,2,4",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3755-composite-monomial-generic-fibre-residue-obstruction")
print("scope=algebraically_closed_characteristic_zero;m>=2;F_in_k[W];F(0)=0")
print("geometry=all_Q_smooth;F_nonzero_reducible_zero_fibre_noncoordinate")
print("rational_mate=degF>=2_none;degF=1_explicit_negative_X_primitive")
print("polynomial_mate=iff_F_zero")
print("mechanism=generic_fibre_residue_and_D_F_integration_lemma")
print(
    f"smooth_controls={smooth_controls};factor_controls={factor_controls};"
    f"residue_controls={residue_controls};invariant_controls={invariant_controls}"
)
print(
    f"linear_boundary_controls={linear_boundary_controls};"
    f"mate_obstructions={mate_obstructions};rank_gap_min={rank_gap_min}"
)
print("multi_charge_control=X+X^3*T+X^6*T^2;Euler_charges=1,2,4")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
