#!/usr/bin/env python3
"""Exact companion for THM-3741's full radial two-charge classification."""

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


X, T, Z = sp.symbols("X T Z")


def jac(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, X) * sp.diff(second, T)
        - sp.diff(first, T) * sp.diff(second, X)
    )


def is_zero(expression: sp.Expr) -> bool:
    return sp.expand(expression) == 0


def is_constant(expression: sp.Expr) -> bool:
    return not sp.expand(expression).has(Z)


def is_nonzero_monomial(expression: sp.Expr) -> tuple[bool, int]:
    polynomial = sp.Poly(sp.expand(expression), Z)
    terms = polynomial.terms()
    if len(terms) != 1 or terms[0][1] == 0:
        return False, -1
    return True, terms[0][0][0]


def is_squarefree_with_nonzero_origin(expression: sp.Expr) -> bool:
    polynomial = sp.Poly(sp.expand(expression), Z)
    if polynomial.eval(0) == 0:
        return False
    return sp.gcd(polynomial, sp.diff(polynomial, Z)).degree() == 0


def predicted_nonsingular(phi: sp.Expr, psi: sp.Expr) -> bool:
    """The exact alternatives in THM-3741."""
    phi = sp.expand(phi)
    psi = sp.expand(psi)
    if is_constant(phi) and is_constant(psi):
        return not (is_zero(phi) and is_zero(psi))
    if is_zero(psi):
        return not is_zero(phi) and is_squarefree_with_nonzero_origin(phi)
    if is_zero(phi):
        return not is_zero(psi) and is_squarefree_with_nonzero_origin(psi)
    psi_monomial, psi_order = is_nonzero_monomial(psi)
    phi_monomial, phi_order = is_nonzero_monomial(phi)
    if is_constant(phi) and not is_zero(phi) and psi_monomial:
        return psi_order >= 2
    if is_constant(psi) and not is_zero(psi) and phi_monomial:
        return phi_order >= 2
    return False


# Universal torus critical determinant, now allowing arbitrary vanishing at
# z=0.  This is the common engine behind THM-3738 and the axis-boundary split.
phi_function = sp.Function("phi")(Z)
psi_function = sp.Function("psi")(Z)
A_profile = sp.diff(Z * phi_function, Z)
B_profile = sp.diff(Z * psi_function, Z)
critical_matrix = sp.Matrix((
    (A_profile, Z * sp.diff(psi_function, Z)),
    (Z * sp.diff(phi_function, Z), B_profile),
))
gate(sp.simplify(
    critical_matrix.det() - sp.diff(Z * phi_function * psi_function, Z)
) == 0, "universal radial critical determinant")

# The single-charge mate equations are literal one-variable derivatives.
f_coefficients = sp.symbols("f_0:7")
F = sum(f_coefficients[i] * Z**i for i in range(len(f_coefficients)))
profile_coefficients = sp.symbols("profile_0:6")
profile = sum(
    profile_coefficients[i] * Z**i for i in range(len(profile_coefficients))
)
lower_Q = sp.expand(X * profile.subs(Z, X * T))
lower_P = sp.expand(T * F.subs(Z, X * T))
upper_Q = sp.expand(T * profile.subs(Z, X * T))
upper_P = sp.expand(X * F.subs(Z, X * T))
gate(sp.expand(
    jac(lower_P, lower_Q)
    + sp.diff(Z * profile * F, Z).subs(Z, X * T)
) == 0, "single lower charge antiderivative")
gate(sp.expand(
    jac(upper_P, upper_Q)
    - sp.diff(Z * profile * F, Z).subs(Z, X * T)
) == 0, "single upper charge antiderivative")

# Exhaust the complete degree-at-most-two coefficient cube over {-1,0,1}.
# The affine critical ideal agrees exactly with the theorem's alternatives,
# including zero profiles, repeated roots, and both axis-vanishing directions.
grid_total = 0
grid_nonsingular = 0
grid_predicted = 0
values = (-1, 0, 1)
for phi_coefficients in product(values, repeat=3):
    phi = sum(phi_coefficients[i] * Z**i for i in range(3))
    for psi_coefficients in product(values, repeat=3):
        psi = sum(psi_coefficients[i] * Z**i for i in range(3))
        if is_zero(phi) and is_zero(psi):
            continue
        Q = sp.expand(
            X * phi.subs(Z, X * T) + T * psi.subs(Z, X * T)
        )
        ideal = sp.groebner((sp.diff(Q, X), sp.diff(Q, T)), X, T, order="lex")
        nonsingular = ideal.contains(sp.Integer(1))
        predicted = predicted_nonsingular(phi, psi)
        gate(nonsingular == predicted, "degree-two full radial classification")
        grid_total += 1
        grid_nonsingular += int(nonsingular)
        grid_predicted += int(predicted)

# Explicit smooth mixed-axis survivors and their monomial-Broughton no-mate
# controls.  The all-degree no-mate assertion is inherited from THM-3716;
# finite systems guard scalars, both orientations, and low degrees.
mixed_survivors = 0
mate_obstructions = 0
for n in range(2, 9):
    for a_value, b_value in ((1, 1), (2, -3), (-1, 2)):
        candidates = (
            a_value * X + b_value * X**n * T ** (n + 1),
            a_value * X ** (n + 1) * T**n + b_value * T,
        )
        for Q in candidates:
            ideal = sp.groebner((sp.diff(Q, X), sp.diff(Q, T)), X, T, order="lex")
            gate(ideal.contains(sp.Integer(1)), "mixed monomial smooth survivor")
            mixed_survivors += 1
            for bound in range(0, 8):
                basis = [
                    X**i * T**j
                    for i in range(bound + 1)
                    for j in range(bound + 1 - i)
                ]
                coefficients = sp.symbols(
                    f"mate_{n}_{mixed_survivors}_{bound}_0:{len(basis)}"
                )
                P = sum(
                    coefficient * monomial
                    for coefficient, monomial in zip(coefficients, basis)
                )
                equation = sp.Poly(jac(P, Q) - 1, X, T)
                system = [coefficient for _, coefficient in equation.terms()]
                gate(sp.linsolve(system, coefficients) == sp.EmptySet,
                     "mixed monomial finite no-mate")
                mate_obstructions += 1

# The n=1 axis boundary is singular, while n=0 is linear.  These are the
# equality/failure flanks around the n>=2 mixed family.
for a_value, b_value in ((1, 1), (2, -3), (-1, 2)):
    for Q in (
        a_value * X + b_value * X * T**2,
        a_value * X**2 * T + b_value * T,
    ):
        ideal = sp.groebner((sp.diff(Q, X), sp.diff(Q, T)), X, T, order="lex")
        gate(not ideal.contains(sp.Integer(1)), "mixed monomial n=1 singular flank")

# Pure-charge smooth and singular controls.
pure_profiles = (
    (1 + Z, True),
    (1 + Z + Z**2, True),
    ((1 - Z) * (1 + 2 * Z), True),
    ((1 - Z) ** 2, False),
    (Z * (1 + Z), False),
)
pure_controls = 0
for profile_value, expected in pure_profiles:
    for Q in (
        X * profile_value.subs(Z, X * T),
        T * profile_value.subs(Z, X * T),
    ):
        ideal = sp.groebner((sp.diff(Q, X), sp.diff(Q, T)), X, T, order="lex")
        gate(ideal.contains(sp.Integer(1)) == expected,
             "pure-charge squarefree classification")
        pure_controls += 1

# Semantic digest covers every structural alternative rather than the finite
# grid enumeration order.
semantic_rows = (
    "linear:phi,psi_constant_not_both_zero",
    "pure:Xphi_or_Tpsi;origin_nonzero;squarefree",
    "mixed:aX+bX^nT^(n+1)_or_dual;n>=2",
    "mate:linear_only;single_charge_antiderivative;THM3716_monomial_tail",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3741-radial-two-charge-Keller-component-classification")
print("scope=algebraically_closed_characteristic_zero;arbitrary_phi,psi_in_k[z]")
print("nonsingular=linear_or_squarefree_pure_charge_or_n>=2_mixed_monomial_flank")
print("jacobian_mate=iff_Q_is_nonzero_linear")
print("mechanisms=torus_eliminant;axis_valuation;charge_antiderivative;THM3716")
print(
    f"degree_two_grid={grid_total};nonsingular={grid_nonsingular};"
    f"predicted={grid_predicted}"
)
print(
    f"mixed_survivors={mixed_survivors};mate_obstructions={mate_obstructions};"
    f"pure_controls={pure_controls}"
)
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
