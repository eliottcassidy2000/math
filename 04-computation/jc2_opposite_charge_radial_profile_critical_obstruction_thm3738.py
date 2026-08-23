#!/usr/bin/env python3
"""Exact companion for THM-3738's opposite-charge critical obstruction."""

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
phi_function = sp.Function("phi")
psi_function = sp.Function("psi")
phi = phi_function(Z)
psi = psi_function(Z)
A = sp.diff(Z * phi, Z)
B = sp.diff(Z * psi, Z)
C = sp.diff(Z * phi * psi, Z)

# Universal elimination identity.  On the torus XT=z, multiplying the two
# critical equations by X and T gives this 2x2 matrix; its determinant is
# exactly the derivative of z*phi*psi.
critical_matrix = sp.Matrix(
    ((A, Z * sp.diff(psi, Z)), (Z * sp.diff(phi, Z), B))
)
gate(sp.simplify(critical_matrix.det() - C) == 0,
     "universal torus critical determinant")

semantic_rows: list[str] = []


def audit_profiles(phi_value: sp.Expr, psi_value: sp.Expr, label: str) -> None:
    """Check the proof's gcd count and an exact critical-point certificate."""
    P = sp.expand(phi_value * psi_value)
    C_value = sp.diff(Z * P, Z)
    P_poly = sp.Poly(P, Z)
    C_poly = sp.Poly(C_value, Z)
    repeated = sp.gcd(P_poly, sp.diff(P_poly, Z))
    common = sp.gcd(P_poly, C_poly)
    squarefree_degree = P_poly.degree() - repeated.degree()
    gate(phi_value.subs(Z, 0) != 0 and psi_value.subs(Z, 0) != 0,
         f"nonzero axis gradients: {label}")
    gate(common.monic() == repeated.monic(), f"gcd identity: {label}")
    gate(C_poly.degree() == P_poly.degree(), f"eliminant degree: {label}")
    quotient = sp.cancel(C_poly.as_expr() / common.as_expr())
    gate(sp.Poly(quotient, Z).degree() == squarefree_degree,
         f"off-product critical-root count: {label}")
    gate(squarefree_degree >= 1, f"positive off-product count: {label}")

    # Over the algebraic closure the quotient has a root z0.  The theorem
    # turns a non-axis kernel vector into a critical point by a common square
    # root scaling.  For the exact companion, independently verify the same
    # conclusion by the affine critical ideal.
    Q = sp.expand(
        X * phi_value.subs(Z, X * T)
        + T * psi_value.subs(Z, X * T)
    )
    ideal = sp.groebner((sp.diff(Q, X), sp.diff(Q, T)), X, T, order="lex")
    gate(not ideal.contains(sp.Integer(1)), f"critical ideal nonempty: {label}")
    semantic_rows.append(
        f"{label}:"
        + hashlib.sha256(
            "|".join(sp.srepr(item) for item in (phi_value, psi_value, quotient)).encode()
        ).hexdigest()
    )


# Hostile profiles include simple, repeated, and shared roots.  They guard the
# precise point at which a root of the eliminant could otherwise yield only
# an axis kernel.
hostile_profiles = (
    (1 + Z, 1 - Z, "linear-opposite"),
    ((1 - Z) ** 4, (1 + 2 * Z) ** 3, "separate-multiple-roots"),
    ((1 - Z) ** 3 * (1 + Z) ** 2,
     (1 - Z) ** 2 * (1 + 3 * Z) ** 4,
     "shared-multiple-root"),
    (1 + Z + Z**3, 2 - 3 * Z + Z**4, "sparse-asymmetric"),
    ((1 + 2 * Z) ** 8, -(1 - 2 * Z) ** 7, "deep-binomial-pair"),
)
for phi_value, psi_value, label in hostile_profiles:
    audit_profiles(sp.expand(phi_value), sp.expand(psi_value), label)

# Exhaust the normalized degree-at-most-two integer box used to discover the
# obstruction.  Every nontrivial mixed profile has a nonempty critical ideal.
grid_total = 0
grid_critical = 0
values = (-2, -1, 0, 1, 2)
for psi_constant in (-2, -1, 1, 2):
    for p1, p2, q1, q2 in product(values, repeat=4):
        if p1 == p2 == q1 == q2 == 0:
            continue
        phi_value = 1 + p1 * Z + p2 * Z**2
        psi_value = psi_constant + q1 * Z + q2 * Z**2
        Q = sp.expand(
            X * phi_value.subs(Z, X * T)
            + T * psi_value.subs(Z, X * T)
        )
        ideal = sp.groebner((sp.diff(Q, X), sp.diff(Q, T)), X, T, order="lex")
        grid_total += 1
        gate(not ideal.contains(sp.Integer(1)), "degree-two grid critical point")
        grid_critical += 1

# The theorem's equality boundary consists exactly of nonzero constant
# profiles, giving a linear coordinate with constant nonzero gradient.
linear_controls = 0
for phi_constant in (-3, -1, 1, 2):
    for psi_constant in (-2, -1, 1, 4):
        Q = phi_constant * X + psi_constant * T
        ideal = sp.groebner((sp.diff(Q, X), sp.diff(Q, T)), X, T, order="lex")
        gate(ideal.contains(sp.Integer(1)), "linear no-critical boundary")
        linear_controls += 1

# Every pair of lower/upper THM-3734 divided-power profiles is a special case,
# even at unequal depths.  Check a rectangular exact sample independently.
binomial_pairs = 0
for r_lower in range(1, 7):
    phi_lower = sp.cancel(
        r_lower * ((1 + 2 * Z / r_lower) ** (r_lower + 1) - 1)
        / (2 * (r_lower + 1) * Z)
    )
    for r_upper in range(1, 7):
        phi_upper = sp.cancel(
            -r_upper * (1 - (1 - 2 * Z / r_upper) ** (r_upper + 1))
            / (2 * (r_upper + 1) * Z)
        )
        Q = sp.expand(
            X * phi_lower.subs(Z, X * T)
            + T * phi_upper.subs(Z, X * T)
        )
        ideal = sp.groebner((sp.diff(Q, X), sp.diff(Q, T)), X, T, order="lex")
        gate(not ideal.contains(sp.Integer(1)),
             "unequal-depth binomial mixture critical point")
        binomial_pairs += 1

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()
print("theorem=THM-3738-opposite-charge-radial-profile-critical-obstruction")
print("scope=algebraically_closed_characteristic_zero;phi(0)psi(0)!=0")
print("classification=Q=Xphi(XT)+Tpsi(XT)_nonsingular_iff_phi,psi_constant")
print("eliminant=det_torus_critical_matrix=(z_phi_psi)prime")
print("mechanism=degree_gcd_count_forces_off-product_root_with_nonaxis_kernel")
print(f"degree_two_grid={grid_critical}/{grid_total}_critical")
print(f"linear_controls={linear_controls}_nonsingular")
print(f"binomial_pairs={binomial_pairs}_critical")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
