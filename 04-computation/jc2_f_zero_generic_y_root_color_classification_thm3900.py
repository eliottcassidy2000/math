#!/usr/bin/env python3
"""Focused exact companion for THM-3900's generic-y root-color classification.

The only theorem-level input is THM-3895's proved generic-y cutoff.  Every
remaining degree-two, degree-one, and degree-zero calculation is replayed
here, including the fifteen constant-T rational-root candidates.  No
elliptic-surface, Mordell--Weil, or integrality claim from THM-3888 is used.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


# ---------------------------------------------------------------------------
# 1. Exact normalization and the two known generic solutions.
# ---------------------------------------------------------------------------

u, z, g, scale, L0 = sp.symbols("u z g scale L0", nonzero=True)
a0 = scale**2
T0 = L0 * u / scale
K0 = scale**3 * L0 * z
G0 = L0**2 * g
quartic0 = sp.expand(
    L0**4 - 6 * a0 * L0**2 * T0**2 - 8 * K0 * T0**3 - 3 * a0**2 * T0**4
)
zero(
    quartic0 - L0**4 * (1 - 6 * u**2 - 8 * z * u**3 - 3 * u**4),
    "normalized f-zero quartic",
)
zero(G0**2 - L0**4 * g**2, "normalized square root")

normalized_rhs = sp.expand(1 - 6 * u**2 - 8 * z * u**3 - 3 * u**4)
zero(
    normalized_rhs - 1 + u**2 * (6 + 8 * z * u + 3 * u**2),
    "root-color factor debt",
)

u_star = -sp.Rational(2, 3) * z
g_star = sp.Rational(4, 3) * z**2 - 1
zero(g_star**2 - normalized_rhs.subs(u, u_star), "normalized hostile T-star")
zero(normalized_rhs.subs(u, 0) - 1, "normalized T-zero branch")

x, y = sp.symbols("x y")
a = x + 1
L = 9 * x + 4
F = 15 * x**2 + 15 * x + 4
K = y**2 - F
T_star = -2 * K / (3 * a**2)
G_star = 4 * K**2 / (3 * a**3) - L**2
quartic = sp.expand(L**4 - 6 * a * L**2 * sp.Symbol("T")**2
                    - 8 * K * sp.Symbol("T")**3
                    - 3 * a**2 * sp.Symbol("T")**4)
zero(G_star**2 - quartic.subs(sp.Symbol("T"), T_star),
     "original-coordinate hostile T-star")


# ---------------------------------------------------------------------------
# 2. Same-color quadratic roots force T-star.
# ---------------------------------------------------------------------------

epsilon, color_scalar = sp.symbols("epsilon color_scalar")
same_color_g = epsilon + color_scalar * u**2
same_color_debt = sp.expand(
    same_color_g**2 - normalized_rhs
)
same_color_reduced = sp.Poly(same_color_debt, epsilon).rem(
    sp.Poly(epsilon**2 - 1, epsilon)
).as_expr()
zero(
    same_color_reduced
    - u**2
    * (
        (color_scalar**2 + 3) * u**2
        + 8 * z * u
        + 2 * epsilon * color_scalar
        + 6
    ),
    "same-color scalar equation",
)
zero(
    same_color_reduced.subs(
        {color_scalar: -3 * epsilon, u: -sp.Rational(2, 3) * z},
        simultaneous=True,
    ).subs(epsilon**2, 1),
    "same-color solution is T-star",
)


# ---------------------------------------------------------------------------
# 3. Split colors: Hermite response, exact ideal, and sharp hostile.
# ---------------------------------------------------------------------------

Y = sp.symbols("Y")
A, q, midpoint, correction, fpar = sp.symbols(
    "A q midpoint correction fpar", nonzero=True
)
X = sp.symbols("X")
u_split = A * (X**2 - q**2)
z_split = (X + midpoint) ** 2 - fpar
g_hermite = 3 * X / (2 * q) - X**3 / (2 * q**3)
g_split = sp.expand(g_hermite + correction * u_split**2)

zero(g_hermite.subs(X, q) - 1, "Hermite color plus")
zero(g_hermite.subs(X, -q) + 1, "Hermite color minus")
zero(sp.diff(g_hermite, X).subs(X, q), "Hermite derivative plus")
zero(sp.diff(g_hermite, X).subs(X, -q), "Hermite derivative minus")

split_residual = sp.Poly(
    sp.expand(
        sp.together(
            g_split**2
            - (1 - 6 * u_split**2 - 8 * z_split * u_split**3 - 3 * u_split**4)
        )
        * 4
        * q**6
    ),
    X,
)

coefficient_7 = sp.factor(split_residual.coeff_monomial(X**7))
coefficient_5 = sp.factor(split_residual.coeff_monomial(X**5))
coefficient_3 = sp.factor(split_residual.coeff_monomial(X**3))
zero(
    coefficient_7 - 4 * A**2 * q**3 * (16 * A * q**3 * midpoint - correction),
    "split degree-seven response",
)
zero(
    coefficient_5
    + 4 * A**2 * q**5 * (48 * A * q**3 * midpoint - 5 * correction),
    "split degree-five response",
)
zero(
    coefficient_3
    - 4 * A**2 * q**7 * (48 * A * q**3 * midpoint - 7 * correction),
    "split degree-three response",
)

# The first two displayed responses imply midpoint=correction=0 because
# A*q is nonzero.  Substitute that forced branch and read the top response.
split_centered = sp.Poly(
    sp.expand(split_residual.as_expr().subs({midpoint: 0, correction: 0})), X
)
zero(
    split_centered.coeff_monomial(X**8)
    - 4 * A**3 * q**6 * (3 * A + 8),
    "split top response forces A=-8/3",
)

A_split = -sp.Rational(8, 3)
split_final = sp.Poly(
    sp.expand(split_centered.as_expr().subs(A, A_split)), X
)
even_equations = [
    sp.together(split_final.coeff_monomial(X**power)).as_numer_denom()[0]
    for power in (6, 4, 2, 0)
]
split_basis = sp.groebner(
    even_equations,
    fpar,
    q,
    order="lex",
    domain=sp.QQ,
)
split_basis_expressions = [polynomial.as_expr() for polynomial in split_basis.polys]
gate(
    split_basis_expressions
    == [fpar + sp.Rational(13, 3) * q**2, q**4 - sp.Rational(9, 512)],
    "complete split-color response ideal",
)

C_split = sp.Rational(8, 3) * q**2
zero(
    sp.Poly(C_split**2 - sp.Rational(1, 8), q).rem(
        sp.Poly(q**4 - sp.Rational(9, 512), q)
    ).as_expr(),
    "split C squared",
)
zero(
    sp.Poly(
        (-sp.Rational(13, 3) * q**2) ** 2 - sp.Rational(169, 512), q
    ).rem(sp.Poly(q**4 - sp.Rational(9, 512), q)).as_expr(),
    "split f squared",
)

u_sharp = A_split * (X**2 - q**2)
f_sharp = -sp.Rational(13, 3) * q**2
g_sharp = 3 * X / (2 * q) - X**3 / (2 * q**3)
sharp_numerator = sp.together(
    g_sharp**2
    - (
        1
        - 6 * u_sharp**2
        - 8 * (X**2 - f_sharp) * u_sharp**3
        - 3 * u_sharp**4
    )
).as_numer_denom()[0]
zero(
    sp.Poly(sharp_numerator, q).rem(
        sp.Poly(q**4 - sp.Rational(9, 512), q)
    ).as_expr(),
    "sharp split-color hostile",
)

# In the actual family f^2=F^2/(a^3*L^2), which is not 169/512.
actual_split_mismatch = sp.expand(512 * F**2 - 169 * a**3 * L**2)
gate(sp.Poly(actual_split_mismatch, x).degree() == 5,
     "actual split mismatch degree")
gate(sp.Poly(actual_split_mismatch, x).LC() == -13689,
     "actual split mismatch leading coefficient")


# ---------------------------------------------------------------------------
# 4. Linear u has odd degree; constant u has fifteen rational-root cells.
# ---------------------------------------------------------------------------

linear_A, linear_B = sp.symbols("linear_A linear_B", nonzero=True)
u_linear = linear_A * Y + linear_B
linear_rhs = sp.Poly(
    sp.expand(
        1
        - 6 * u_linear**2
        - 8 * (Y**2 - fpar) * u_linear**3
        - 3 * u_linear**4
    ),
    Y,
)
gate(linear_rhs.degree() == 5, "linear channel has degree five")
zero(linear_rhs.LC() + 8 * linear_A**3, "linear odd-degree leading form")

tconst, croot = sp.symbols("tconst croot")
quartic_tconst = sp.Poly(
    sp.expand(
        L**4
        - 6 * a * L**2 * tconst**2
        - 8 * K * tconst**3
        - 3 * a**2 * tconst**4
    ),
    y,
)
zero(quartic_tconst.coeff_monomial(y**2) + 8 * tconst**3,
     "constant-T y-square coefficient")
zero(quartic_tconst.coeff_monomial(y), "constant-T missing y coefficient")
C_tconst = sp.expand(
    L**4 - 6 * a * L**2 * tconst**2 + 8 * F * tconst**3 - 3 * a**2 * tconst**4
)
zero(quartic_tconst.coeff_monomial(1) - C_tconst,
     "constant-T residual constant")
zero(C_tconst.subs(tconst, 0) - L**4, "constant-T zero boundary")

C_coefficients = [
    sp.Poly(coefficient, x, domain=sp.QQ)
    for coefficient in sp.Poly(C_tconst, tconst).all_coeffs()
]
C_content = C_coefficients[0]
for coefficient in C_coefficients[1:]:
    C_content = C_content.gcd(coefficient)
gate(C_content.degree() == 0, "constant-T polynomial primitive in k[x]")
gate(sp.gcd(a, L) == 1, "rational-root numerator and denominator primes")

rational_root_candidates: list[tuple[int, int, sp.Expr]] = []
for numerator_power in range(5):
    for denominator_power in range(3):
        candidate = croot * L**numerator_power / a**denominator_power
        candidate_numerator = sp.expand(
            sp.together(C_tconst.subs(tconst, candidate)).as_numer_denom()[0]
        )
        coefficient_expressions = [
            coefficient
            for coefficient in sp.Poly(candidate_numerator, x).all_coeffs()
            if coefficient != 0
        ]
        coefficient_gcd = sp.Poly(
            coefficient_expressions[0], croot, domain=sp.QQ
        )
        for coefficient in coefficient_expressions[1:]:
            coefficient_gcd = coefficient_gcd.gcd(
                sp.Poly(coefficient, croot, domain=sp.QQ)
            )
        gate(
            coefficient_gcd.degree() == 0,
            f"constant-T rational-root cell ({numerator_power},{denominator_power})",
        )
        rational_root_candidates.append(
            (
                numerator_power,
                denominator_power,
                sp.monic(coefficient_gcd.as_expr(), croot),
            )
        )
gate(len(rational_root_candidates) == 15, "all constant-T cells exhausted")


semantic = {
    "theorem": "THM-3900",
    "dependency": "THM3895 generic-y cutoff only",
    "normalized": "g^2=1-6u^2-8zu^3-3u^4,z=Y^2-f",
    "quadratic_same": "u=-2z/3 gives T-star",
    "quadratic_split": "A=-8/3,B=0,C^2=1/8,f^2=169/512",
    "split_hostile": "exact over the sharp constant f, excluded by actual nonconstant f",
    "linear": "unique odd degree five",
    "constant": "15 rational-root cells all unit gcd",
    "classification": "T=0 or T=-2K/(3a^2)",
    "elliptic_dependency": "none",
    "scope": "generic k(x)[y] f-zero chart; global polynomial descent separate",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3900-f-zero-generic-y-root-color-classification")
print("status=RESERVED_PROVISIONAL_PROOF_CANDIDATE_VERIFIED_EXACT_AWAITING_INDEPENDENT_HOSTILE_AUDIT")
print("dependency=THM-3895_generic_y_cutoff")
print("normalized_equation=g^2=1-6u^2-8zu^3-3u^4")
print("same_color=T_star")
print("split_response=A=-8/3;B=0;C^2=1/8;f^2=169/512")
print("sharp_split_hostile=EXISTS_IN_CONSTANT_FAMILY")
print("actual_split_family=EXCLUDED_BY_NONCONSTANT_f^2")
print("linear_channel=EMPTY_ODD_DEGREE_FIVE")
print("constant_T_candidates=15_ALL_EMPTY")
print(
    "constant_T_candidate_gcds="
    + ",".join(
        f"({i},{j}):{gcd_value}"
        for i, j, gcd_value in rational_root_candidates
    )
)
print("classification=T_ZERO_OR_T_STAR")
print("elliptic_dependencies=NONE")
print("global_polynomial_descent=SEPARATE")
print("JC2_status=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
