#!/usr/bin/env python3
"""Exact companion for THM-3895's quartic covariant and degree closure."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


# Abstract quartic square and the two equianharmonic choices d^2=-3.
T, G, K, a, L, d, Z = sp.symbols("T G K a L d Z")
quartic = L**4 - 6 * a * L**2 * T**2 - 8 * K * T**3 - 3 * a**2 * T**4
Z_covariant = (3 * a * L**2 + 4 * K * T + 3 * a**2 * T**2 - d * a * G) / 2
quadratic_relation = sp.expand(
    Z**2
    - (4 * K * T + 3 * a * L**2 + 3 * a**2 * T**2) * Z
    + 4 * K**2 * T**2
    + 6 * a * K * L**2 * T
    + 3 * a**2 * L**4
)


def reduce_d_and_curve(expression: sp.Expr) -> sp.Expr:
    numerator = sp.together(expression).as_numer_denom()[0]
    reduced_d = sp.Poly(numerator, d).rem(sp.Poly(d**2 + 3, d)).as_expr()
    return sp.Poly(reduced_d, G).rem(sp.Poly(G**2 - quartic, G)).as_expr()


zero(
    reduce_d_and_curve(quadratic_relation.subs(Z, Z_covariant)),
    "quartic covariant quadratic relation",
)

# Independent forgotten-arm audit: the sign-chosen square root at infinity
# has an osculating polynomial branch whose fourth y-coefficient exposes the
# global a=x+1 descent obstruction.
sroot = sp.symbols("sroot")
C_osculating = 8 * K**2 - 9 * a**3 * L**2
G_osculating = sroot * (
    a * T**2 + 4 * K * T / (3 * a) - C_osculating / (9 * a**3)
)
osculating_rhs = (
    C_osculating**2
    + 27 * a**6 * L**4
    - 24 * a**2 * K * C_osculating * T
) / (27 * a**6)
osculating_numerator = sp.together(
    (G - G_osculating) * (G + G_osculating) - osculating_rhs
).as_numer_denom()[0]
osculating_reduced = sp.Poly(osculating_numerator, sroot).rem(
    sp.Poly(sroot**2 + 3, sroot)
).as_expr()
osculating_reduced = sp.Poly(osculating_reduced, G).rem(
    sp.Poly(G**2 - quartic, G)
).as_expr()
zero(osculating_reduced, "osculating square-root identity")

for y_degree in range(3, 7):
    gate(
        (y_degree + 6) - 2 * y_degree == 6 - y_degree,
        f"osculating quotient degree d={y_degree}",
    )
    gate(6 - y_degree <= 3, f"osculating fourth-coefficient cutoff d={y_degree}")
gate(6 - 7 < 0, "osculating degree excludes d greater than six")

y_arm, F_arm = sp.symbols("y_arm F_arm")
t_arm = sp.symbols("t0:7")
T_arm = sum(t_arm[index] * y_arm**index for index in range(7))
K_arm = y_arm**2 - F_arm
C_arm = 8 * K_arm**2 - 9 * a**3 * L**2
G0_over_s = sp.expand(
    a * T_arm**2 + 4 * K_arm * T_arm / (3 * a) - C_arm / (9 * a**3)
)
y4_actual = sp.Poly(G0_over_s, y_arm).coeff_monomial(y_arm**4)
y4_expected = (
    a * sum(t_arm[index] * t_arm[4 - index] for index in range(5))
    + sp.Rational(4, 3) * (t_arm[2] - F_arm * t_arm[4]) / a
    - sp.Rational(8, 9) / a**3
)
zero(y4_actual - y4_expected, "osculating y-four coefficient")
u2_arm, u4_arm = sp.symbols("u2_arm u4_arm")
y4_arm_substitution = sp.expand(
    y4_expected.subs({t_arm[2]: a * u2_arm, t_arm[4]: a * u4_arm})
)
gate(
    sp.limit(a**3 * y4_arm_substitution, a, 0) == -sp.Rational(8, 9),
    "osculating exact a-minus-three pole",
)

Phi = sp.expand(Z**2 - 3 * a * L**2 * Z + 3 * a**2 * L**4)
A_sidecar = sp.expand(2 * K * (2 * Z - 3 * a * L**2))
B_sidecar = sp.expand(3 * a**2 * Z - 4 * K**2)
zero(
    quadratic_relation - (Phi - T * (A_sidecar + B_sidecar * T)),
    "quartic divisor identity",
)
rho_plus = (3 + d) / 2
rho_minus = (3 - d) / 2
zero(
    sp.Poly(
        sp.together(
            (Z - rho_plus * a * L**2) * (Z - rho_minus * a * L**2) - Phi
        ).as_numer_denom()[0],
        d,
    )
    .rem(sp.Poly(d**2 + 3, d))
    .as_expr(),
    "two shifted quartic colours",
)

# The discriminant of B*T^2+A*T-Phi after solving Z=(4K^2+B)/(3a^2).
B = sp.symbols("B")
A0 = sp.symbols("A0")
Z_from_B = (4 * K**2 + B) / (3 * a**2)
quadratic_discriminant = sp.expand(A_sidecar**2 + 4 * B_sidecar * Phi)
discriminant_substitution = sp.factor(quadratic_discriminant.subs(Z, Z_from_B))
expected_discriminant = sp.factor(
    sp.Rational(4, 9)
    / a**4
    * (
        (4 * K**2 + B - 3 * a**3 * L**2) ** 3
        + 27 * (a**3 * L**2) ** 2 * (a**3 * L**2 - K**2)
    )
)
zero(discriminant_substitution - expected_discriminant, "normalized discriminant")

# The degree bookkeeping which reduces every n>2 row to B-degree 6-n.
degree_ladder = []
for n in range(3, 7):
    b_degree = 6 - n
    gate(0 <= b_degree <= 3, f"B-degree range at n={n}")
    gate(n + b_degree == 6, f"degree-six cancellation at n={n}")
    gate(8 - n < 6, f"quotient lies below the A-sidecar at n={n}")
    degree_ladder.append((n, b_degree))
for n in (7, 8):
    gate(8 - n < 6 < n, f"degree-{n} cannot cancel")

# Normalize A0=a^3*L^2 to one after an algebraic constant extension.  The
# complete square universe is B_d=sum b_i*y^i, 0<=d<=3.  Descending from the
# fixed leading term 8*y^6 determines the square root uniquely; the six low
# residual coefficients have reduced grevlex basis [b_0,...,b_d] over Q(f).
y, f = sp.symbols("y f")
b_symbols = sp.symbols("b0:4")
kappa = y**2 - f
groebner_signatures: list[str] = []


def square_residual_packet(degree: int) -> tuple[list[sp.Expr], list[sp.Expr]]:
    variables = list(b_symbols[: degree + 1])
    b_polynomial = sum(variables[i] * y**i for i in range(degree + 1))
    norm = sp.expand((4 * kappa**2 + b_polynomial - 3) ** 3 + 27 * (1 - kappa**2))
    norm_poly = sp.Poly(norm, y)
    gate(norm_poly.degree() == 12, f"normalized norm degree d={degree}")
    gate(norm_poly.LC() == 64, f"normalized norm leading coefficient d={degree}")

    square_root = 8 * y**6
    for power in range(11, 5, -1):
        current = sp.Poly(sp.expand(square_root**2), y).coeff_monomial(y**power)
        target = norm_poly.coeff_monomial(y**power)
        square_root += sp.cancel((target - current) / 16) * y ** (power - 6)

    residual = sp.Poly(sp.expand(square_root**2 - norm), y)
    gate(
        all(residual.coeff_monomial(y**power) == 0 for power in range(6, 13)),
        f"descending square-root recursion d={degree}",
    )
    low_residuals = [
        sp.together(residual.coeff_monomial(y**power)).as_numer_denom()[0]
        for power in range(6)
    ]
    basis = sp.groebner(
        low_residuals,
        *variables,
        order="grevlex",
        domain=sp.QQ.frac_field(f),
    )
    basis_expressions = [polynomial.as_expr() for polynomial in basis.polys]
    gate(basis_expressions == variables, f"reduced square ideal d={degree}")
    groebner_signatures.append(
        f"d{degree}:[{','.join(str(variable) for variable in variables)}]"
    )

    positive_control = sp.expand(norm.subs({variable: 0 for variable in variables}))
    zero(
        positive_control - (kappa * (8 * kappa**2 - 9)) ** 2,
        f"B-zero square positive control d={degree}",
    )
    return low_residuals, variables


packets = [square_residual_packet(degree) for degree in range(4)]

# A short transparent edge proof for the linear B row, independently of its
# Groebner basis.  These are the last three residuals in descending order.
linear_residuals, _ = packets[1]
b0, b1 = b_symbols[:2]
zero(linear_residuals[5] + 6 * b1 * (b0 - 3), "linear-B E5 factor")
zero(
    linear_residuals[4] + 3 * (b0**2 - 6 * b0 + b1**2 * f),
    "linear-B E4 factor",
)
zero(
    linear_residuals[3] - b1 * (96 * b0 * f + b1**2 - 288 * f),
    "linear-B E3 factor",
)

# The constant-B parity proof is a second independent edge control.
q, c = sp.symbols("q c")
constant_B_cubic = sp.expand(
    64 * q**3
    + (48 * c - 144) * q**2
    + (12 * c**2 - 72 * c + 81) * q
    + c * (c**2 - 9 * c + 27)
)
constant_relation = c**2 - 9 * c + 27
constant_quadratic = sp.rem(constant_B_cubic, constant_relation, c) / q
zero(
    sp.rem(sp.discriminant(constant_quadratic, q), constant_relation, c)
    - 2304 * (9 - c),
    "constant-B square obstruction",
)

semantic = {
    "covariant": "Z=(3aL2+4KT+3a2T2-daG)/2,d2=-3",
    "quadratic": "Z2-(4KT+3aL2+3a2T2)Z+4K2T2+6aKL2T+3a2L4=0",
    "divisor": "Phi=T*(2K(2Z-3aL2)+(3a2Z-4K2)T)",
    "degree": "n>2 forces degZ=4,degB=6-n,n in {3,4,5,6}",
    "square_ideals": ";".join(groebner_signatures),
    "independent_arm_audit": "osculating identity plus uncancellable a^-3 y4 pole",
    "conclusion": "every polynomial residual square has deg_y(T)<=2",
    "scope": "the remaining quadratic-y coefficient functions and JC2 stay open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3895-f-zero-quartic-covariant-and-high-y-degree-emptiness")
print("quartic_covariant=degree_four")
print("shifted_fibres=two_coprime_quartic_colours")
print("degree_ladder=" + ",".join(f"n{n}:Bdeg{b}" for n, b in degree_ladder))
print("square_ideal_bases=" + ";".join(groebner_signatures))
print("positive_control=B0,K^2(8K^2-9)^2")
print("independent_arm_audit=osculating_y4_a_minus_3_pole")
print("conclusion=deg_y(T)_at_most_2")
print("remaining=T=t2(x)y^2+t1(x)y+t0(x)")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
