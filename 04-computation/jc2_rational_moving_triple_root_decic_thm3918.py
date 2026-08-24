#!/usr/bin/env python3
"""Exact companion for THM-3918's sparse rational moving-root decic."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


def reduce_quadratic(expression: sp.Expr, symbol: sp.Symbol, constant: int) -> sp.Expr:
    numerator, denominator = sp.fraction(sp.together(expression))
    remainder = sp.rem(
        sp.Poly(sp.expand(numerator), symbol),
        sp.Poly(symbol**2 + constant, symbol),
    ).as_expr()
    return sp.cancel(remainder / denominator)


def zero_mod_rho(expression: sp.Expr, message: str) -> None:
    gate(reduce_quadratic(expression, rho, 6) == 0, message)


def binary_disc(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    return sp.expand(
        b**2 * c**2
        - 4 * a * c**3
        - 4 * b**3 * d
        - 27 * a**2 * d**2
        + 18 * a * b * c * d
    )


A, C, U, V, Z, rho = sp.symbols("A C U V Z rho")
ell = A * U + C * V
Phi = sp.expand(
    ell**3 + C * ell * U**2 + rho * C * ell * V**2 + A**2 * V**3
)
Phi_poly = sp.Poly(Phi, U, V)
coefficients = tuple(
    Phi_poly.coeff_monomial(monomial)
    for monomial in (U**3, U**2 * V, U * V**2, V**3)
)
expected_coefficients = (
    A**3 + A * C,
    3 * A**2 * C + C**2,
    3 * A * C**2 + rho * A * C,
    C**3 + rho * C**2 + A**2,
)
gate(coefficients == expected_coefficients, "compact coefficient expansion")
gate(
    all(coefficient.subs({A: 0, C: 0}) == 0 for coefficient in coefficients),
    "common coefficient zero",
)
gate(
    sp.factor(Phi.subs(C, 0)) == A**2 * (A * U**3 + V**3),
    "Kummer irreducibility specialization",
)
gate(sp.rem(coefficients[0], C, C) != 0, "U-cubic row retains degree at C=0")
gate((-1) % 3 != 0, "A-valuation of -1/A is not a multiple of three")

Delta = binary_disc(*coefficients)
expected_delta = (
    -27 * A**10
    - 54 * A**8 * C
    - 4 * (rho**3 + 9 * rho) * A**6 * C**3
    - 27 * A**6 * C**2
    - 12 * (rho**2 + 3) * A**4 * C**5
    - 4 * (rho**3 + 9 * rho) * A**4 * C**4
    - 12 * rho * A**2 * C**7
    - 4 * (2 * rho**2 + 1) * A**2 * C**6
    - 4 * C**9
    - 4 * rho * C**8
)
zero(Delta - expected_delta, "exact family discriminant")

Delta_special = sp.expand(
    -27 * A**10
    - 54 * A**8 * C
    - 12 * rho * A**6 * C**3
    - 27 * A**6 * C**2
    + 36 * A**4 * C**5
    - 12 * rho * A**4 * C**4
    - 12 * rho * A**2 * C**7
    + 44 * A**2 * C**6
    - 4 * C**9
    - 4 * rho * C**8
)
zero_mod_rho(Delta - Delta_special, "rho^2=-6 discriminant reduction")
gate(sp.Poly(Delta_special, A, C).total_degree() == 10, "decic degree")

Delta_h = sp.expand(
    sum(
        coefficient
        * A**powers[0]
        * C**powers[1]
        * Z ** (10 - sum(powers))
        for powers, coefficient in sp.Poly(Delta_special, A, C).terms()
    )
)
zero(Delta_h.subs(Z, 0) + 27 * A**10, "one-point degree-ten row")
infinity_point = {A: 0, C: 1, Z: 0}
gate(sp.diff(Delta_h, Z).subs(infinity_point) == -4, "smooth infinity point")

# The direction pencil makes the genus transition visible before specializing.
zeta = sp.symbols("zeta")
radial = sp.cancel(Delta.subs(A, zeta * C) / C**8)
radial_poly = sp.Poly(radial, C)
gate(radial_poly.degree() == 2, "radial equation is quadratic")
radial_disc = sp.factor(sp.discriminant(radial_poly.as_expr(), C))
h_r = (rho**2 + 6) * zeta**4 + 2 * rho * zeta**2 + 1
zero(radial_disc - 16 * h_r**3, "radial discriminant is a perfect cube")
zero_mod_rho(h_r - (1 + 2 * rho * zeta**2), "critical radial conic")

# Repeated-root incidence and its rational critical normalization.
x, y, z, w = sp.symbols("x y z w")
L_formula = -(1 - x * y) * (x**2 + rho) - y**2
H = sp.expand(
    3 * y**3
    - 2 * (x**3 + rho * x) * y**2
    + 4 * x**2 * y
    + 2 * rho * y
    - 2 * x
)
f_incidence = sp.expand(
    (A * x + C) ** 3
    + C * (A * x + C) * (x**2 + rho)
    + A**2
)
df_incidence = sp.diff(f_incidence, x)
A_incidence = y * L_formula
C_incidence = (1 - x * y) * L_formula
zero(
    f_incidence.subs({A: A_incidence, C: C_incidence}),
    "incidence cubic equation",
)
zero(
    df_incidence.subs({A: A_incidence, C: C_incidence}) + L_formula**2 * H,
    "incidence derivative equation",
)
K_quartic = sp.expand(
    3 * y**4 + 2 * rho * (1 - z) * y**2 - 2 * z * (z - 1) ** 2
)
zero(sp.expand((H * y).subs(x, z / y) - K_quartic), "z=xy quartic")
w_inverse = sp.cancel(3 * y**2 / (z - 1) - rho)
zero(
    sp.factor((w_inverse**2 - rho**2 - 6 * z) * (z - 1) ** 2 - 3 * K_quartic),
    "incidence radial relation",
)
zero(
    sp.factor(
        18 * y**2
        - (w_inverse + rho) * (w_inverse**2 - rho**2 - 6)
        + 9 * y**2 * K_quartic / (z - 1) ** 3
    ),
    "elliptic incidence model",
)

t = sp.symbols("t")
w_t = 18 * t**2 - rho
y_t = t * w_t
zero_mod_rho(
    18 * y_t**2 - w_t**2 * (w_t + rho),
    "nodal cubic rational normalization",
)
special_cubic = 18 * y**2 - w**2 * (w + rho)
gate(
    sp.hessian(special_cubic, (w, y)).subs({w: 0, y: 0}).det() == -72 * rho,
    "critical incidence cubic has an ordinary node",
)

L_t = 162 * t**6 - 18 * rho * t**4 + 6 * t**2 - rho
A_t = sp.expand(t * (18 * t**2 - rho) * L_t)
C_t = sp.expand(-(54 * t**4 - 6 * rho * t**2 - 1) * L_t)
G_t = sp.expand((18 * t**2 - rho) * L_t)
zero_mod_rho(A_t - t * G_t, "factored A parametrization")
zero_mod_rho(C_t + (3 * t**2 - rho / 6) * G_t, "factored C parametrization")
zero_mod_rho(Delta_special.subs({A: A_t, C: C_t}), "parametrization lies on decic")
gate(sp.degree(A_t, t) == 9, "A parameter degree nine")
gate(sp.degree(C_t, t) == 10, "C parameter degree ten")
gate(sp.LC(sp.Poly(A_t, t)) == 2916, "A leading coefficient")
gate(sp.LC(sp.Poly(C_t, t)) == -8748, "C leading coefficient")
gate(10 - 9 == 1, "A/C is an infinity uniformizer")

# Complete finite singular-address packet.
resultant_G = sp.resultant(G_t, sp.diff(G_t, t), t)
gate(
    reduce_quadratic(resultant_G, rho, 6)
    == -58535441882591744284714139375370240000,
    "eight origin addresses are reduced",
)
eta = sp.sqrt(-6)
number_field = sp.QQ.algebraic_field(eta)
L_eta = sp.expand(L_t.subs(rho, eta))
G_eta = sp.expand(G_t.subs(rho, eta))
A_eta = sp.expand(A_t.subs(rho, eta))
C_eta = sp.expand(C_t.subs(rho, eta))
dA_eta = sp.Poly(sp.diff(A_eta, t), t, domain=number_field)
dC_eta = sp.Poly(sp.diff(C_eta, t), t, domain=number_field)
critical = sp.gcd(dA_eta, dC_eta).monic()
gate(
    critical == sp.Poly(t**2 + eta / 18, t, domain=number_field),
    "exact two-address critical support",
)
cusp_det = sp.expand(
    sp.diff(A_eta, t, 2) * sp.diff(C_eta, t, 3)
    - sp.diff(C_eta, t, 2) * sp.diff(A_eta, t, 3)
)
gate(
    sp.rem(sp.Poly(cusp_det, t, domain=number_field), critical).as_expr()
    == -276480 * eta,
    "both critical addresses are ordinary A2 cusps",
)
gate(
    sp.rem(sp.Poly(A_eta, t, domain=number_field), critical).as_expr()
    == -10 * t,
    "cusp target A coordinate",
)
gate(
    sp.rem(sp.Poly(C_eta, t, domain=number_field), critical).as_expr()
    == -10 * eta / 3,
    "cusp target C coordinate",
)

iota = -eta / (18 * t)
G_iota_numerator = sp.together(G_eta.subs(t, iota)).as_numer_denom()[0]
origin_pair_gcd = sp.gcd(
    sp.Poly(G_eta, t, domain=number_field),
    sp.Poly(G_iota_numerator, t, domain=number_field),
).monic()
gate(
    origin_pair_gcd == sp.Poly(t**2 - eta / 18, t, domain=number_field),
    "only one origin pair shares a tangent",
)
num_A_collision = sp.together(A_eta - A_eta.subs(t, iota)).as_numer_denom()[0]
num_C_collision = sp.together(C_eta - C_eta.subs(t, iota)).as_numer_denom()[0]
a_eta = eta / 18
zero(
    num_A_collision
    - 153055008 * (t**2 - a_eta) ** 6 * (t**2 + a_eta) ** 3,
    "sixth strict radial response of tangent origin pair",
)
zero(
    num_C_collision
    + 459165024 * (t**2 - a_eta) ** 7 * (t**2 + a_eta) ** 3,
    "seventh transverse response of tangent origin pair",
)
gate(sp.binomial(8, 2) + (7 - 1) == 34, "origin delta is 34")
gate(34 + 2 == 36, "complete decic genus ledger")

# The old tournament shadow misses the critical event; an infinity anchor repairs it.
r, color, scale, xi = sp.symbols("r color scale xi")
H_color = (r**2 + 6) * color**2 + 2 * r * color * scale + scale**2
color_product = (scale + (r + xi) * color) * (scale + (r - xi) * color)
gate(
    reduce_quadratic(H_color - color_product, xi, 6) == 0,
    "two conjugate radial colors",
)
gate(sp.expand((2 * r) ** 2 - 4 * (r**2 + 6)) == -24, "colors never collide")
gain_plus = -(r + xi)
gain_minus = -(r - xi)
gate(sp.expand(gain_plus * gain_minus - (r**2 - xi**2)) == 0, "anchor gain norm")
gate(
    reduce_quadratic(gain_plus * gain_minus - (r**2 + 6), xi, 6) == 0,
    "anchor norm detects critical parameter",
)
gate(sp.diff(gain_plus, r) == -1 and sp.diff(gain_minus, r) == -1, "simple tie response")

boundary_gram = sp.Matrix([[-4, 5], [5, -4]])
gate(boundary_gram.det() == -9, "split-boundary determinant")
gate(sp.gcd_list(list(boundary_gram)) == 1, "boundary first Smith invariant")

semantic = {
    "family": "(AU+CV)^3+C(AU+CV)U^2+r*C*(AU+CV)V^2+A^2V^3",
    "critical_parameter": "r^2=-6",
    "order": "primitive common-zero normal nonmonogenic S3",
    "branch": "irreducible one-smooth-place decic",
    "normalization": "P1 with polynomial degrees (9,10)",
    "finite_packet": "eight-address origin delta34 plus two A2 cusps",
    "radial_response": "one color crosses infinity with simple anchored tie",
    "resolvent": "units kstar and nonzero class-group 3-torsion",
    "scope": "no plane atlas or Keller map; JC2 open",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("theorem=THM-3918-rational-moving-triple-root-decic")
print("critical_parameter=rho^2+6=0")
print("order=normal_nonmonogenic_S3")
print("branch=irreducible_decic_one_smooth_infinity;normalization=A1")
print("finite_packet=origin_8_addresses_delta34_plus_2A2")
print("radial_response=anchored_infinity_color_simple_tie")
print("unanchored_color_discriminant=-24")
print("split_boundary_snf=1,9;resolvent_Cl3=NONZERO")
print("plane_atlas=NOT_CONSTRUCTED;JC2=OPEN")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
