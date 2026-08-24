#!/usr/bin/env python3
"""Exact companion for THM-3920's depressed-cubic chart obstruction.

Reproduction:
  python3 04-computation/jc2_affine_plane_boundary_unibranch_depressed_cubic_chart_thm3920.py
  python3 -O 04-computation/jc2_affine_plane_boundary_unibranch_depressed_cubic_chart_thm3920.py
"""

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


s, z, x, t, w, b = sp.symbols("s z x t w b")
aa, cc = sp.symbols("a c")
P, H = sp.symbols("P H")


# ---------------------------------------------------------------------------
# Universal depressed-cubic and rational-chart identities.
# ---------------------------------------------------------------------------

F_symbolic = z**3 - 3 * P * z + 2 * H
r_symbolic = P**3 - H**2
zero(
    sp.discriminant(F_symbolic, z) - 108 * r_symbolic,
    "depressed-cubic discriminant",
)
zero(
    sp.resultant(F_symbolic, z**2 - P, z) + 4 * r_symbolic,
    "derivative-curve intersection resultant",
)
zero(
    sp.rem(F_symbolic, z**2 - P, z) - 2 * (H - P * z),
    "restriction to the derivative curve",
)

p_fun = sp.Function("p")(s)
h_fun = sp.Function("h")(s)
F_fun = z**3 - 3 * p_fun * z + 2 * h_fun
A_fun = F_fun / 4
C_fun = s * A_fun
chart_jacobian = sp.factor(
    sp.diff(A_fun, s) * sp.diff(C_fun, z)
    - sp.diff(A_fun, z) * sp.diff(C_fun, s)
)
zero(chart_jacobian + F_fun * sp.diff(F_fun, z) / 16, "universal chart Jacobian")


# ---------------------------------------------------------------------------
# The constant-square Chebyshev seam and its third-Veronese normalization.
# ---------------------------------------------------------------------------

F_cheb = z**3 - 3 * z + 2 * s
A_cheb = F_cheb / 4
C_cheb = s * A_cheb
b_cheb = sp.factor(4 * A_cheb**2 - 2 * C_cheb)
zero(b_cheb - A_cheb * (z**3 - 3 * z), "Chebyshev triangular target coordinate")
zero(F_cheb.subs(z, 1) - 2 * (s - 1), "positive derivative line attachment")
zero(F_cheb.subs(z, -1) - 2 * (s + 1), "negative derivative line attachment")

b_cover = sp.cancel((aa * (z**3 - 3 * z)).subs({aa: t**3, z: w / t}))
zero(b_cover - (w**3 - 3 * w * t**2), "Chebyshev mu3 cover")
cover_jacobian = sp.factor(
    sp.diff(t**3, t) * sp.diff(b_cover, w)
    - sp.diff(t**3, w) * sp.diff(b_cover, t)
)
zero(cover_jacobian - 9 * t**2 * (w**2 - t**2), "three ramification rays")

veronese = (t**3, t**2 * w, t * w**2, w**3)
zero(b_cover - (veronese[3] - 3 * veronese[1]), "target lies in third Veronese")
zero(veronese[1] ** 2 - veronese[0] * veronese[2], "Veronese relation 1")
zero(veronese[1] * veronese[2] - veronese[0] * veronese[3], "Veronese relation 2")
zero(veronese[2] ** 2 - veronese[1] * veronese[3], "Veronese relation 3")


# ---------------------------------------------------------------------------
# Nonreduced p=0 seam: exact Kummer valuation controls.
# ---------------------------------------------------------------------------

for exponent in (1, 2, 4, 5, 7, 8):
    numerator = 4 * aa ** (exponent + 1) - cc**exponent
    gate(sp.Poly(numerator, aa).coeff_monomial(aa**0) == -cc**exponent,
         f"Kummer a-line numerator is a unit for m={exponent}")
    gate(exponent % 3 != 0, f"Kummer exponent is nontrivial mod three for m={exponent}")


# ---------------------------------------------------------------------------
# Exact symmetric polynomial-part family used as a quantitative control.
# ---------------------------------------------------------------------------

f1 = 4 * b + 1
f2 = 8 * b**2 + 2 * b + 1
A_factor = 48 * b**2 + 16 * b + 3
B_factor = 64 * b**3 + 48 * b**2 + 24 * b + 5
K_factor = 2304 * b**5 + 10176 * b**4 + 4064 * b**3 + 996 * b**2 + 84 * b + 5

exceptional = (b, f1, f2, A_factor, B_factor, K_factor)
for index, left in enumerate(exceptional):
    gate(sp.gcd(left, sp.diff(left, b)) == 1,
         f"exceptional factor {index} is squarefree")
    for right_index in range(index):
        gate(sp.gcd(left, exceptional[right_index]) == 1,
             f"exceptional factors {right_index},{index} are disjoint")

p = (s**2 - 1) * (s**2 + b) ** 2
h = (
    s**9
    + (3 * b - sp.Rational(3, 2)) * s**7
    + (3 * b**2 - sp.Rational(9, 2) * b + sp.Rational(3, 8)) * s**5
    + (b**3 - sp.Rational(9, 2) * b**2 + sp.Rational(9, 8) * b + sp.Rational(1, 16)) * s**3
    + (-sp.Rational(3, 2) * b**3 + sp.Rational(9, 8) * b**2 + sp.Rational(3, 16) * b + sp.Rational(3, 128)) * s
)

scaled_p = sp.expand(t**6 * p.subs(s, 1 / t))
scaled_h = sp.expand(t**9 * h.subs(s, 1 / t))
polynomial_part = sp.series(
    (1 + b * t**2) ** 3 * (1 - t**2) ** sp.Rational(3, 2),
    t,
    0,
    10,
).removeO()
zero(scaled_p - (1 - t**2) * (1 + b * t**2) ** 2, "scaled symmetric p")
zero(scaled_h - polynomial_part, "symmetric h is the polynomial part")

r = sp.factor(p**3 - h**2)
Sx = (
    384 * f1 * f2 * x**4
    + 32 * (1152 * b**4 + 64 * b**3 - 36 * b - 11) * x**3
    + 48 * (768 * b**5 - 640 * b**4 - 240 * b**3 - 120 * b**2 - 26 * b - 1) * x**2
    + 3 * (4096 * b**6 - 14336 * b**5 - 3840 * b**4 - 1920 * b**3 - 480 * b**2 - 48 * b - 3) * x
    - 16384 * b**6
)
zero(16384 * r - Sx.subs(x, s**2), "even residual quartic")
zero(Sx.subs(x, 0) + 16384 * b**6, "origin address factor")
zero(Sx.subs(x, -b) - b * A_factor**2, "double-factor address")
zero(Sx.subs(x, 1) + B_factor**2, "simple-factor address")
zero(
    sp.discriminant(Sx, x) + 27 * 2**14 * A_factor**6 * B_factor**3 * K_factor,
    "residual quartic discriminant",
)
zero(
    sp.discriminant(Sx.subs(x, s**2), s)
    + 2**57 * 3**7 * b**6 * f1 * f2 * A_factor**12 * B_factor**6 * K_factor**2,
    "full finite-address discriminant",
)


def reduce_mod(expression: sp.Expr, modulus: sp.Expr) -> sp.Expr:
    """Reduce a rational expression in b modulo a squarefree Q-polynomial."""

    modulus_poly = sp.Poly(modulus, b, domain=sp.QQ)
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    numerator_remainder = sp.rem(sp.Poly(numerator, b, domain=sp.QQ), modulus_poly).as_expr()
    denominator_remainder = sp.rem(sp.Poly(denominator, b, domain=sp.QQ), modulus_poly).as_expr()
    denominator_inverse = sp.invert(sp.Poly(denominator_remainder, b, domain=sp.QQ), modulus_poly).as_expr()
    return sp.factor(
        sp.rem(
            sp.Poly(numerator_remainder * denominator_inverse, b, domain=sp.QQ),
            modulus_poly,
        ).as_expr()
    )


def reduce_poly_x(expression: sp.Expr, modulus: sp.Expr) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), x)
    return sp.factor(
        sum(
            reduce_mod(polynomial.coeff_monomial(x**power), modulus) * x**power
            for power in range(polynomial.degree() + 1)
        )
    )


# b=0: one bad x=0 address and three simple good x-roots.
S_b0 = sp.factor(Sx.subs(b, 0))
zero(S_b0 - x * (384 * x**3 - 352 * x**2 - 48 * x - 9), "b=0 factorization")
C_b0 = sp.cancel(S_b0 / x)
gate(sp.discriminant(C_b0, x) == -2488320000, "b=0 cubic is squarefree")
gate(C_b0.subs(x, 0) == -9 and C_b0.subs(x, 1) == -25, "b=0 good roots")

# f2=0: the root at infinity leaves a squarefree good cubic.
S_f2 = reduce_poly_x(Sx, f2)
expected_f2 = 448 * b * x**3 - (672 * b + 144) * x**2 + (144 * b + 159) * x + 80 * b - 8
zero(S_f2 - expected_f2, "f2 cubic reduction")
zero(reduce_mod(sp.discriminant(S_f2, x), f2) - 8297856 * (20 * b - 1), "f2 cubic discriminant")
zero(reduce_mod(S_f2.subs(x, 0), f2) - 8 * (10 * b - 1), "f2 avoids zero")
zero(reduce_mod(S_f2.subs(x, 1), f2) - 7, "f2 avoids x=1")
zero(reduce_mod(S_f2.subs(x, -b), f2) - sp.Rational(7, 2) * f1, "f2 avoids x=-b")

# A=0: an exact double bad root x=-b and two simple good roots.
QA = 144 * (10 * b + 3) * x**2 - (1392 * b + 531) * x + 84 - 128 * b
zero(reduce_poly_x(Sx - sp.Rational(16, 27) * (x + b) ** 2 * QA, A_factor), "A-factor quotient")
zero(reduce_mod(sp.discriminant(QA, x), A_factor) - 10125 * (32 * b - 3), "A quotient discriminant")
zero(reduce_mod(QA.subs(x, 0), A_factor) - 4 * (21 - 32 * b), "A quotient avoids zero")
zero(reduce_mod(QA.subs(x, 1), A_factor) + 5 * (16 * b + 3), "A quotient avoids x=1")
zero(reduce_mod(QA.subs(x, -b), A_factor) + 135 * b, "A double root is exact")

# B=0: an exact double bad root x=1 and two simple good roots.
QB = (
    (192 * b**2 + 144 * b + 36) * x**2
    + (48 * b**2 - 60 * b - 31) * x
    - (60 * b**2 + 39 * b + 5)
)
zero(reduce_poly_x(Sx + 16 * (x - 1) ** 2 * QB, B_factor), "B-factor quotient")
zero(reduce_mod(sp.discriminant(QB, x), B_factor) - 2 * (288 * b**2 + 6 * b - 37), "B quotient discriminant")
zero(reduce_mod(QB.subs(x, 0), B_factor) + 60 * b**2 + 39 * b + 5, "B quotient avoids zero")
zero(reduce_mod(QB.subs(x, 1), B_factor) - 45 * b * f1, "B double root is exact")
zero(reduce_mod(4 * QB.subs(x, -b), B_factor) + 5 * f1, "B quotient avoids x=-b")

# K=0: one exact double good root and a disjoint squarefree good quadratic.
xi = (
    95616 * b**4 + 421728 * b**3 + 162368 * b**2 + 16726 * b + 965
) / sp.Integer(15120)
q2 = 384 * f1 * f2
q1 = sp.Rational(32, 5) * (4224 * b**4 + 1472 * b**3 + 192 * b**2 - 156 * b - 55)
q0 = -sp.Rational(16, 45) * (188352 * b**4 + 111696 * b**3 + 34556 * b**2 + 5272 * b + 215)
QK = q2 * x**2 + q1 * x + q0
factor_difference = sp.Poly(sp.expand(Sx - (x - xi) ** 2 * QK), x)
for power in range(5):
    zero(reduce_mod(factor_difference.coeff_monomial(x**power), K_factor), f"K factor coefficient {power}")

K_nonzero = (
    sp.cancel(xi).as_numer_denom()[0],
    sp.cancel(xi - 1).as_numer_denom()[0],
    sp.cancel(xi + b).as_numer_denom()[0],
    sp.cancel(QK.subs(x, 0)).as_numer_denom()[0],
    sp.cancel(QK.subs(x, 1)).as_numer_denom()[0],
    sp.cancel(QK.subs(x, -b)).as_numer_denom()[0],
    sp.cancel(QK.subs(x, xi)).as_numer_denom()[0],
    sp.cancel(sp.discriminant(QK, x)).as_numer_denom()[0],
)
for index, residue in enumerate(K_nonzero):
    gate(sp.gcd(K_factor, residue) == 1, f"K noncollision residue {index}")

# The collapsed cubic is irreducible exactly off f1=0.
F_family = z**3 - 3 * p * z + 2 * h
c3, c2, c1, c0 = sp.symbols("c3 c2 c1 c0")
root = c3 * s**3 + c2 * s**2 + c1 * s + c0
root_equation = sp.Poly(sp.expand(F_family.subs(z, root)), s)


def root_coefficient(power: int, substitutions: dict[sp.Symbol, sp.Expr]) -> sp.Expr:
    return sp.factor(root_equation.coeff_monomial(s**power).subs(substitutions))


zero(root_coefficient(9, {}) - (c3 - 1) ** 2 * (c3 + 2), "root leading choices")
minus = {c3: -2}
zero(root_coefficient(8, minus) - 9 * c2, "minus c2")
minus[c2] = 0
zero(root_coefficient(7, minus) - 9 * (2 * b + c1 - 1), "minus c1")
minus[c1] = 1 - 2 * b
zero(root_coefficient(6, minus) - 9 * c0, "minus c0")
minus[c0] = 0
zero(root_coefficient(5, minus) + sp.Rational(9, 4) * f1, "minus boundary")
gate(root_coefficient(1, minus).subs(b, -sp.Rational(1, 4)) == sp.Rational(27, 64), "minus terminal contradiction")

plus = {c3: 1}
zero(root_coefficient(7, plus) - 3 * c2**2, "plus c2")
plus[c2] = 0
zero(root_coefficient(5, plus) - sp.Rational(3, 4) * (2 * b - 2 * c1 - 1) ** 2, "plus c1")
plus[c1] = b - sp.Rational(1, 2)
zero(root_coefficient(3, plus) - 3 * c0**2, "plus c0")
plus[c0] = 0
zero(root_coefficient(1, plus) - sp.Rational(3, 64) * f1**2, "plus reducible boundary")

boundary_factorization = sp.factor(F_family.subs(b, -sp.Rational(1, 4)))
boundary_root = s**3 - sp.Rational(3, 4) * s
zero(F_family.subs({b: -sp.Rational(1, 4), z: boundary_root}), "reducible boundary root")
gate(boundary_factorization != F_family.subs(b, -sp.Rational(1, 4)), "reducible boundary factors")


semantic_payload = {
    "boundary": "normal_affine_completion_of_A2_has_rational_unibranch_boundary_curves",
    "address_cap": "finite_flat_degree_d_boundary_normalization_fibre_has_at_most_d_points",
    "nonsquare_case": "support_residual_at_most_one_then_C3_character_forces_p_zero",
    "nonconstant_square_case": "two_pure_power_attachments_Mason_or_reducible",
    "constant_nonzero_square_case": "Chebyshev_third_Veronese_deleted_total_ray_makes_a_unit",
    "zero_square_case": "pure_Kummer_C3_excluded_by_THM3801_with_low_exponent_valuation_controls",
    "symmetric_family": "generic8_b_f2_K6_A_B4_good_addresses_f1_reducible",
    "scope": "full_irreducible_depressed_cubic_radial_chart_closed_JC2_open",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("theorem=THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction")
print("boundary=irreducible_curves_rational_and_unibranch")
print("address_cap=normalization_addresses_over_y_at_most_finite_flat_degree")
print("nonsquare_derivative=empty_by_support_one_C3_character_contradiction")
print("square_nonconstant=empty_by_two_attachments_Mason_and_same_center_reducibility")
print("square_constant_nonzero=empty_by_Chebyshev_third_Veronese_unit")
print("square_zero=empty_by_pure_Kummer_C3_THM3801;low_exponent_valuation_controls")
print("symmetric_family=generic8;b0_f2_K6;A_B4;f1_reducible")
print("scope=full_irreducible_depressed_cubic_radial_chart_closed;JC2_open")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
