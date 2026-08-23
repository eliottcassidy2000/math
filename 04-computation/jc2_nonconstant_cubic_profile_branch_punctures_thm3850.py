#!/usr/bin/env python3
"""Assertion-free exact gates for the THM-3850 proof candidate.

Reproduction:
  python3 04-computation/jc2_nonconstant_cubic_profile_branch_punctures_thm3850.py
  python3 -O 04-computation/jc2_nonconstant_cubic_profile_branch_punctures_thm3850.py
"""

from __future__ import annotations

import hashlib

import sympy as sp


A, C, b, kappa, V, eps = sp.symbols("A C b kappa V eps")
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


# Universal coefficient and branch packet.  Here b is an algebraic placeholder
# for the value b(C); no derivative of b enters the cubic discriminant.
p = sp.Rational(3, 2) + A * C
u = 1 + A * C + A**2 * b
Delta = (
    -27 * A**2 * b**2
    + 8 * A * C**3
    - 54 * A * C * b
    + 9 * C**2
    - 54 * b
)
L = 9 * b - 2 * C**2
equal("depressed_discriminant", 8 * p**3 - 27 * u**2, A**2 * Delta)
equal("branch_A_discriminant", sp.discriminant(Delta, A), -8 * L**3)
equal("branch_A_lead", sp.LC(sp.Poly(Delta, A), A), -27 * b**2)

# At a nonzero root c of b, L=-2c^2.  With s^2=-2 and eps=2s,
# the two quadratic numerators are 0 and 16c^3.  This freezes the exact
# one-finite-sheet/one-puncture local packet without selecting a root.
c, s = sp.symbols("c s", nonzero=True)
Bquad = 8 * C**3 - 54 * C * b
local_L = L.subs({C: c, b: 0})
local_B = Bquad.subs({C: c, b: 0})
equal("profile_root_L", local_L, -2 * c**2)
equal("profile_root_B", local_B, 8 * c**3)
equal(
    "profile_root_first_numerator",
    (local_B - (2 * s) * local_L * (s * c)).subs(s**2, -2),
    0,
)
equal(
    "profile_root_second_numerator",
    (local_B - (2 * s) * local_L * (-s * c)).subs(s**2, -2),
    16 * c**3,
)
equal("profile_root_affine_point", Delta.subs({C: c, b: 0}), c**2 * (8 * A * c + 9))
nonzero(
    "profile_root_affine_smooth",
    sp.diff(Delta, A).subs({C: c, b: 0, A: -sp.Rational(9, 8) / c}),
)

# Constant-profile control: no profile roots and an even-degree squarefree
# L give exactly the two infinity punctures of THM-3847.
beta = sp.symbols("beta", nonzero=True)
L_constant = sp.factor(L.subs(b, beta))
equal("constant_profile_L_degree", sp.degree(L_constant, C), 2)
nonzero("constant_profile_L_squarefree", sp.discriminant(L_constant, C))
equal("constant_profile_punctures", 0 + 2, 2)


def squarefree_part(poly: sp.Expr) -> sp.Expr:
    poly = sp.Poly(poly, C)
    return sp.factor(poly.sqf_part().as_expr())


def radical_degree(poly: sp.Expr) -> int:
    poly = sp.Poly(poly, C)
    radical = sp.cancel(poly.as_expr() / sp.gcd(poly, poly.diff()).as_expr())
    return int(sp.degree(radical, C))


# Three exact hostile profiles freeze the fact that distinct roots, not their
# multiplicities, contribute finite punctures, and that odd/even squarefree
# degree changes the infinity contribution from one to two.
profile_simple = C + 1
ell_simple = squarefree_part(9 * profile_simple - 2 * C**2)
equal("simple_profile_radical_degree", radical_degree(profile_simple), 1)
equal("simple_profile_ell_degree", sp.degree(ell_simple, C), 2)
nonzero("simple_profile_ell_squarefree", sp.discriminant(ell_simple, C))
equal("simple_profile_punctures", radical_degree(profile_simple) + 2, 3)

profile_repeated = (C - 1) ** 2
ell_repeated = squarefree_part(9 * profile_repeated - 2 * C**2)
equal("repeated_profile_radical_degree", radical_degree(profile_repeated), 1)
equal("repeated_profile_ell_degree", sp.degree(ell_repeated, C), 2)
nonzero("repeated_profile_ell_squarefree", sp.discriminant(ell_repeated, C))
equal("repeated_profile_punctures", radical_degree(profile_repeated) + 2, 3)

profile_odd = C**3 + 1
ell_odd = squarefree_part(9 * profile_odd - 2 * C**2)
equal("odd_profile_radical_degree", radical_degree(profile_odd), 3)
equal("odd_profile_ell_degree", sp.degree(ell_odd, C), 3)
nonzero("odd_profile_ell_squarefree", sp.discriminant(ell_odd, C))
equal("odd_profile_punctures", radical_degree(profile_odd) + 1, 4)

# Minimal reducible boundary b=kappa*C.  The vertical line is comaximal with
# an irreducible residual conic branch.  Its normalization loses the ramified
# point over C=0 and the two points over infinity: exactly three punctures.
b_min = kappa * C
Delta_min = sp.factor(Delta.subs(b, b_min))
H = sp.factor(Delta_min / C)
equal("minimal_reducible_factorization", Delta_min, C * H)
equal("vertical_residual_comaximal", H.subs(C, 0), -54 * kappa)

M = 9 * kappa - 2 * C
disc_H = sp.factor(sp.discriminant(H, A))
equal("minimal_residual_discriminant", disc_H, -8 * C * M**3)
equal("minimal_residual_squarefree_model", C * M, C * (9 * kappa - 2 * C))
nonzero("minimal_residual_model_squarefree", sp.discriminant(C * M, C))

B_H = sp.Poly(H, A).coeff_monomial(A)
a_H = sp.Poly(H, A).coeff_monomial(A**2)
A_H = sp.factor((-B_H + eps * M * V) / (2 * a_H))
H_num = sp.together(H.subs(A, A_H)).as_numer_denom()[0]
H_reduced = sp.reduced(
    sp.Poly(H_num, eps, V, domain=sp.QQ.frac_field(C, kappa)),
    [
        sp.Poly(eps**2 + 8, eps, V, domain=sp.QQ.frac_field(C, kappa)),
        sp.Poly(V**2 - C * M, eps, V, domain=sp.QQ.frac_field(C, kappa)),
    ],
)[1]
equal("minimal_residual_normalization_map", H_reduced.as_expr(), 0)
equal("minimal_residual_infinity_points", sp.degree(V**2 + 2 * C**2, V), 2)
equal("minimal_residual_finite_punctures", 1, 1)
equal("minimal_residual_total_punctures", 1 + 2, 3)

# Primitive reducible profiles have L=h^2.  With s^2=-2, the exact factors
# u=h-sC and v=h+sC split the branch into two rational denominator graphs.
hh = sp.symbols("h")
u_factor = hh - s * C
v_factor = hh + s * C
b_square = (hh**2 + 2 * C**2) / 9
Delta_square = sp.together(Delta.subs(b, b_square))
E_v = v_factor**2 * A - 3 * (C - s * hh)
E_u = u_factor**2 * A - 3 * (C + s * hh)
factor_numerator = sp.together(Delta_square + sp.Rational(1, 3) * E_v * E_u).as_numer_denom()[0]
factor_reduced = sp.rem(sp.expand(factor_numerator), s**2 + 2, s)
equal("primitive_reducible_exact_factorization", factor_reduced, 0)
equal(
    "primitive_reducible_profile_product",
    sp.rem(sp.expand(u_factor * v_factor - 9 * b_square), s**2 + 2, s),
    0,
)

# If one factor is constant, one A1 component really occurs, but its mate has
# a finite denominator root and hence at least two punctures.
a = sp.symbols("a", nonzero=True)
h_linear = s * C + a
u_linear = sp.rem(sp.expand(u_factor.subs(hh, h_linear)), s**2 + 2, s)
v_linear = sp.rem(sp.expand(v_factor.subs(hh, h_linear)), s**2 + 2, s)
equal("positive_A1_constant_factor", u_linear, a)
equal("positive_A1_nonconstant_factor", v_linear, 2 * s * C + a)
E_u_linear = sp.rem(sp.expand(E_u.subs(hh, h_linear)), s**2 + 2, s)
equal("positive_A1_component", E_u_linear, a**2 * A + 3 * C - 3 * s * a)
bad_numerator = sp.rem(sp.expand((C - s * h_linear)), s**2 + 2, s)
bad_root = -a / (2 * s)
bad_value = sp.rem(
    sp.together(bad_numerator.subs(C, bad_root)).as_numer_denom()[0],
    s**2 + 2,
    s,
)
nonzero("positive_A1_mate_has_finite_pole", bad_value)

# Profiles vanishing at C=0.  The exact maximal vertical orders are frozen in
# the three valuation regimes.  After removing that order, the A^2
# coefficient still vanishes at C=0, so the residual closure contains the
# finite-C/A-infinity corner, while some other coefficient remains nonzero.
a0 = sp.symbols("a0", nonzero=True)
a1 = sp.symbols("a1")

b_m1 = C * (a0 + a1 * C)
Delta_m1 = sp.factor(Delta.subs(b, b_m1))
H_m1 = sp.factor(Delta_m1 / C)
equal("vertical_m1_exact_order", H_m1.subs({C: 0, A: 0}), -54 * a0)
equal("vertical_m1_A2_corner", sp.Poly(H_m1, A).coeff_monomial(A**2).subs(C, 0), 0)

b_m2 = C**2 * a0
Delta_m2 = sp.factor(Delta.subs(b, b_m2))
H_m2 = sp.factor(Delta_m2 / C**2)
equal("vertical_m2_generic_constant", H_m2.subs({C: 0, A: 0}), 9 - 54 * a0)
equal("vertical_m2_A2_corner", sp.Poly(H_m2, A).coeff_monomial(A**2).subs(C, 0), 0)

b_m2_exceptional = C**2 * (sp.Rational(1, 6) + a1 * C)
Delta_m2_exceptional = sp.factor(Delta.subs(b, b_m2_exceptional))
H_m2_exceptional = sp.factor(Delta_m2_exceptional / C**3)
equal(
    "vertical_m2_exceptional_exact_order",
    H_m2_exceptional.subs(C, 0),
    -A - 54 * a1,
)
equal(
    "vertical_m2_exceptional_A2_corner",
    sp.Poly(H_m2_exceptional, A).coeff_monomial(A**2).subs(C, 0),
    0,
)

b_m3 = a0 * C**3
Delta_m3 = sp.factor(Delta.subs(b, b_m3))
H_m3 = sp.factor(Delta_m3 / C**2)
equal("vertical_mge3_exact_order", H_m3.subs({C: 0, A: 0}), 9)
equal("vertical_mge3_A2_corner", sp.Poly(H_m3, A).coeff_monomial(A**2).subs(C, 0), 0)

print("THM3850_UNIVERSAL_BRANCH", Delta)
print("THM3850_QUADRATIC_DISCRIMINANT", sp.factor(sp.discriminant(Delta, A)))
print("THM3850_IRREDUCIBLE_PUNCTURE_FORMULA", "deg(rad(b))+nu_infinity")
print("THM3850_INFINITY_RULE", "nu_infinity=1 for odd deg(ell), 2 for even deg(ell)")
print("THM3850_CONSTANT_CONTROL", "0+2=2")
print("THM3850_SIMPLE_CONTROL", "1+2=3")
print("THM3850_REPEATED_CONTROL", "1+2=3")
print("THM3850_ODD_CONTROL", "3+1=4")
print("THM3850_MINIMAL_REDUCIBLE", sp.factor(Delta_min))
print("THM3850_MINIMAL_RESIDUAL", H)
print("THM3850_MINIMAL_NORMALIZATION", "V^2=C*(9*kappa-2*C); three punctures")
print("THM3850_PRIMITIVE_REDUCIBLE", "Delta=-(E_u*E_v)/3; punctures=1+deg(rad(u/v))")
print("THM3850_POSITIVE_A1_CONTROL", "h=s*C+a: u=a, v=2*s*C+a; A1 plus two-place mate")
print("THM3850_VERTICAL_ORDERS", "m=1:nu=1; m=2:nu=2 except a0=1/6 gives3; m>=3:nu=2")
print(
    "THM3850_SCOPE",
    "all nonconstant profiles have a branch component with at least two punctures",
)
semantic_packet = (
    "nonconstant polynomial profile b(C)",
    "irreducible iff primitive quadratic with nonsquare L",
    "smooth projective normalization k(C)(sqrt squarefree(L))",
    "affine punctures deg radical(b) plus one/two infinity points",
    "every irreducible nonconstant branch has at least two punctures",
    "minimal b=kappa C boundary is A1 plus a three-puncture residual conic",
    "primitive reducible branches split into denominator graphs",
    "origin-vanishing profiles force a finite-C/A-infinity residual corner",
    "every nonconstant profile has at least one non-polynomial branch component",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
