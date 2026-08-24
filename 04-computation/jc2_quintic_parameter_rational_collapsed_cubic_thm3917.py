#!/usr/bin/env python3
"""Exact companion for THM-3917's genus-zero collapsed cubic.

Reproduction:
  python3 04-computation/jc2_quintic_parameter_rational_collapsed_cubic_thm3917.py
  python3 -O 04-computation/jc2_quintic_parameter_rational_collapsed_cubic_thm3917.py
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


u, x, z, t, y, beta = sp.symbols("u x z t y beta")
AA, CC, ZZ = sp.symbols("A C Z")
c3, c2, c1, c0 = sp.symbols("c3 c2 c1 c0")

K = (
    2304 * beta**5
    + 10176 * beta**4
    + 4064 * beta**3
    + 996 * beta**2
    + 84 * beta
    + 5
)
K_poly = sp.Poly(K, beta, domain=sp.QQ)


def reduce_mod_K(expression: sp.Expr) -> sp.Expr:
    """Return the degree-<5 representative in Q[beta]/(K)."""

    numerator, denominator = sp.cancel(expression).as_numer_denom()
    numerator_remainder = sp.rem(
        sp.Poly(numerator, beta, domain=sp.QQ), K_poly
    ).as_expr()
    denominator_remainder = sp.rem(
        sp.Poly(denominator, beta, domain=sp.QQ), K_poly
    ).as_expr()
    denominator_inverse = sp.invert(
        sp.Poly(denominator_remainder, beta, domain=sp.QQ), K_poly
    ).as_expr()
    return sp.rem(
        sp.Poly(numerator_remainder * denominator_inverse, beta, domain=sp.QQ),
        K_poly,
    ).as_expr()


# ---------------------------------------------------------------------------
# The quintic parameter and the polynomial-part construction.
# ---------------------------------------------------------------------------

gate(sp.Poly(K, beta, modulus=11).degree() == 5, "quintic keeps degree mod 11")
gate(sp.Poly(K, beta, modulus=11).is_irreducible, "quintic is irreducible mod 11")
gate(
    sp.discriminant(K, beta) == 1170238183833600000000000,
    "quintic is squarefree",
)

exceptional_factors = (
    beta,
    4 * beta + 1,
    8 * beta**2 + 2 * beta + 1,
    48 * beta**2 + 16 * beta + 3,
    64 * beta**3 + 48 * beta**2 + 24 * beta + 5,
)
for factor in exceptional_factors:
    gate(sp.gcd(K, factor) == 1, f"quintic avoids {factor}")

p = (u**2 - 1) * (u**2 + beta) ** 2
h = (
    u**9
    + (3 * beta - sp.Rational(3, 2)) * u**7
    + (3 * beta**2 - sp.Rational(9, 2) * beta + sp.Rational(3, 8)) * u**5
    + (
        beta**3
        - sp.Rational(9, 2) * beta**2
        + sp.Rational(9, 8) * beta
        + sp.Rational(1, 16)
    )
    * u**3
    + (
        -sp.Rational(3, 2) * beta**3
        + sp.Rational(9, 8) * beta**2
        + sp.Rational(3, 16) * beta
        + sp.Rational(3, 128)
    )
    * u
)

P_t = sp.expand(t**6 * p.subs(u, 1 / t))
H_t = sp.expand(t**9 * h.subs(u, 1 / t))
expected_polynomial_part = sp.series(
    (1 + beta * t**2) ** 3 * (1 - t**2) ** sp.Rational(3, 2),
    t,
    0,
    10,
).removeO()
zero(H_t - expected_polynomial_part, "h is the polynomial part of p^(3/2)")
zero(P_t - (1 - t**2) * (1 + beta * t**2) ** 2, "scaled p at infinity")

r = sp.factor(p**3 - h**2)
S = (
    384 * (4 * beta + 1) * (8 * beta**2 + 2 * beta + 1) * x**4
    + 32 * (1152 * beta**4 + 64 * beta**3 - 36 * beta - 11) * x**3
    + 48
    * (768 * beta**5 - 640 * beta**4 - 240 * beta**3 - 120 * beta**2 - 26 * beta - 1)
    * x**2
    + 3
    * (
        4096 * beta**6
        - 14336 * beta**5
        - 3840 * beta**4
        - 1920 * beta**3
        - 480 * beta**2
        - 48 * beta
        - 3
    )
    * x
    - 16384 * beta**6
)
zero(16384 * r - S.subs(x, u**2), "even residual quartic identity")
gate(sp.degree(r, u) == 8, "generic residual degree eight")

A_factor = 48 * beta**2 + 16 * beta + 3
B_factor = 64 * beta**3 + 48 * beta**2 + 24 * beta + 5
zero(
    sp.discriminant(S, x) + 27 * 2**14 * A_factor**6 * B_factor**3 * K,
    "quartic discriminant factorization",
)
zero(
    sp.resultant(p, h, u)
    + beta**2 * A_factor**4 * B_factor**2 / 2**42,
    "p and h have no common zero on the quintic locus",
)

# ---------------------------------------------------------------------------
# The exact repeated root and the four finite ramification points.
# ---------------------------------------------------------------------------

x0 = (
    95616 * beta**4
    + 421728 * beta**3
    + 162368 * beta**2
    + 16726 * beta
    + 965
) / sp.Integer(15120)
q2 = 384 * (4 * beta + 1) * (8 * beta**2 + 2 * beta + 1)
q1 = sp.Rational(32, 5) * (
    4224 * beta**4 + 1472 * beta**3 + 192 * beta**2 - 156 * beta - 55
)
q0 = -sp.Rational(16, 45) * (
    188352 * beta**4
    + 111696 * beta**3
    + 34556 * beta**2
    + 5272 * beta
    + 215
)
Q = q2 * x**2 + q1 * x + q0

for power in range(5):
    coefficient = sp.Poly(sp.expand(S - (x - x0) ** 2 * Q), x).coeff_monomial(
        x**power
    )
    zero(reduce_mod_K(coefficient), f"K-quotient factorization coefficient {power}")

E_factor = 35712 * beta**4 + 41952 * beta**3 + 13408 * beta**2 + 2438 * beta + 85
D_factor = (
    114974016 * beta**4
    + 47144272 * beta**3
    + 11949168 * beta**2
    + 978444 * beta
    + 58837
)
zero(reduce_mod_K(Q.subs(x, x0)) + sp.Rational(2, 3) * E_factor, "double root is not a root of Q")
zero(
    reduce_mod_K(sp.discriminant(Q, x)) + sp.Rational(1024, 27) * D_factor,
    "remaining quadratic discriminant",
)

nonzero_residues = (
    sp.cancel(x0).as_numer_denom()[0],
    E_factor,
    D_factor,
    sp.cancel(q0).as_numer_denom()[0],
)
for residue in nonzero_residues:
    gate(sp.gcd(K, residue) == 1, "required K-residue is nonzero")

zero(S.subs(x, 1) + B_factor**2, "residual roots avoid u^2=1")
zero(S.subs(x, -beta) - beta * A_factor**2, "residual roots avoid u^2=-beta")

subresultants = sp.subresultants(S, sp.diff(S, x), x)
gate(sp.degree(subresultants[-2], x) == 1, "penultimate subresultant is linear")
gate(
    any(
        reduce_mod_K(sp.Poly(subresultants[-2], x).coeff_monomial(x**power)) != 0
        for power in (0, 1)
    ),
    "linear subresultant survives modulo K",
)
gate(reduce_mod_K(subresultants[-1]) == 0, "resultant vanishes modulo K")

# ---------------------------------------------------------------------------
# Irreducibility of the collapsed cubic.
# ---------------------------------------------------------------------------

F = z**3 - 3 * p * z + 2 * h
zero(sp.discriminant(F, z) - 108 * r, "collapsed-cubic discriminant")
zero(sp.resultant(F, sp.diff(F, z) / 3, z) + 4 * r, "F--F_z intersection eliminant")

root_ansatz = c3 * u**3 + c2 * u**2 + c1 * u + c0
root_equation = sp.Poly(sp.expand(F.subs(z, root_ansatz)), u)

def root_coefficient(power: int, substitutions: dict[sp.Symbol, sp.Expr]) -> sp.Expr:
    return sp.factor(root_equation.coeff_monomial(u**power).subs(substitutions))


zero(root_coefficient(9, {}) - (c3 - 1) ** 2 * (c3 + 2), "root leading alternatives")

minus_chain = {c3: -2}
zero(root_coefficient(8, minus_chain) - 9 * c2, "minus branch c2 gate")
minus_chain[c2] = 0
zero(root_coefficient(7, minus_chain) - 9 * (2 * beta + c1 - 1), "minus branch c1 gate")
minus_chain[c1] = 1 - 2 * beta
zero(root_coefficient(6, minus_chain) - 9 * c0, "minus branch c0 gate")
minus_chain[c0] = 0
zero(root_coefficient(5, minus_chain) + sp.Rational(9, 4) * (4 * beta + 1), "minus branch parameter gate")
gate(
    root_coefficient(1, minus_chain).subs(beta, -sp.Rational(1, 4))
    == sp.Rational(27, 64),
    "minus branch terminal contradiction",
)

plus_chain = {c3: 1}
zero(root_coefficient(7, plus_chain) - 3 * c2**2, "plus branch c2 gate")
plus_chain[c2] = 0
zero(
    root_coefficient(5, plus_chain)
    - sp.Rational(3, 4) * (2 * beta - 2 * c1 - 1) ** 2,
    "plus branch c1 gate",
)
plus_chain[c1] = beta - sp.Rational(1, 2)
zero(root_coefficient(3, plus_chain) - 3 * c0**2, "plus branch c0 gate")
plus_chain[c0] = 0
zero(
    root_coefficient(1, plus_chain) - sp.Rational(3, 64) * (4 * beta + 1) ** 2,
    "plus branch reducible boundary",
)
gate(K.subs(beta, -sp.Rational(1, 4)) == sp.Rational(81, 4), "quintic avoids reducible boundary")

# ---------------------------------------------------------------------------
# Three unramified infinity places and the rational ramification curve.
# ---------------------------------------------------------------------------

D_infinity = (4 * beta + 1) * (8 * beta**2 + 2 * beta + 1)
s_formal = (1 + beta * t**2) * sp.sqrt(1 - t**2)
E_formal = sp.series(H_t - s_formal**3, t, 0, 13).removeO().expand()
gate(
    sp.factor(E_formal.coeff(t, 10)) == -sp.Rational(3, 256) * D_infinity,
    "first split-sheet error at order ten",
)
gate(sp.gcd(K, D_infinity) == 1, "split-sheet leading coefficient nonzero")

# For y=s+t^5 q, the residual equation modulo t is
# 3(q^2-D/128); for y=-2s+t^10 q it is 9q-3D/128.
q = sp.symbols("q")
double_residual = 3 * (q**2 - D_infinity / 128)
simple_residual = 9 * q - 3 * D_infinity / 128
zero(sp.discriminant(double_residual, q) - sp.Rational(9, 32) * D_infinity, "two split infinity roots")
gate(sp.diff(simple_residual, q) == 9, "simple infinity root")
zero(simple_residual.subs(q, D_infinity / 384), "simple infinity correction")

# F_z=0 normalizes to the conic v^2=u^2-1 after z=(u^2+beta)v.
v = sp.symbols("v")
zero(
    (sp.diff(F, z) / 3).subs(z, (u**2 + beta) * v)
    - (u**2 + beta) ** 2 * (v**2 - u**2 + 1),
    "rational ramification-curve normalization",
)

# The universal rational chart and its two coefficient-one Jacobian divisors.
A_chart = F / 4
C_chart = u * F / 4
chart_jacobian = sp.factor(
    sp.diff(A_chart, u) * sp.diff(C_chart, z)
    - sp.diff(A_chart, z) * sp.diff(C_chart, u)
)
zero(chart_jacobian + F * sp.diff(F, z) / 16, "rational-chart Jacobian")

# The associated integral cubic over the target plane.  Its A^10 factor is
# power-basis index debt, while Delta is the reduced field-branch candidate.
P_global = sp.cancel(AA**6 * p.subs(u, CC / AA))
H_global = sp.cancel(AA**9 * h.subs(u, CC / AA))
Q_global = 2 * AA**10 - H_global
f_global = ZZ**3 - 3 * P_global * ZZ - 2 * Q_global
Delta_global = sp.cancel(AA**8 * r.subs(u, CC / AA) + 4 * H_global - 4 * AA**10)
gate(P_global.is_polynomial(AA, CC), "global P is polynomial")
gate(H_global.is_polynomial(AA, CC), "global H is polynomial")
gate(Delta_global.is_polynomial(AA, CC), "global branch candidate is polynomial")
zero(
    f_global.subs({CC: u * AA, ZZ: AA**3 * z}) - AA**9 * (F - 4 * AA),
    "global cubic rational chart",
)
zero(
    sp.discriminant(f_global, ZZ) - 108 * AA**10 * Delta_global,
    "global power-basis discriminant",
)

semantic_payload = {
    "parameter": "K(b)=2304b5+10176b4+4064b3+996b2+84b+5_irreducible",
    "collapsed_divisor": "F=z3-3pz+2h_irreducible_genus0",
    "finite_packet": "four_simple_ramification_plus_two_unramified_nodes",
    "infinity_packet": "three_unramified_places",
    "ramification_divisor": "rational_with_six_normalization_branches_at_one_X_point",
    "general_gate": "irreducible_boundary_curve_of_A2_open_in_normal_surface_is_unibranch",
    "conclusion": "same_field_Keller_atlas_empty",
    "scope": "specific_quintic_parameter_design_closed_JC2_open",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("theorem=THM-3917-quintic-parameter-rational-collapsed-cubic")
print("parameter_quintic=2304*b^5+10176*b^4+4064*b^3+996*b^2+84*b+5")
print("quintic=irreducible_mod_11;squarefree")
print("residual=one_double_x_root_plus_two_simple_x_roots")
print("finite_discriminant=two_double_u_roots_plus_four_simple_u_roots")
print("finite_ramification=4;double_roots=two_unramified_nodes")
print("infinity_places=3;infinity_ramification=0")
print("collapsed_cubic=irreducible;normalization_genus=0;K_rational_infinity_point=1")
print("ramification_curve=rational;collapsed_intersection_support=6")
print("chart_jacobian=-F*F_z/16")
print("normalization_ramification_image=at_least_6_branches_at_one_point")
print("scope=genus_obstruction_escaped_but_boundary_unibranch_kills_same_field_Keller_atlas;JC2_open")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
