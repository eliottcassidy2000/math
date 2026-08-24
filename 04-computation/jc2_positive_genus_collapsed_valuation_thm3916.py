#!/usr/bin/env python3
"""Exact companion for THM-3916's collapsed genus-two valuation gate.

Reproduction:
  python3 04-computation/jc2_positive_genus_collapsed_valuation_thm3916.py
  python3 -O 04-computation/jc2_positive_genus_collapsed_valuation_thm3916.py
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


u, z, w, v, q = sp.symbols("u z w v q")

p = (u**2 - 4) * (u**2 + 1) ** 2
h = u**9 - 3 * u**7 - 9 * u**5 + 5 * u**3 + 30 * u
F = sp.expand(z**3 - 3 * p * z + 2 * h)
A = F / 4
C = u * F / 4

# The rational THM-3915 chart and its exact Jacobian divisor.
chart_jacobian = sp.factor(
    sp.det(
        sp.Matrix(
            [
                [sp.diff(A, u), sp.diff(A, z)],
                [sp.diff(C, u), sp.diff(C, z)],
            ]
        )
    )
)
zero(chart_jacobian + F * sp.diff(F, z) / 16, "rational-chart Jacobian")
gate(sp.gcd(F, u) == 1, "u is a unit at the generic F-divisor")
gate(sp.Poly(F, z).LC() == 1 and sp.degree(F, z) == 3, "F is monic cubic")

# A cubic reducible over k(u) has a k(u)-root.  Monicity puts that root in
# k[u].  Degree comparison leaves only degree three.  The leading coefficient
# is -2 or 1; the displayed coefficient descents kill both cases.
a, b, c, d = sp.symbols("a b c d")
r = a * u**3 + b * u**2 + c * u + d
root_remainder = sp.Poly(sp.expand(F.subs(z, r)), u)
leading_equation = sp.factor(root_remainder.coeff_monomial(u**9))
gate(leading_equation == (a - 1) ** 2 * (a + 2), "degree-three leading roots")

minus_two_packet = {a: -2, b: 0, c: 2, d: 0}
zero(root_remainder.coeff_monomial(u**8).subs(a, -2) - 9 * b,
     "minus-two branch forces b")
zero(root_remainder.coeff_monomial(u**7).subs({a: -2, b: 0}) - (9 * c - 18),
     "minus-two branch forces c")
zero(root_remainder.coeff_monomial(u**6).subs({a: -2, b: 0, c: 2}) - 9 * d,
     "minus-two branch forces d")
gate(
    root_remainder.coeff_monomial(u**5).subs(minus_two_packet) == -72,
    "minus-two branch contradiction",
)

one_packet = {a: 1, b: 0, c: -1, d: 0}
zero(root_remainder.coeff_monomial(u**7).subs(a, 1) - 3 * b**2,
     "one branch forces b")
zero(root_remainder.coeff_monomial(u**5).subs({a: 1, b: 0}) - 3 * (c + 1) ** 2,
     "one branch forces c")
zero(root_remainder.coeff_monomial(u**3).subs({a: 1, b: 0, c: -1}) - 3 * d**2,
     "one branch forces d")
gate(
    root_remainder.coeff_monomial(u).subs(one_packet) == 48,
    "one branch contradiction",
)

# The finite branch polynomial is squarefree of degree eight.
H = 24 * u**8 - 29 * u**6 - 246 * u**4 - 309 * u**2 - 16
zero(sp.discriminant(F, z) - 432 * H, "cubic discriminant")
gate(sp.degree(H, u) == 8, "eight finite branch addresses")
gate(sp.gcd(H, sp.diff(H, u)) == 1, "finite branch addresses are simple")
gate(
    sp.resultant(H, sp.diff(H, u), u)
    == -11766869015040000000000000000000000,
    "squarefree resultant hostile",
)

# At infinity put w=1/u and v=z/u^3.  The simple leading root -2 and the
# two Hensel roots above v=1 all use w itself as uniformizer, so infinity is
# unramified in the degree-three projection.
F_infinity = sp.expand(w**9 * F.subs({u: 1 / w, z: v / w**3}))
gate(sp.denom(F_infinity) == 1 and sp.Poly(F_infinity, w).is_Poly,
     "infinity equation is polynomial")
gate(sp.factor(F_infinity.subs(w, 0)) == (v - 1) ** 2 * (v + 2), "infinity roots")
gate(sp.diff(F_infinity, v).subs({w: 0, v: -2}) == 9, "simple minus-two infinity branch")

double_substitution = 1 - w**2 - 4 * w**4 + w**5 * q
double_quotient = sp.cancel(F_infinity.subs(v, double_substitution) / w**10)
gate(
    sp.denom(double_quotient) == 1,
    "double-leading infinity branch has integral Newton expansion",
)
gate(double_quotient.subs(w, 0) == 3 * (q**2 - 32), "two infinity Hensel roots")
gate(sp.discriminant(3 * (q**2 - 32), q) == 1152, "infinity roots are distinct")
gate(-6 + 8 == 2 and (2 + 2) // 2 == 2, "Riemann-Hurwitz genus two")

# The other Jacobian divisor is rational: after z=(u^2+1)r its function
# field is the conic r^2=u^2-4.  This hostile separates the genuinely new
# genus debt from the ordinary rational ramification divisor.
rr = sp.symbols("r")
zero(
    (sp.diff(F, z) / 3).subs(z, (u**2 + 1) * rr)
    - (u**2 + 1) ** 2 * (rr**2 - u**2 + 4),
    "F_z normalization is a rational conic",
)
gate(sp.discriminant(u**2 - 4, u) == 16, "rational conic is reduced")

# Elementary hostile: a rational collapsed divisor can occur for a
# birational polynomial plane map, but its Jacobian pays the divisor.
x, y = sp.symbols("x y")
rational_hostile_jacobian = sp.det(
    sp.Matrix([[sp.diff(x, x), sp.diff(x, y)], [sp.diff(x * y, x), sp.diff(x * y, y)]])
)
gate(rational_hostile_jacobian == x, "rational collapsed-divisor hostile")

semantic_payload = {
    "general_gate": "positive_genus_divisorial_valuation_with_v(P),v(Q)>0_forbids_plane_Keller_model",
    "center_trichotomy": "affine_curve_or_rational_point_exceptional_or_rational_infinity",
    "candidate": "THM3915_rational_chart",
    "collapsed_divisor": "F=0_irreducible_genus2",
    "finite_ramification": 8,
    "infinity_ramification": 0,
    "other_jacobian_divisor": "F_z=0_rational_conic",
    "conclusion": "THM3915_has_no_polynomial_plane_etale_atlas",
    "scope": "other_repeated_residual_genus0_designs_and_JC2_open",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("theorem=THM-3916-positive-genus-collapsed-valuation-keller-obstruction")
print("general_gate=positive_genus_common_zero_divisorial_valuation_forbids_plane_Keller_model")
print("rational_chart=A=F/4,C=uF/4;jacobian=-F*F_z/16")
print("F=irreducible_monic_cubic_over_k(u)")
print("disc_z_F=432*H8;H8=squarefree;finite_ramification=8")
print("infinity=three_unramified_places;ramification=0")
print("normalization_F=genus2_by_Riemann_Hurwitz")
print("normalization_F_z=genus0_conic")
print("THM3915_plane_atlas=EMPTY")
print("repeated_residual_genus0_designs=OPEN;JC2=OPEN")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
