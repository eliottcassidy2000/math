#!/usr/bin/env python3
"""THM-3905 third-response hostile for the named THM-3902 n=1 lift."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


x, y, d = sp.symbols("x y d")
modulus = sp.Poly(d**2 + 3, d)


def reduce_d(expr: sp.Expr) -> sp.Expr:
    """Reduce a rational expression in Q(x)[d]/(d^2+3)."""
    numerator, denominator = sp.fraction(sp.cancel(expr))
    numerator = sp.rem(sp.Poly(sp.expand(numerator), d), modulus).as_expr()
    denominator = sp.rem(sp.Poly(sp.expand(denominator), d), modulus).as_expr()
    conjugate = denominator.xreplace({d: -d})
    norm = sp.rem(sp.Poly(sp.expand(denominator * conjugate), d), modulus).as_expr()
    require(norm != 0, "nonzero quadratic-field denominator")
    reduced = sp.rem(sp.Poly(sp.expand(numerator * conjugate), d), modulus).as_expr()
    return sp.cancel(reduced / norm)


a = x + 1
L = 9 * x + 4
kappa = -15 * x**2 - 15 * x - 4
K = y**2 + kappa
P = a * L**2

h = (a + 3 * L**2) / 2
u = -d * (3 * L**2 - a) / 6  # (3L^2-a)/(2d), since d^2=-3
p = -2 * (1 - 4 * d) / 7203  # -2/[147(1+4d)]
v1 = reduce_d(p * h)
j1 = reduce_d(a * v1 + sp.Rational(1, 3))
u1 = reduce_d(d * j1 / 3)
j2 = reduce_d(
    (kappa + a * u) * h + (L**2 * p / 2) * (3 * a * p * h + 2)
)
g1 = reduce_d(h * v1 + j1)
g2 = reduce_d(j2 + v1 * j1)

require(reduce_d(h**2 - 3 * (a * L**2 - u**2)) == 0, "leading norm")
require(reduce_d(u1.subs(x, 0) - 4 * v1.subs(x, 0)) == 0, "origin address")

T = u * y + u1
f = y + v1
r = a * T + K * f
A = K * T + a * P * f
S = sp.expand(
    L**4
    + 2 * (3 * A + 3 * P + r**2) * L**2 * f
    + (8 * A + 6 * P + 3 * r**2) * (P * f**2 - T**2)
)
S5 = reduce_d(sp.Poly(S, y).coeff_monomial(y**5))
g3_direct = reduce_d((S5 - 2 * g1 * g2) / (2 * h))

# Independent extraction through the epsilon^3 color coefficient.
E0 = reduce_d(a * L**2 - u**2)
E1 = reduce_d(2 * a * L**2 * v1 - 2 * u * u1)
ratio1 = reduce_d(u1 - u * v1)
B3 = reduce_d(
    4 * kappa * L**2
    + 6 * ((kappa + a * u) * E1 + a * ratio1 * E0)
    + 12 * a * L**2 * u
    - 8 * u**3
)
j3 = reduce_d((B3 - 2 * j1 * j2) / (2 * h))
g3_color = reduce_d(j3 + v1 * j2)
require(reduce_d(g3_direct - g3_color) == 0, "direct/color third response")

numerator, denominator = sp.fraction(sp.cancel(g3_direct))
numerator = sp.rem(sp.Poly(sp.expand(numerator), d), modulus).as_expr()
components = sp.Poly(numerator, d)
h_polynomial = 2 * h
rem_one = sp.factor(
    sp.rem(sp.Poly(components.coeff_monomial(1), x), sp.Poly(h_polynomial, x)).as_expr()
)
rem_d = sp.factor(
    sp.rem(sp.Poly(components.coeff_monomial(d), x), sp.Poly(h_polynomial, x)).as_expr()
)
require(rem_one != 0 and rem_d != 0, "nonpolynomial third response")
require(sp.expand(rem_d + 4 * rem_one) == 0, "quadratic-component remainder ratio")

semantic = {
    "control": "THM-3902 address-compatible n=1 positive two-jet",
    "third_equation": "[y^5]S=2*h*g3+2*g1*g2",
    "denominator": str(sp.factor(denominator)),
    "remainder_mod_2h": str(sp.factor(rem_one + d * rem_d)),
    "consequence": "the named two-jet control has no polynomial third lift",
    "scope": "does not close the n=1 seam or JC2",
}
blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("THM3905_N1_NAMED_CONTROL_THIRD_RESPONSE_HOSTILE")
print("status=PASS;named_two_jet_control_dies_at_third_polynomial_response")
print(f"g3_denominator={sp.factor(denominator)}")
print(f"numerator_remainder_one={rem_one}")
print(f"numerator_remainder_d={rem_d}")
print(f"semantic_sha256={hashlib.sha256(blob).hexdigest()}")
print("ALL CHECKS PASSED")
