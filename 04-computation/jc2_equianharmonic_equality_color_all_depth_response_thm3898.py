#!/usr/bin/env python3
"""Exact all-depth equality-color response companion for THM-3898."""

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


# q is the exact placeholder epsilon^n.  Keeping it independent of epsilon
# proves the identity at every n rather than sampling response depths.
epsilon, q = sp.symbols("epsilon q", nonzero=True)
U, V = sp.symbols("U V", nonzero=True)
a, L, kappa = sp.symbols("a L kappa")
P = a * L**2

K = (1 + kappa * epsilon**2) / epsilon**2
T = U / q
f = V / q
r = a * T + K * f
A = K * T + a * P * f
B = P * f**2 - T**2

Rhat = 1 + (kappa + a * U / V) * epsilon**2
Ahat = (1 + kappa * epsilon**2) * U + a * P * V * epsilon**2
E = P * V**2 - U**2
scale = q**4 * epsilon**4 / V**2

macro_terms = {
    "L4": L**4,
    "6AL2f": 6 * A * L**2 * f,
    "6PL2f": 6 * P * L**2 * f,
    "2r2L2f": 2 * r**2 * L**2 * f,
    "8AB": 8 * A * B,
    "6PB": 6 * P * B,
    "3r2B": 3 * r**2 * B,
}

normalized_terms = {
    "L4": L**4 * q**4 * epsilon**4 / V**2,
    "6AL2f": 6 * L**2 * q**2 * epsilon**2 * Ahat / V,
    "6PL2f": 6 * P * L**2 * q**3 * epsilon**4 / V,
    "2r2L2f": 2 * L**2 * q * V * Rhat**2,
    "8AB": 8 * q * epsilon**2 * Ahat * E / V**2,
    "6PB": 6 * P * q**2 * epsilon**4 * E / V**2,
    "3r2B": 3 * Rhat**2 * E,
}

for label in macro_terms:
    zero(scale * macro_terms[label] - normalized_terms[label], f"normalized macro {label}")

S = sp.expand(sum(macro_terms.values()))
C_product = sp.cancel(scale * S + 3 * U**2)
all_depth_rhs = sp.cancel(3 * U**2 + sum(normalized_terms.values()))
zero(C_product - all_depth_rhs, "universal all-depth color identity")

# The base arrival schedule is affine in n.  The first marked source is
# always q=epsilon^n; every other nonpersistent macro begins strictly later.
n = sp.symbols("n", integer=True, positive=True)
arrival_schedule = {
    "persistent": 0,
    "marked_r2L2f": n,
    "AB": n + 2,
    "AL2f": 2 * n + 2,
    "PB": 2 * n + 4,
    "PL2f": 3 * n + 4,
    "L4": 4 * n + 4,
}
for label, exponent in arrival_schedule.items():
    if label == "persistent":
        continue
    gate(sp.simplify(exponent - n) >= 0, f"arrival {label} is not before n")
gate(
    all(
        sp.simplify(exponent - n) > 0
        for label, exponent in arrival_schedule.items()
        if label not in ("persistent", "marked_r2L2f")
    ),
    "the marked source is uniquely first",
)

# Canonical equianharmonic leading payment.
d = sp.symbols("d", nonzero=True)
d_relation = sp.Poly(d**2 + 3, d)


def reduce_d(expression: sp.Expr) -> sp.Expr:
    """Reduce a rational expression in the quadratic field d^2=-3."""
    numerator, denominator = sp.together(expression).as_numer_denom()
    numerator = sp.Poly(numerator, d).rem(d_relation).as_expr()
    denominator = sp.Poly(denominator, d).rem(d_relation).as_expr()
    denominator_poly = sp.Poly(denominator, d)
    constant = denominator_poly.coeff_monomial(1)
    linear = denominator_poly.coeff_monomial(d)
    inverse = (constant - linear * d) / (constant**2 + 3 * linear**2)
    reduced = sp.Poly(sp.together(numerator * inverse).as_numer_denom()[0], d).rem(
        d_relation
    ).as_expr()
    reduced_denominator = sp.together(numerator * inverse).as_numer_denom()[1]
    return sp.cancel(reduced / reduced_denominator)


h = (a + 3 * L**2) / 2
u = (3 * L**2 - a) / (2 * d)
zero(h - d * u - a, "canonical minus color")
zero(h + d * u - 3 * L**2, "canonical plus color")
norm_numerator = sp.together(h**2 - 3 * (P - u**2)).as_numer_denom()[0]
zero(sp.Poly(norm_numerator, d).rem(d_relation).as_expr(), "canonical leading norm")

Q = kappa + a * u
Rhat_can = 1 + Q * epsilon**2
subs_can = {U: u, V: 1}
canonical_product = sp.cancel(all_depth_rhs.subs(subs_can))
canonical_base = h**2 * Rhat_can**2 + 3 * u**2
canonical_product = reduce_d(canonical_product)
canonical_base = reduce_d(canonical_base)
zero(reduce_d(canonical_product - canonical_base).subs(q, 0), "canonical persistent square")

# For every sampled n, independently expand the exact identity.  The first
# departure from the persistent square is at epsilon^n with coefficient
# 2L^2.  The symbolic schedule above upgrades the samples to all n>=1.
for degree in range(1, 9):
    extra = reduce_d((canonical_product - canonical_base).subs(q, epsilon**degree))
    numerator = sp.together(extra).as_numer_denom()[0]
    poly = sp.Poly(sp.expand(numerator), epsilon)
    minimum = min(monomial[0] for monomial, coefficient in poly.terms() if coefficient != 0)
    gate(minimum == degree, f"canonical first source depth n={degree}")
    leading = poly.coeff_monomial(epsilon**degree)
    denominator = sp.together(extra).as_numer_denom()[1]
    zero(leading / denominator - 2 * L**2, f"canonical source coefficient n={degree}")

# At n=4, the persistent fourth response cancels internally and the marked
# source forces J4=L^2/h.  The same computation describes the first marked
# response at arbitrary n, but makes no claim about arbitrary sidecars.
J4 = sp.symbols("J4")
rhs_n4 = sp.Poly(
    sp.cancel(canonical_product.subs(q, epsilon**4)), epsilon
).coeff_monomial(epsilon**4)
lhs_n4 = sp.Poly(
    sp.expand((h * Rhat_can + J4 * epsilon**4) ** 2 + 3 * u**2), epsilon
).coeff_monomial(epsilon**4)
zero(reduce_d(rhs_n4 - (h**2 * Q**2 + 2 * L**2)), "n=4 exact fourth source")
zero(
    reduce_d(lhs_n4 - (h**2 * Q**2 + 2 * h * J4)),
    "n=4 square-root response",
)
zero(
    reduce_d((rhs_n4 - lhs_n4).subs(J4, L**2 / h)),
    "rational response positive control",
)

x = sp.symbols("x")
a_x = x + 1
L_x = 9 * x + 4
h_x = sp.expand((a_x + 3 * L_x**2) / 2)
gate(sp.Poly(h_x, x).degree() == 2, "canonical h is a nonunit quadratic")
gate(sp.gcd(sp.Poly(h_x, x), sp.Poly(L_x, x)).degree() == 0, "gcd(h,L)=1")
gate(sp.denom(sp.cancel(L_x**2 / h_x)).has(x), "forced response is nonpolynomial")

# The old abandoned cubic-carrier reservation exposed the same two constants.
# Preserve only the exact constant-field bridge, not the uncanonized carrier.
lam = sp.symbols("lam")
lam_relation = sp.Poly(lam**2 - lam + 1, lam)
d_lam = 2 * lam - 1
zero(sp.Poly(d_lam**2 + 3, lam).rem(lam_relation).as_expr(), "C3 square root")
zero((3 + d_lam) / 2 - (1 + lam), "C3 first color")
zero((3 - d_lam) / 2 - (2 - lam), "C3 second color")
zero(
    sp.Poly((2 - lam) * (1 + lam) - 3, lam).rem(lam_relation).as_expr(),
    "C3 color product",
)

semantic = {
    "identity": "seven normalized macros give exact all-depth color product",
    "schedule": "0,n,n+2,2n+2,2n+4,3n+4,4n+4",
    "canonical": "first marked response forces L2/h at depth n",
    "n4": "fourth response is nonpolynomial for canonical zero sidecars",
    "hostile": "L2/h is a valid rational response over k(x)",
    "scope": "arbitrary sidecars and equality-seam existence remain open",
    "lineage": "d=2lambda-1 and (2-lambda)(1+lambda)=3",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3898-equianharmonic-equality-color-all-depth-response")
print("identity=seven_normalized_macros_exact_at_all_depths")
print("arrival_schedule=0,n,n+2,2n+2,2n+4,3n+4,4n+4")
print("canonical_first_marked_response=J_n=L2/h")
print("n4_zero_sidecar_control=NONPOLYNOMIAL_FOURTH_RESPONSE")
print("rational_x_hostile=J_n=L2/h_exists")
print("arbitrary_sidecars_and_JC2=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
