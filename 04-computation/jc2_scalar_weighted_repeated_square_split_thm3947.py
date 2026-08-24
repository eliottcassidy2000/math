#!/usr/bin/env python3
"""Exact companion for THM-3947's scalar-weighted repeated-square split."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


P, G, x, t, Z, r = sp.symbols("P G x t Z r")
d, lam, j = sp.symbols("delta lambda iota")


def dreduce(expression: sp.Expr) -> sp.Expr:
    """Reduce a polynomial expression modulo delta^2+3."""
    numerator, denominator = sp.fraction(sp.cancel(expression))
    reduced = sp.rem(sp.Poly(sp.expand(numerator), d), sp.Poly(d**2 + 3, d))
    return sp.cancel(reduced.as_expr() / denominator)


def dzero(expression: sp.Expr, message: str) -> None:
    gate(dreduce(expression) == 0, message)


omega = (-1 + d) / 2
omega2 = (-1 - d) / 2
dzero(omega**2 + omega + 1, "omega quadratic relation")
dzero(omega * omega2 - 1, "conjugate roots multiply to one")
dzero(omega - omega2 - d, "delta is omega-omega^2")


# ---------------------------------------------------------------------------
# Arbitrarily weighted internal split of p1-p0=G^2.
# Q_i=lambda*q_i clears the only parameter denominator, with t=lambda^2.
# ---------------------------------------------------------------------------

p0 = P
p1 = P + G**2
L1 = p1 - omega * p0
L2 = p1 - omega2 * p0
Q0 = G * (L2 - t * L1)
Q1 = G * (L2 + t * L1)

gate(p1 - p0 == G**2, "repeated-square difference")
dzero(Q1 - Q0 - 2 * t * G * L1,
      "scaled q1-q0 internal factor")
dzero(Q1 + Q0 - 2 * G * L2,
      "scaled q1+q0 internal factor")
dzero(L1 * L2 * G**2 - (p1**3 - p0**3),
      "difference-of-cubes factor packet")
dzero(Q1**2 - Q0**2 - 4 * t * (p1**3 - p0**3),
      "two discriminant rows agree")

Hscaled0 = sp.expand(Q0**2 - 4 * t * p0**3)
Hscaled1 = sp.expand(Q1**2 - 4 * t * p1**3)
dzero(Hscaled0 - Hscaled1, "common scaled discriminant")


# ---------------------------------------------------------------------------
# Weighted one-variable root cubic and its complete parameter discriminant.
# ---------------------------------------------------------------------------

h = sp.expand(
    (1 + (1 - omega2) * x - t * (1 + (1 - omega) * x)) ** 2
    - 4 * t * x**3
)
Hhom = sp.expand(
    G**2 * (G**2 + (1 - omega2) * P
            - t * (G**2 + (1 - omega) * P)) ** 2
    - 4 * t * P**3
)
dzero(Hscaled0 - Hhom,
      "weighted homogenization of h_t")
dzero(Hhom - sp.expand(G**6 * h.subs(x, P / G**2)),
      "ratio P/G^2 recovers h_t")
dzero(sp.Poly(h, x).LC() + 4 * t,
      "root cubic leading coefficient is -4t")
gate(sp.expand(h.subs(x, 0) - (1 - t) ** 2) == 0,
     "root cubic constant coefficient")

disc_h = sp.discriminant(h, x)
disc_expected = 48 * d * t * (t - 1) ** 3 * (t + omega) ** 3
dzero(disc_h - disc_expected,
      "complete root-cubic discriminant factorization")


# ---------------------------------------------------------------------------
# The two and only two repeated-root endpoints.
# ---------------------------------------------------------------------------

h_one = dreduce(h.subs(t, 1))
h_other = dreduce(h.subs(t, -omega))
dzero(h_one + x**2 * (4 * x + 3),
      "t=1 has doubled p0 root")
dzero(h_other - omega * (x + 1) ** 2 * (4 * x + 1),
      "t=-omega has doubled p1 root")

gate(h_one.subs(x, 0) == 0, "t=1 double root lies on p0")
gate(sp.diff(h_one, x).subs(x, 0) == 0,
     "t=1 root has multiplicity at least two")
gate(sp.diff(h_one, x, 2).subs(x, 0) != 0,
     "t=1 root is not triple")
gate(h_one.subs(x, sp.Rational(-3, 4)) == 0,
     "t=1 remaining root")
gate(sp.diff(h_one, x).subs(x, sp.Rational(-3, 4)) != 0,
     "t=1 remaining root is simple")

dzero(h_other.subs(x, -1), "t=-omega double root lies on p1")
dzero(sp.diff(h_other, x).subs(x, -1),
      "t=-omega root has multiplicity at least two")
gate(dreduce(sp.diff(h_other, x, 2).subs(x, -1)) != 0,
     "t=-omega root is not triple")
dzero(h_other.subs(x, sp.Rational(-1, 4)),
      "t=-omega remaining root")
gate(dreduce(sp.diff(h_other, x).subs(x, sp.Rational(-1, 4))) != 0,
     "t=-omega remaining root is simple")

H_one = dreduce(Hscaled0.subs(t, 1))
# Since H=Hscaled/t and 1/(-omega)=-omega^2 at the second endpoint.
H_other = dreduce(-omega2 * Hscaled0.subs(t, -omega))
dzero(H_one + P**2 * (4 * P + 3 * G**2),
      "first endpoint full discriminant")
dzero(H_other + (P + G**2) ** 2 * (4 * P + G**2),
      "second endpoint full discriminant")

# Generic hostile control: t=2 is away from both exceptional values.
disc_two = dreduce(disc_h.subs(t, 2))
gate(disc_two != 0, "generic specialization has three distinct roots")

# Forbidden boundary: t=0 removes the cubic term but is unavailable because
# lambda^-1 occurs in the original split.
h_zero = dreduce(h.subs(t, 0))
dzero(h_zero - (1 + (1 - omega2) * x) ** 2,
      "excluded t=0 boundary is only a square quadratic")
gate(sp.Poly(h_zero, x).degree() == 2,
     "excluded t=0 loses cubic degree")


# ---------------------------------------------------------------------------
# Every individual root factor is A1/one-place, while their union is not.
# ---------------------------------------------------------------------------

component = P * Z - r * G**2
gate(component.subs({P: 1, G: 0, Z: 0}) == 0,
     "nonzero-r parabola has the advertised infinity point")
gate(sp.diff(component, Z).subs({P: 1, G: 0, Z: 0}) == 1,
     "parabola is smooth at infinity")
gate(sp.diff(P - r * G**2, P) == 1,
     "every affine root component is smooth")
r1, r2 = sp.symbols("r1 r2")
gate(sp.expand((P - r1 * G**2) - (P - r2 * G**2)
               - (r2 - r1) * G**2) == 0,
     "distinct root parabolas meet only at the origin")


# ---------------------------------------------------------------------------
# The t=-omega endpoint is the t=1 endpoint after swapping presentations.
# Reduce modulo delta^2+3, iota^2+1 and lambda^2+omega.
# ---------------------------------------------------------------------------

def endpoint_reduce(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    reduced = sp.rem(
        sp.Poly(sp.expand(numerator), lam),
        sp.Poly(lam**2 + omega, lam),
    ).as_expr()
    reduced = sp.rem(
        sp.Poly(sp.expand(reduced), j),
        sp.Poly(j**2 + 1, j),
    ).as_expr()
    reduced = sp.rem(
        sp.Poly(sp.expand(reduced), d),
        sp.Poly(d**2 + 3, d),
    ).as_expr()
    return sp.cancel(reduced / denominator)


def endpoint_zero(expression: sp.Expr, message: str) -> None:
    gate(endpoint_reduce(expression) == 0, message)


# At lambda^2=-omega, lambda^-1=-omega^2*lambda.
lambda_inv = -omega2 * lam
q0_endpoint = G * (lambda_inv * L2 - lam * L1)
q1_endpoint = G * (lambda_inv * L2 + lam * L1)

pp0, pp1 = p1, p0
qq0, qq1 = q1_endpoint, q0_endpoint
Gp = j * G
omegap, omegap2 = omega2, omega
lambdap = -j * omega * lam
LL1 = pp1 - omegap * pp0
LL2 = pp1 - omegap2 * pp0

endpoint_zero(pp1 - pp0 - Gp**2,
              "swapped repeated square uses G'=iota G")
endpoint_zero(lambdap**2 - 1,
              "second endpoint scalar becomes t'=1")
endpoint_zero(LL1 + omega2 * L1,
              "swapped first cube-difference factor")
endpoint_zero(LL2 + omega * L2,
              "swapped second cube-difference factor")
endpoint_zero(qq1 - qq0 - 2 * lambdap * Gp * LL1,
              "endpoint swapped q-difference row")
# Since lambdap^2=1, its inverse equals itself.
endpoint_zero(qq1 + qq0 - 2 * lambdap * Gp * LL2,
              "endpoint swapped q-sum row")
endpoint_zero(
    H_other + pp0**2 * (4 * pp0 + 3 * Gp**2),
    "endpoint full discriminant is THM-3944 after the swap",
)


summary = {
    "checks": CHECKS,
    "grammar": "arbitrary reciprocal scalar on repeated square G^2",
    "root_cubic_discriminant": "48*delta*t*(t-1)^3*(t+omega)^3",
    "generic_support": "three distinct smooth A1 parabolas",
    "endpoint_t_1": "-P^2(4P+3G^2)",
    "endpoint_t_minus_omega": "-(P+G^2)^2(4P+G^2)",
    "triple_root_parameters": 0,
    "full_discriminant_irreducible_parameters": 0,
    "endpoint_equivalence": "swap torus rows; omega->omega^2; G->iota*G",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3947 scalar-weighted repeated-square split companion")
print(f"CHECKS={CHECKS}")
print("ROOT_CUBIC_DISC=48*delta*t*(t-1)^3*(t+omega)^3")
print("GENERIC=THREE_DISTINCT_A1_COMPONENTS")
print("T=1:-P^2(4P+3G^2);MULTIPLICITIES=2+1")
print("T=-omega:-(P+G^2)^2(4P+G^2);MULTIPLICITIES=2+1")
print("TRIPLE_ROOT_PARAMETERS=0")
print("FULL_DISCRIMINANT_IRREDUCIBLE_PARAMETERS=0")
print("ENDPOINTS=EQUIVALENT_BY_ROW_SWAP_AND_G_RESCALING")
print("HOSTILE_T=0:EXCLUDED_AND_CUBIC_DEGREE_DROPS_TO_TWO")
print("SCOPE=INDIVIDUAL_COMPONENTS_ONE_PLACE;FULL_BRANCH_ALWAYS_REDUCIBLE")
print(f"SEMANTIC_SHA256={semantic}")
