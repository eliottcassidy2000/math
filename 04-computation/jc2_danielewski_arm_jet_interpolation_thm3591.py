#!/usr/bin/env python3
"""Finite exact companion for provisional THM-3591.

The universal arm-jet theorem is proof-driven.  This script checks exact
Poisson restrictions, inverse-derivative sidecars, separated hostiles, the
calibrated bulk defect, and sharp boundaries without Python assert gates.
"""

from math import prod

import sympy as sp


b, c, e = sp.symbols("b c e")
CHECKS = 0


def require(label, condition):
    """Record one active exact gate and fail with a stable label."""
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError("FAILED: " + label)


def zero(expr):
    return sp.cancel(sp.expand(expr)) == 0


def bracket(F, G, n, Sigma):
    """Poisson bracket on c^n e=Sigma before quotient reduction."""
    return sp.expand(
        c**n * (sp.diff(F, b) * sp.diff(G, c) - sp.diff(F, c) * sp.diff(G, b))
        - sp.diff(Sigma, b)
        * (sp.diff(F, c) * sp.diff(G, e) - sp.diff(F, e) * sp.diff(G, c))
        - n
        * c ** (n - 1)
        * e
        * (sp.diff(F, b) * sp.diff(G, e) - sp.diff(F, e) * sp.diff(G, b))
    )


def arm_det(F, G):
    return sp.expand(
        sp.diff(F, c).subs(c, 0) * sp.diff(G, e).subs(c, 0)
        - sp.diff(F, e).subs(c, 0) * sp.diff(G, c).subs(c, 0)
    )


print("THM-3591 exact companion")
print("SECTION arm restriction and inverse-derivative sidecars")

ARM_ROWS = 0
for d in range(2, 8):
    roots = tuple(range(d))
    Sigma = sp.expand(prod(b - beta for beta in roots))
    Sigma_p = sp.diff(Sigma, b)
    J = sp.rem(sp.invert(Sigma_p, Sigma), Sigma, domain=sp.QQ)

    require(f"squarefree d={d}", sp.gcd(sp.Poly(Sigma, b), sp.Poly(Sigma_p, b)).degree() == 0)
    require(f"inverse derivative d={d}", sp.rem(J * Sigma_p - 1, Sigma, domain=sp.QQ) == 0)
    require(f"sidecar degree d={d}", 1 <= sp.degree(J, b) < d)
    require(
        f"reciprocal sum d={d}",
        sum(sp.Rational(1, Sigma_p.subs(b, beta)) for beta in roots) == 0,
    )

    for n in range(2, 8):
        # A deliberately broad polynomial pair; no Darboux assumption is used
        # for the universal arm-restriction identity itself.
        F = b**2 + c + c**3 + e**2 + Sigma * (b * c + e**2)
        G = b**3 + 2 * c + c**2 + e + e**3 + Sigma * (c**2 + b * e)
        Delta = arm_det(F, G)
        for beta in roots:
            lhs = bracket(F, G, n, Sigma).subs({b: beta, c: 0})
            rhs = -Sigma_p.subs(b, beta) * Delta.subs(b, beta)
            require(f"arm bracket n={n} d={d} beta={beta}", zero(lhs - rhs))
            require(
                f"J value n={n} d={d} beta={beta}",
                zero(J.subs(b, beta) * Sigma_p.subs(b, beta) - 1),
            )
            ARM_ROWS += 2

print(f"PASS {ARM_ROWS} arm restrictions and sidecar values")


print("SECTION explicit sidecars and calibrated bulk debt")
Sigma_a13 = b * (b**2 + 1)
J_a13 = 1 + sp.Rational(3, 2) * b**2
require("A13 inverse", sp.rem(J_a13 * sp.diff(Sigma_a13, b) - 1, Sigma_a13) == 0)

Sigma_quad = b * (b + 4)
J_quad = (b + 2) / 8
require("quadratic inverse", sp.rem(J_quad * sp.diff(Sigma_quad, b) - 1, Sigma_quad) == 0)

BULK_ROWS = 0
for Sigma, J, label in ((Sigma_a13, J_a13, "A13"), (Sigma_quad, J_quad, "quadratic")):
    for n in range(2, 8):
        raw = bracket(c, -J * e, n, Sigma)
        reduced = sp.expand(raw.subs(c**n * e, Sigma))
        expected = sp.diff(Sigma * J, b)
        require(f"calibrated bracket {label} n={n}", zero(reduced - expected))
        require(f"bulk divisible {label} n={n}", sp.rem(expected - 1, Sigma) == 0)
        require(f"bulk nonconstant {label} n={n}", sp.degree(expected, b) >= 1)
        BULK_ROWS += 3

print(f"PASS explicit inverse jets and {BULK_ROWS} calibrated bulk controls")


print("SECTION separated and Sigma-corrected arm-blind controls")
BLIND_ROWS = 0
for d in range(2, 7):
    roots = tuple(range(d))
    Sigma = sp.expand(prod(b - beta for beta in roots))
    Sigma_p = sp.diff(Sigma, b)
    R = b * c + e**2 + b * c * e
    S = c**2 + b * e + c * e**2
    F = (b**2 + b + 1) + (c + c**4) + (e**2 + e**3) + Sigma * R
    G = (b**3 - b) + (2 * c + c**3) + (e + e**4) + Sigma * S
    Delta = arm_det(F, G)
    values = [sp.expand(Delta.subs(b, beta)) for beta in roots]
    require(f"blind common jet d={d}", all(zero(value - values[0]) for value in values))
    slopes = [Sigma_p.subs(b, beta) for beta in roots]
    require(f"distinct arm slopes d={d}", len(set(slopes)) > 1)
    products = [sp.expand(-slope * value) for slope, value in zip(slopes, values)]
    require(f"blind bracket mismatch d={d}", len(set(map(str, products))) > 1)
    BLIND_ROWS += 3

print(f"PASS {BLIND_ROWS} separated/Sigma-corrected arm-blind gates")


print("SECTION sharp boundaries")
for n in range(2, 8):
    require(f"degree-one plane n={n}", zero(bracket(c, -e, n, b) - 1))

Sigma_n1 = b * (b - 1)
require("exponent-one surviving term", zero(bracket(b, e, 1, Sigma_n1) + e))

CHAR_ROWS = 0
for prime in (2, 3, 5, 7):
    Sigma_p = b**prime - b
    deriv = sp.diff(Sigma_p, b)
    require(
        f"char-p squarefree p={prime}",
        sp.gcd(sp.Poly(Sigma_p, b, modulus=prime), sp.Poly(deriv, b, modulus=prime)).degree() == 0,
    )
    raw = bracket(c, e, 2, Sigma_p)
    coeffs = sp.Poly(sp.expand(raw - 1), b, c, e).coeffs()
    require(f"char-p Darboux hostile p={prime}", all(int(coef) % prime == 0 for coef in coeffs))
    CHAR_ROWS += 2

print(f"PASS degree-one, exponent-one, and {CHAR_ROWS} characteristic-p hostile gates")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
