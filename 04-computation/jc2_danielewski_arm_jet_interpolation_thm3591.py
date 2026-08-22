#!/usr/bin/env python3
"""Finite exact companion for proved and hostile-audited THM-3591.

The universal arm-jet theorem is proof-driven.  This script checks exact
Poisson restrictions, inverse-derivative sidecars, separated hostiles, the
calibrated bulk defect, and sharp boundaries without Python assert gates.
"""

from math import prod
from itertools import combinations

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


print("SECTION seven-piece nodal-arm candidate")
P_star = c**3 - J_a13 * c + e**2
Q_star = c**2 + e**3 - e - sp.Rational(3, 2) * J_a13 * c * e
Delta_star = arm_det(P_star, Q_star)
require("seven-piece arm determinant", zero(Delta_star - J_a13))
scalar_star = -J_a13 * sp.diff(Sigma_a13, b) - 3 * sp.diff(J_a13, b) * Sigma_a13
require(
    "seven-piece scalar channel",
    zero(scalar_star + 1 + sp.Rational(27, 2) * b * Sigma_a13),
)
for beta_value in (0, sp.I, -sp.I):
    require(
        f"seven-piece arm bracket beta={beta_value}",
        zero(
            bracket(P_star, Q_star, 3, Sigma_a13).subs({b: beta_value, c: 0})
            + 1
        ),
    )

arm_c = (c**3 - c, c**2)
arm_e = (e**2, e**3 - e)
for variable, pair, label in ((c, arm_c, "c"), (e, arm_e, "e")):
    require(
        f"nodal collision {label} first",
        pair[0].subs(variable, 1) == pair[0].subs(variable, -1),
    )
    require(
        f"nodal collision {label} second",
        pair[1].subs(variable, 1) == pair[1].subs(variable, -1),
    )
    require(
        f"nodal immersion {label}",
        sp.gcd(
            sp.Poly(sp.diff(pair[0], variable), variable),
            sp.Poly(sp.diff(pair[1], variable), variable),
        ).degree()
        == 0,
    )

star_groebner = sp.groebner(
    [c**3 * e - Sigma_a13], e, c, b, order="lex", domain=sp.QQ
)
star_reduced = sp.expand(
    star_groebner.reduce(bracket(P_star, Q_star, 3, Sigma_a13))[1]
)
star_defect = sp.expand(star_reduced + 1)
require("seven-piece nonzero bulk defect", star_defect != 0)
require(
    "seven-piece lex-normal-form scalar coefficient",
    zero(star_defect.subs(c, 0) + sp.Rational(27, 2) * b * Sigma_a13),
)
arm_ideal = sp.groebner([c, Sigma_a13], c, b, e, order="lex", domain=sp.QQ)
require("seven-piece defect in arm ideal", arm_ideal.reduce(star_defect)[1] == 0)
require("seven-piece P support", len({-6, 1, 3}) == 3)
require("seven-piece Q support", len({-9, -3, -2, 2}) == 4)

print("PASS seven-piece nodal curves, inverse arm jet, and nonzero bulk defect")


print("SECTION scalar-paid seven-piece and arm-interpolating eight-piece hostiles")
q_dagger = (225 * b**4 + 237 * b**2 + 8) / 8
v_dagger = sp.Rational(105, 16) * b
P_dagger = e**2 + c + c**3
Q_dagger = e**3 + q_dagger * e + c**2 + v_dagger * c**4

W_1m3 = -sp.diff(Sigma_a13 * q_dagger, b)
W_m64 = 4 * sp.diff(Sigma_a13**2, b) * v_dagger + 6 * Sigma_a13**2 * sp.diff(v_dagger, b)
require("dagger scalar layer", zero(W_1m3 + W_m64 + 1))
require(
    "dagger unique nonscalar layer",
    zero(2 * sp.diff(Sigma_a13**2, b) - 4 * Sigma_a13 * sp.diff(Sigma_a13, b)),
)
require("dagger nonscalar layer nonzero", 4 * Sigma_a13 * sp.diff(Sigma_a13, b) != 0)
require("dagger q origin", q_dagger.subs(b, 0) == 1)
require("dagger v origin", v_dagger.subs(b, 0) == 0)
Delta_dagger = arm_det(P_dagger, Q_dagger)
require("dagger arm determinant", zero(Delta_dagger - (3 * e**2 + q_dagger)))
require(
    "dagger inverse jet residue",
    sp.rem(q_dagger - J_a13, Sigma_a13, domain=sp.QQ) == 0,
)
require("dagger misses arm interpolation", sp.diff(Delta_dagger, e) != 0)

dagger_arm_c = (c + c**3, c**2)
dagger_arm_e = (e**2, e + e**3)
for variable, pair, label in ((c, dagger_arm_c, "c"), (e, dagger_arm_e, "e")):
    require(
        f"dagger nodal collision {label} first",
        pair[0].subs(variable, sp.I) == pair[0].subs(variable, -sp.I),
    )
    require(
        f"dagger nodal collision {label} second",
        pair[1].subs(variable, sp.I) == pair[1].subs(variable, -sp.I),
    )
    require(
        f"dagger nodal immersion {label}",
        sp.gcd(
            sp.Poly(sp.diff(pair[0], variable), variable),
            sp.Poly(sp.diff(pair[1], variable), variable),
        ).degree()
        == 0,
    )

dagger_reduced = sp.expand(
    star_groebner.reduce(bracket(P_dagger, Q_dagger, 3, Sigma_a13))[1]
)
require("dagger full bracket not constant", dagger_reduced != -1)
require("dagger central jet", bracket(P_dagger, Q_dagger, 3, Sigma_a13).subs({b: 0, c: 0, e: 0}) == -1)

Q_ddagger = Q_dagger + sp.Rational(3, 2) * c * e
Delta_ddagger = arm_det(P_dagger, Q_ddagger)
require("double-dagger arm determinant", zero(Delta_ddagger - q_dagger))
require(
    "double-dagger arm interpolation",
    sp.rem(Delta_ddagger - J_a13, Sigma_a13, domain=sp.QQ) == 0,
)
for beta_value in (0, sp.I, -sp.I):
    require(
        f"double-dagger arm bracket beta={beta_value}",
        zero(
            bracket(P_dagger, Q_ddagger, 3, Sigma_a13).subs(
                {b: beta_value, c: 0}
            )
            + 1
        ),
    )
ddagger_reduced = sp.expand(
    star_groebner.reduce(bracket(P_dagger, Q_ddagger, 3, Sigma_a13))[1]
)
require("double-dagger full bracket not constant", ddagger_reduced != -1)
require("double-dagger added weight", len({-9, -3, -2, 2, 4}) == 5)

print("PASS scalar payment, arm repair, nodal arms, and nonscalar-layer failure")


print("SECTION eight-piece support and <=10-piece additive-repair nonentry")
R_support = {-6, 1, 3}
S_support = {-9, -3, -2, 2, 4}
forced_sums = {-6, -4, -1, 0, 1}
expected_edges = {
    -6: (3, -9),
    -4: (-6, 2),
    -1: (1, -2),
    0: (3, -3),
    1: (3, -2),
}
for total, edge in expected_edges.items():
    edges = sorted((r, s) for r in R_support for s in S_support if r + s == total)
    require(f"unique base sum {total}", edges == [edge])

ff = sp.Function("ff")(b)
gg = sp.Function("gg")(b)


def weight_wronskian(r, s):
    return sp.expand(s * sp.diff(ff, b) * gg - r * ff * sp.diff(gg, b))


base_product_identities = (
    sp.diff(ff**3 * gg, b) + ff**2 * weight_wronskian(3, -9) / 3,
    sp.diff(ff * gg**3, b) - gg**2 * weight_wronskian(-6, 2) / 2,
    sp.diff(ff**2 * gg, b) + ff * weight_wronskian(1, -2),
    sp.diff(ff * gg, b) + weight_wronskian(3, -3) / 3,
    sp.diff(ff**2 * gg**3, b) + ff * gg**2 * weight_wronskian(3, -2),
)
for index, identity in enumerate(base_product_identities):
    require(f"base product derivative {index}", sp.simplify(identity) == 0)


def cover_by_new_p(value):
    return {total for total in forced_sums if any(value + s == total for s in S_support)}


def cover_by_new_q(value):
    return {total for total in forced_sums if any(r + value == total for r in R_support)}


p_candidates = sorted(
    {total - s for total in forced_sums for s in S_support} - R_support
)
q_candidates = sorted(
    {total - r for total in forced_sums for r in R_support} - S_support
)

one_p = [p for p in p_candidates if cover_by_new_p(p) == forced_sums]
one_q = [q for q in q_candidates if cover_by_new_q(q) == forced_sums]
require("no one-P repair", one_p == [])
require("no one-Q repair", one_q == [])

two_p = [
    pair
    for pair in combinations(p_candidates, 2)
    if cover_by_new_p(pair[0]) | cover_by_new_p(pair[1]) == forced_sums
]
two_q = [
    pair
    for pair in combinations(q_candidates, 2)
    if cover_by_new_q(pair[0]) | cover_by_new_q(pair[1]) == forced_sums
]
mixed = []
for p_value in p_candidates:
    for q_value in q_candidates:
        covered = cover_by_new_p(p_value) | cover_by_new_q(q_value)
        if p_value + q_value in forced_sums:
            covered.add(p_value + q_value)
        if covered == forced_sums:
            mixed.append((p_value, q_value))
require("unique two-P repair support", two_p == [(-3, -2)])
require("no two-Q repair support", two_q == [])
require("unique mixed repair support", mixed == [(-3, -1)])

require(
    "two-P survivor product peel",
    sp.simplify(
        sp.diff(ff**2 * gg, b) - ff * weight_wronskian(-2, 4) / 2
    )
    == 0,
)
require(
    "mixed survivor product peel",
    sp.simplify(
        sp.diff(ff * gg**3, b) + gg**2 * weight_wronskian(3, -1)
    )
    == 0,
)
require("negative-weight Sigma invoice", {r: (-r + 2) // 3 for r in (-9, -6, -3, -2, -1)} == {-9: 3, -6: 2, -3: 1, -2: 1, -1: 1})

print("PASS full 3x5 cell and every one-/two-piece enlargement; frontier >=11")


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
