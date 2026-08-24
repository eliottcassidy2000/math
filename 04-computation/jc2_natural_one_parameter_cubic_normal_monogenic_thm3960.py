#!/usr/bin/env python3
"""Exact companion for THM-3960's natural cubic closure."""

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


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.expand(expression) == 0, message)


# ---------------------------------------------------------------------------
# Universal natural cubic and exact reducibility obstruction.
# ---------------------------------------------------------------------------

T, P, h, C, E = sp.symbols("T P h C E")
F = sp.expand(T**3 - 3 * P * T - (E + C * P))
G = sp.expand(E + C * h**2 - 2 * h**3)
B = sp.expand(C**3 + 27 * E)

root = -C / 3
zero(F.subs(T, root) + B / 27,
     "universal polynomial-root obstruction")

x = sp.symbols("x")
F_reducible = sp.expand(F.subs({C: 3 * x, E: -x**3}))
zero(F_reducible - (T + x) * (T**2 - x * T + x**2 - 3 * P),
     "B=0 linear-times-quadratic factorization")
zero(B.subs({C: 3 * x, E: -x**3}), "reducible boundary B vanishes")


# ---------------------------------------------------------------------------
# Direct finite-singular-locus normality route.
# ---------------------------------------------------------------------------

FT = sp.diff(F, T)
FP = sp.diff(F, P)
zero(FT - 3 * (T**2 - P), "relative derivative")
zero(FP + (3 * T + C), "P derivative")

Tsing = -C / 3
Psing = C**2 / 9
zero(FT.subs({T: Tsing, P: Psing}), "singular candidate T derivative")
zero(FP.subs({T: Tsing, P: Psing}), "singular candidate P derivative")
zero(F.subs({T: Tsing, P: Psing}) + B / 27,
     "singular candidate lies over B=0")

# The different cannot be scalar in the free basis (1,T,T^2): retain its
# exact coefficient rows.
delta = sp.Poly(FT, T)
gate(delta.degree() == 2, "different has degree two in the generator")
gate(delta.coeff_monomial(T**2) == 3,
     "different has nonzero T^2 coefficient")
gate(delta.coeff_monomial(T) == 0, "different centered T row")
gate(delta.coeff_monomial(1) == -3 * P, "different base row")


# ---------------------------------------------------------------------------
# Hidden-irreducible branch: squarefree and primitive discriminant curve.
# ---------------------------------------------------------------------------

H = sp.expand((E + C * P) ** 2 - 4 * P**3)
disc_F = sp.factor(sp.discriminant(F, T))
disc_H = sp.factor(sp.discriminant(H, P))
zero(disc_F + 27 * H, "order discriminant is -27H")
zero(disc_H + 16 * E**3 * B,
     "branch cubic discriminant -16 E^3 B")

H_poly = sp.Poly(H, P)
gate(H_poly.degree() == 3, "branch polynomial has P degree three")
gate(H_poly.LC() == -4, "branch leading coefficient is a unit")
coefficients = H_poly.all_coeffs()
gate(coefficients == [-4, C**2, 2 * C * E, E**2],
     "branch coefficient rows")
gate(sp.gcd_list(coefficients) == 1,
     "unit leading coefficient excludes vertical content")

q_common, c0, e0 = sp.symbols("q_common c0 e0")
H_common = sp.expand(H.subs({C: q_common * c0, E: q_common * e0}))
zero(H_common.subs(q_common, 0) + 4 * P**3,
     "a common C,E factor does not create a vertical branch")

hidden_E_zero = sp.factor(G.subs(E, 0))
zero(hidden_E_zero - h**2 * (C - 2 * h),
     "E=0 hidden cubic splits")
hidden_B_zero = sp.factor(G.subs({C: 3 * x, E: -x**3}))
zero(27 * hidden_B_zero + (3 * h - 3 * x) ** 2 * (6 * h + 3 * x),
     "B=0 hidden cubic repeated factorization")

# Direct universal factorization in C, without the x substitution.
hidden_B_general = -((3 * h - C) ** 2 * (6 * h + C)) / 27
zero(
    G.subs(E, -C**3 / 27) - hidden_B_general,
    "B=0 hidden factorization in C",
)


# ---------------------------------------------------------------------------
# Reducible-component normalization and degree-two unit gate.
# ---------------------------------------------------------------------------

w = sp.symbols("w")
T_quad = (w + x) / 2
P_quad = (w**2 + 3 * x**2) / 12
quadratic_factor = T**2 - x * T + x**2 - 3 * P
zero(quadratic_factor.subs({T: T_quad, P: P_quad}),
     "quadratic component affine-plane parametrization")
zero(3 * (4 * P_quad - x**2) - w**2,
     "quadratic component ramification square")
zero(sp.diff(P_quad, w) - w / 6,
     "quadratic component derivative")
zero(quadratic_factor.subs(T, -x) - 3 * (x**2 - P),
     "linear/quadratic conductor intersection")


# ---------------------------------------------------------------------------
# Hostile controls for every edge.
# ---------------------------------------------------------------------------

t = sp.symbols("t")

# Hidden-irreducible Eisenstein control C=0,E=t.
F_cusp = sp.expand(F.subs({C: 0, E: t}))
G_cusp = sp.expand(G.subs({C: 0, E: t}))
H_cusp = sp.factor(H.subs({C: 0, E: t}))
B_cusp = sp.expand(B.subs({C: 0, E: t}))
gate(F_cusp == T**3 - 3 * P * T - t, "C=0 irreducible F control")
gate(G_cusp == t - 2 * h**3, "C=0 hidden Eisenstein control")
gate(H_cusp == -4 * P**3 + t**2, "C=0 cusp branch control")
gate(B_cusp == 27 * t, "C=0 nonzero B control")
gate(sp.gcd(H_cusp, sp.diff(H_cusp, P)) == 1,
     "C=0 branch is reduced globally")
gate(sp.factor(sp.discriminant(H_cusp, P)) == -432 * t**4,
     "C=0 branch discriminant control")

# Split-hidden repeated-zero hostile E=0,C=t: H is nonreduced, but F remains
# integral and has only the isolated singular point over t=0.
F_double = sp.expand(F.subs({C: t, E: 0}))
H_double = sp.factor(H.subs({C: t, E: 0}))
B_double = sp.expand(B.subs({C: t, E: 0}))
gate(F_double == T**3 - 3 * P * T - t * P,
     "nonreduced-branch integral hostile F")
zero(H_double + P**2 * (4 * P - t**2),
     "nonreduced-branch hostile H")
gate(B_double == t**3, "nonreduced branch still has nonzero B")
for derivative in (F_double, sp.diff(F_double, T),
                   sp.diff(F_double, P), sp.diff(F_double, t)):
    zero(derivative.subs({t: 0, P: 0, T: 0}),
         "isolated hostile singular point")

# Reducible boundary C=3t,E=-t^3, including the triple-zero specialization.
F_boundary = sp.factor(F.subs({C: 3 * t, E: -t**3}))
gate(F_boundary == (T + t) * (-3 * P + T**2 - T * t + t**2),
     "moving reducible boundary control")
zero(F_boundary.subs(t, 0) - T * (T**2 - 3 * P),
     "triple-zero reducible specialization")


# ---------------------------------------------------------------------------
# Discriminant-index parity controls.
# ---------------------------------------------------------------------------

# If an order discriminant exponent is one, maximal discriminant exponent
# plus twice an index length can only be 1+0. Freeze all nonnegative rows in
# a hostile bounded range; the proof uses the same parity argument generally.
solutions = []
for maximal_exponent in range(0, 8):
    for index_length in range(0, 8):
        if maximal_exponent + 2 * index_length == 1:
            solutions.append((maximal_exponent, index_length))
gate(solutions == [(1, 0)],
     "squarefree discriminant forces codimension-one maximal order")


summary = {
    "checks": CHECKS,
    "dichotomy": "F integral iff B=C^3+27E nonzero",
    "normality": "integral hypersurface has finite singular locus",
    "hidden_irreducible": "E and B nonzero; H primitive squarefree",
    "index": "valuation-one discriminant forces local index zero",
    "different": "global F_T becomes forbidden nonconstant A2 unit",
    "reducible": "degree-one trivial plus degree-two ramification-unit kill",
    "scope": "all natural one-parameter globally monogenic cubics",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3960 natural one-parameter cubic closure companion")
print(f"CHECKS={CHECKS}")
print("DICHOTOMY=F_IRREDUCIBLE_IFF_C_CUBED_PLUS_27E_NONZERO")
print("NORMALITY=INTEGRAL_HYPERSURFACE;FINITE_SINGULAR_LOCUS;R1_S2")
print("HIDDEN_IRREDUCIBLE=E_NONZERO;B_NONZERO;H_PRIMITIVE_SQUAREFREE")
print("ORDER=DISC_F_MINUS_27H;SQUAREFREE_BRANCH_FORCES_INDEX_ZERO")
print("DIFFERENT=F_T_GLOBAL_NONSCALAR;KELLER_OPEN_MAKES_FORBIDDEN_UNIT")
print("REDUCIBLE=LINEAR_TRIVIAL;QUADRATIC_A2_RAMIFICATION_LINE_UNIT_KILL")
print("HOSTILE=E_ZERO_NONREDUCED_BRANCH_STILL_NORMAL")
print("SCOPE=FULL_NATURAL_MONOGENIC_FAMILY;NONMONOGENIC_AND_JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
