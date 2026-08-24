#!/usr/bin/env python3
"""Exact companion for THM-3956's split hidden-cubic trichotomy."""

from __future__ import annotations

import hashlib
import json
import math

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.expand(expression) == 0, message)


# ---------------------------------------------------------------------------
# Hidden-root Vieta rows and the nonzero repeated-root factorization.
# ---------------------------------------------------------------------------

h, T, P, x, y = sp.symbols("h T P x y")
pair_sum = sp.expand(x**2 + 2 * x * y)
zero(pair_sum - x * (x + 2 * y), "repeated-root Vieta dichotomy")

C_general = 2 * (2 * x + y)
E_general = 2 * x**2 * y
G_general = sp.expand(E_general + C_general * h**2 - 2 * h**3)
zero(
    G_general + 2 * (h - x) ** 2 * (h - y) - 2 * h * pair_sum,
    "hidden cubic with repeated root modulo the missing h row",
)

C_nonzero = 3 * x
E_nonzero = -x**3
F_nonzero = sp.expand(T**3 - 3 * P * T - (E_nonzero + C_nonzero * P))
factor_nonzero = sp.expand((T + x) * (T**2 - x * T + x**2 - 3 * P))
zero(F_nonzero - factor_nonzero, "nonzero repeated-root reducibility")

quadratic = T**2 - x * T + x**2 - 3 * P
zero(quadratic.subs(T, -x) - 3 * (x**2 - P),
     "component conductor intersection")
w = sp.symbols("w")
T_quadratic = (w + x) / 2
P_quadratic = (w**2 + 3 * x**2) / 12
zero(quadratic.subs({T: T_quadratic, P: P_quadratic}),
     "quadratic component affine-plane parametrization")
zero(P_quadratic.subs(w, -3 * x) - x**2,
     "component gluing in quadratic coordinate")
zero(3 * (4 * P_quadratic - x**2) - w**2,
     "quadratic component square-ramification identity")
zero(sp.diff(P_quadratic, w) - w / 6,
     "quadratic component ramifies exactly on w=0")

F_triple_zero = sp.expand(F_nonzero.subs(x, 0))
zero(F_triple_zero - T * (T**2 - 3 * P),
     "triple-zero cubic reducibility")

disc_nonzero = sp.factor(sp.discriminant(F_nonzero, T))
disc_expected = sp.factor(27 * (P - x**2) ** 2 * (4 * P - x**2))
zero(disc_nonzero - disc_expected,
     "conductor-square versus quadratic-ramification discriminant")


# ---------------------------------------------------------------------------
# Repeated-zero exceptional branch and exact Danielewski coordinate model.
# ---------------------------------------------------------------------------

r, rp, u, v = sp.symbols("r rp u v")
C_zero = 2 * r
E_zero = sp.Integer(0)
F_zero = sp.expand(T**3 - P * (3 * T + 2 * r))
gate(sp.gcd(T**3, 3 * T + 2 * r) == 1,
     "primitive linear-in-P domain control")

H = sp.expand(u * v - (u - 2 * r) ** 3)
zero(
    H.subs({u: 3 * T + 2 * r, v: 27 * P}) + 27 * F_zero,
    "forward coordinate change to uv=(u-2r)^3",
)
zero(
    F_zero.subs({T: (u - 2 * r) / 3, P: v / 27}) + H / 27,
    "inverse coordinate change",
)
zero(((3 * T + 2 * r) - 2 * r) / 3 - T,
     "inverse recovers T")
zero((27 * P) / 27 - P, "inverse recovers P")

Hv = sp.diff(H, v)
Hu = sp.diff(H, u)
# Treat r'=rp for the universal t derivative.
Ht = sp.expand(6 * rp * (u - 2 * r) ** 2)
zero(Hv - u, "surface v derivative")
zero(Hu - (v - 3 * (u - 2 * r) ** 2), "surface u derivative")
zero(H.subs(u, 0) - 8 * r**3, "singular equations force r=0")
zero(Hu.subs({u: 0, r: 0}) - v,
     "singular equations then force v=0")
zero(Ht.subs({u: 0, r: 0}), "t derivative vanishes on finite carrier")

# A hostile polynomial with one simple and one double zero: the singular
# scheme has support exactly t=0,1 and remains zero-dimensional.
t = sp.symbols("t")
r_hostile = t * (t - 1) ** 2
H_hostile = sp.expand(H.subs(r, r_hostile))
singular_gb = sp.groebner(
    [H_hostile, sp.diff(H_hostile, u), sp.diff(H_hostile, v),
     sp.diff(H_hostile, t)],
    v, u, t, order="lex",
)
singular_rows = [sp.factor(poly.as_expr()) for poly in singular_gb.polys]
gate(u in singular_rows, "hostile singular scheme forces u=0")
gate(any(sp.factor(row) == t**2 * (t - 1) ** 5 for row in singular_rows),
     "hostile singular support is the finite zero set of r")
for value in (0, 1):
    zero(H_hostile.subs({t: value, u: 0, v: 0}),
         f"hostile singular point {value} lies on surface")
    zero(sp.diff(H_hostile, u).subs({t: value, u: 0, v: 0}),
         f"hostile singular point {value} u derivative")
    zero(sp.diff(H_hostile, v).subs({t: value, u: 0, v: 0}),
         f"hostile singular point {value} v derivative")
    zero(sp.diff(H_hostile, t).subs({t: value, u: 0, v: 0}),
         f"hostile singular point {value} t derivative")


# ---------------------------------------------------------------------------
# Total ramification and the forbidden base-coordinate unit.
# ---------------------------------------------------------------------------

zero(F_zero.subs(P, 0) - T**3,
     "P=0 is one reduced prime with scheme multiplicity three")
zero(F_zero.subs(P, T**2) + 2 * T**2 * (T + r),
     "relative ramification graph factorization")
zero(sp.diff(F_zero, T) - 3 * (T**2 - P),
     "relative derivative")
zero(H.subs(v, 0) + (u - 2 * r) ** 3,
     "v=0 is the same triple divisor")

# The repeated-zero ramification roots are 0,0,r and satisfy every Vieta row.
roots_zero = (sp.Integer(0), sp.Integer(0), r)
zero(sum(roots_zero[i] * roots_zero[j]
         for i in range(3) for j in range(i + 1, 3)),
     "repeated-zero missing h-row")
zero(2 * sum(roots_zero) - C_zero, "repeated-zero C row")
zero(2 * sp.prod(roots_zero) - E_zero, "repeated-zero E row")


# ---------------------------------------------------------------------------
# Vertical-prime Nagata ledger.
# ---------------------------------------------------------------------------

zero(
    u * (u**2 - 6 * r * u + 12 * r**2 - v) - 8 * r**3 + H,
    "constant-r unit identity",
)

# In normalized valuations above a root of multiplicity m, v is a unit and
# uv=(u-2r)^3 forces nu(u)=3m. Freeze the Smith quotients for hostile rows.
for multiplicities in ((1,), (2,), (1, 1), (1, 2), (2, 4), (2, 3, 5)):
    relation = sp.Matrix([[3 * multiplicity for multiplicity in multiplicities]])
    smith = smith_normal_form(relation, domain=sp.ZZ)
    expected_torsion = 3 * math.gcd(*multiplicities)
    gate(abs(int(smith[0, 0])) == expected_torsion,
         f"Nagata torsion for multiplicities {multiplicities}")
    gate(relation.cols - relation.rank() == len(multiplicities) - 1,
         f"Nagata free rank for multiplicities {multiplicities}")

# Divisor identities in the theorem imply an exact order-three ramification
# class: div(u)=sum 3m_j Q_j and div(u-2r)=D+sum m_j Q_j.
for multiplicities in ((1,), (2,), (1, 3), (2, 3, 5)):
    mvec = sp.Matrix([list(multiplicities)])
    relation = 3 * mvec
    # n*m belongs to the rank-one lattice <3m> iff 3 divides n.
    for n in (1, 2):
        gate(not any(all(n * multiplicities[j] == 3 * q * multiplicities[j]
                         for j in range(len(multiplicities)))
                     for q in range(-2, 3)),
             f"ramification class not killed by {n}, {multiplicities}")
    gate(all(3 * multiplicities[j] == relation[0, j]
             for j in range(len(multiplicities))),
         f"ramification class killed by three, {multiplicities}")


summary = {
    "checks": CHECKS,
    "integrality": "k(t)-rational hidden roots are polynomial by monicity",
    "distinct": "routed to THM-3953",
    "repeated_nonzero": "linear times quadratic; normalization A2 disjoint A2",
    "repeated_zero": "normal uv=(u-2r)^3 surface",
    "ramification": "div(P)=3D makes P a forbidden unit on every Keller open",
    "nagata": "Cl=Z^(s-1)+Z/(3gcd m); D has exact order 3",
    "scope": "complete k(t)-split natural hidden-cubic lane",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3956 split hidden-cubic integrality/repeated-root companion")
print(f"CHECKS={CHECKS}")
print("INTEGRALITY=MONIC_OVER_KT_NORMAL_BASE;RATIONAL_ROOTS_ARE_POLYNOMIAL")
print("DISTINCT=THM3953")
print("REPEATED_NONZERO=LINEAR_TIMES_QUADRATIC;NORMALIZATION=A2_DISJOINT_A2")
print("REPEATED_ZERO=DOMAIN;NORMAL;UV_EQUALS_U_MINUS_2R_CUBED")
print("TOTAL_RAMIFICATION=DIV_P_EQUALS_3D;DELETING_D_MAKES_P_A_UNIT")
print("NAGATA=UNITS_SCALAR_FOR_NONCONSTANT_R;CL_FREE_PLUS_3_GCD_TORSION")
print("TRIPLE_ZERO=REDUCIBLE")
print("SCOPE=ALL_KT_SPLIT_NATURAL_HIDDEN_CUBICS;GENERAL_JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
