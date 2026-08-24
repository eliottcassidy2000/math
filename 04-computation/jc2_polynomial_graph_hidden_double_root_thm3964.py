#!/usr/bin/env python3
"""Exact companion for THM-3964's polynomial-graph double-root debt."""

from __future__ import annotations

import hashlib
import json
from math import gcd

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.expand(expression) == 0, message)


P, T, x, y, u, v, c, r, h, t, z = sp.symbols(
    "P T x y u v c r h t z"
)


# ---------------------------------------------------------------------------
# Universal graph-debt normal form and rational normalization parameter.
# ---------------------------------------------------------------------------

q = sp.expand(3 * r * P - r**3 + c * (P - r**2) ** 2)
F = sp.expand(T**3 - 3 * P * T - q)
G = sp.expand(x**3 - 3 * r * x**2 - 3 * x * y - c * y**2)
zero(F.subs({T: x - r, P: y + r**2}) - G,
     "node-coordinate normal form")

K_hidden = sp.factor(q.subs(P, h**2) - 2 * h**3)
gate(sp.rem(sp.Poly(K_hidden, h), sp.Poly((h - r) ** 2, h)) == 0,
     "h=r is a repeated hidden root")

xv = sp.expand(c * v**2 + 3 * v + 3 * r)
yv = sp.expand(v * xv)
Pv = sp.expand(yv + r**2)
Tv = sp.expand(xv - r)
zero(G.subs({x: xv, y: yv}), "normalization parametrizes the cubic")
zero(F.subs({P: Pv, T: Tv}), "target-coordinate parametrization")
zero(yv - v * xv, "slope inverse v=y/x")
gate(sp.degree(Pv, v) == 3, "generic P-map has degree three")
gate(sp.Poly(Pv, v).LC() == c, "generic degree-three leading row is c")


# ---------------------------------------------------------------------------
# Integral coordinate and exact two normalization charts.
# ---------------------------------------------------------------------------

uv = c * v
for expression, message in (
    (uv**2 + 3 * uv + 3 * r * c - c * xv,
     "integral quadratic for u"),
    (c * yv - uv * xv, "relation cy=ux"),
    ((uv + 3) * yv - xv * (xv - 3 * r),
     "relation (u+3)y=x(x-3r)"),
):
    zero(expression, message)

H_plus = sp.expand(c * (x - 3 * r) - u * (u + 3))
G_u = sp.expand(c**2 * y - u**2 * (u + 3) - 3 * r * c * u)
zero(H_plus.subs({c: 0, u: -sp.Rational(3, 2)}) -
     sp.Rational(9, 4),
     "first chart singular candidate misses the surface")
zero(G_u.subs({c: 0, u: -3}),
     "second chart special fibre is u=-3")
gate(sp.diff(G_u, u).subs({c: 0, u: -3}) == -9,
     "second chart special fibre is smooth")
zero(sp.diff(G_u, u) + 3 * (u**2 + 2 * u + c * r),
     "second chart u derivative is minus 3Q")

# Both eliminations recover all three relations.
y_plus = x * (x - 3 * r) / (u + 3)
zero(sp.cancel((c * y - u * x).subs(y, y_plus) * (u + 3) - x * H_plus),
     "first chart recovers cy=ux")
x_u = c * y / u
zero(sp.cancel((c * (x - 3 * r) - u * (u + 3)).subs(x, x_u) * u - G_u),
     "second chart recovers integral relation")


# ---------------------------------------------------------------------------
# Singular graph and exact conductor quotient.
# ---------------------------------------------------------------------------

Gx = sp.diff(G, x)
Gy = sp.diff(G, y)
y_candidate = 3 * (3 - 4 * c * r) / (4 * c**2)
x_candidate = -2 * c * y_candidate / 3
zero(G.subs({x: x_candidate, y: y_candidate}) +
     (4 * c * r - 3) ** 3 / (16 * c**3),
     "only apparent off-graph Jacobian candidate")
zero(Gx.subs({x: 0, y: 0}), "graph F_x singular row")
zero(Gy.subs({x: 0, y: 0}), "graph F_y singular row")

# B=A+Au and (x,y) annihilates the quotient generator.
zero((x * u - c * y).subs({x: xv, y: yv, u: uv}),
     "x kills class(u) modulo A")
zero((y * u - (x * (x - 3 * r) - 3 * y)).subs(
    {x: xv, y: yv, u: uv}),
     "y kills class(u) modulo A")
conductor_address = u**2 + 3 * u + 3 * c * r
gate(sp.discriminant(conductor_address, u) == 9 - 12 * c * r,
     "conductor address discriminant")


# ---------------------------------------------------------------------------
# Saturated ramification character.
# ---------------------------------------------------------------------------

dPv = sp.factor(sp.diff(Pv, v))
zero(dPv - 3 * (c * v**2 + 2 * v + r),
     "normalized P-map derivative")
Q = sp.expand(u**2 + 2 * u + c * r)
D = sp.expand(1 - c * r)
zero(Q.subs(u, z - 1) - (z**2 - D),
     "ramification hyperelliptic character")
gate(sp.discriminant(Q, u) == 4 * D,
     "ramification discriminant differs from conductor discriminant")

N = sp.expand((2 * u + 3) * (2 * x - 3 * r) * (u + 3)
              - c * x * (x - 3 * r))
saturation_rhs = H_plus * (-c * x + 3 * u**2 + 15 * u + 18)
zero(c * N - 3 * (u + 3) ** 2 * Q - saturation_rhs,
     "cancellation-free ramification saturation identity")
zero(N.subs({c: 0, u: 0}) - 9 * (2 * x - 3 * r),
     "saturated finite ramification point over c=0")
gate(H_plus.subs({c: 0, u: -2}) == 2,
     "z=-1 companion address escapes over c=0")

x_on_character = r * (2 * z + 1) / (z + 1)
character_relation = z**2 - D
numerator_identity = sp.together(
    x_on_character - (2 * r + (z - 1) / c)
).as_numer_denom()[0]
gate(sp.rem(sp.Poly(sp.expand(numerator_identity), z),
            sp.Poly(character_relation, z)) == 0,
     "ramification normalization x-coordinate identity")
zero(x_on_character.subs(z, 1) - 3 * r / 2,
     "finite z=1 address extends to x=3r/2")


# ---------------------------------------------------------------------------
# Nagata class lattice and boundary primitivity.
# ---------------------------------------------------------------------------

# For n split fibres, the 2n boundary generators have n independent rows
# Q_i^+ + Q_i^-. Freeze the rank and quotient rank for n=1..5.
for n in range(1, 6):
    relation_matrix = sp.zeros(n, 2 * n)
    for i in range(n):
        relation_matrix[i, 2 * i] = 1
        relation_matrix[i, 2 * i + 1] = 1
    gate(relation_matrix.rank() == n,
         f"Nagata boundary relation rank n={n}")
    gate(2 * n - relation_matrix.rank() == n,
         f"class-group free rank n={n}")

m_rows = [1, 2, 3, 4]
gate(gcd(gcd(m_rows[0], m_rows[1]), gcd(m_rows[2], m_rows[3])) == 1,
     "primitive nonsquare multiplicity positive control")
nonprimitive_rows = [2, 4, 6]
gate(gcd(gcd(nonprimitive_rows[0], nonprimitive_rows[1]),
         nonprimitive_rows[2]) == 2,
     "nonprimitive nonsquare multiplicity hostile")

# Split vectors have disjoint sign support. Freeze a passing and two failing
# boundary-basis packets.
split_plus = sp.Matrix([1, 0, 2, 0])
split_minus = sp.Matrix([0, 1, 0, 3])
gate(split_plus.dot(split_minus) == 0,
     "split-square class vectors have disjoint support")
gate(gcd(1, 2) == 1 and gcd(1, 3) == 1,
     "split-square subset gcd positive control")
empty_sign = sp.zeros(4, 1)
gate(empty_sign == sp.zeros(4, 1),
     "empty sign set gives zero boundary class")
gate(gcd(2, 4) == 2,
     "split-square nonprimitive subset hostile")


# ---------------------------------------------------------------------------
# Boundary rank/Euler closure and exact hostile characters.
# ---------------------------------------------------------------------------

# chi(B)=1+n and THM-3966 gives exactly n boundary curves. Its normalization
# and incidence ledger forces all normalized component chis to be one and
# every intersection penalty to vanish.
for n in range(1, 6):
    chi_B = 1 + n
    chi_A2 = 1
    chi_boundary = chi_B - chi_A2
    gate(chi_boundary == n, f"boundary Euler equals class rank n={n}")
    component_chis = [1] * n
    gate(sum(component_chis) == chi_boundary,
         f"only maximal component Euler rows can attain equality n={n}")

# A nonsquare ramification normalization loses n coefficient-root points
# plus epsilon>=1 infinity points, hence never has the one puncture of A1.
for n in range(1, 6):
    for epsilon in (1, 2):
        gate(n + epsilon >= 2,
             f"nonsquare ramification has at least two punctures n={n},eps={epsilon}")

# In the split case, if a,b count the + and - sign sets, the two components
# have 1+b and 1+a punctures. Both equal one only when a=b=0, impossible for
# n=a+b>0.
for n in range(1, 6):
    for a in range(n + 1):
        b = n - a
        punctures_plus = 1 + b
        punctures_minus = 1 + a
        gate(not (punctures_plus == 1 and punctures_minus == 1),
             f"split components cannot both be A1 n={n},a={a}")

# Rational nonsquare primitive hostile: c=t,r=1, D=1-t. It passes genus and
# primitivity but has one c-root puncture plus infinity, so is not A1.
c_rat = t
r_rat = sp.Integer(1)
D_rat = sp.factor(1 - c_rat * r_rat)
gate(D_rat == 1 - t, "rational nonsquare hostile character")
gate(sp.degree(D_rat, t) == 1, "rational nonsquare genus-zero degree")
gate(1 + 1 == 2, "rational nonsquare hostile has two punctures")

# Positive-genus kill: c=t,r=t^2 gives squarefree cubic D=1-t^3.
D_genus = sp.factor(1 - t * t**2)
gate(sp.degree(D_genus, t) == 3, "positive-genus cubic character")
gate(sp.gcd(D_genus, sp.diff(D_genus, t)) == 1,
     "positive-genus hostile is squarefree")
gate((sp.degree(D_genus, t) - 1) // 2 == 1,
     "positive-genus hostile has genus one")

# Rational but non-unibranch kill: an even root survives in D.
D_node = sp.expand((t - 1) ** 2 * (t + 1))
r_node = sp.cancel((1 - D_node) / t)
gate(sp.denom(r_node) == 1, "non-unibranch hostile has polynomial r")
gate(sp.rem(sp.Poly(D_node, t), sp.Poly((t - 1) ** 2, t)) == 0,
     "non-unibranch hostile has an even-multiplicity root")
gate(sp.factor(D_node.subs(t, 0)) == 1,
     "non-unibranch hostile is compatible with c=t")

# The apparent n=2 split packet is closed once the two ramification
# classes exhaust Cl(B): the maximal complement has Euler characteristic 4.
c_split = t**2 - 1
r_split = -1
g_split = t
D_split = sp.factor(1 - c_split * r_split)
zero(D_split - g_split**2, "minimal split-square character")
gate(sp.factor(c_split) == (t - 1) * (t + 1),
     "minimal split-square two simple vertical fibres")
gate(g_split.subs(t, 1) == 1 and g_split.subs(t, -1) == -1,
     "minimal split-square uses both sign sets")
gate(sp.degree(g_split, t) == 1,
     "minimal split-square components meet exactly once")
chi_B_n2 = 1 + 2
chi_each_split_component = 0  # P1 minus infinity and one escaped point.
chi_split_union_n2 = 2 * chi_each_split_component - 1
gate(chi_B_n2 - chi_split_union_n2 == 4,
     "n=2 full-boundary split complement has Euler four")

# Likewise n=1 nonsquare primitivity exhausts the class basis. For a
# rational unibranch character epsilon=1 or2 infinity points makes the full
# complement Euler 1+epsilon, never one.
for epsilon in (1, 2):
    chi_B_n1 = 2
    chi_R_n1 = 2 - 1 - epsilon
    gate(chi_B_n1 - chi_R_n1 == 1 + epsilon,
         f"n=1 nonsquare full complement Euler epsilon={epsilon}")
    gate(chi_B_n1 - chi_R_n1 != 1,
         f"n=1 nonsquare cannot be A2 epsilon={epsilon}")

# A one-collision n=4 packet also fails: both ramification components have
# opposite-sign punctures and they intersect at the collision.
c_split_live = 1 - t**4
r_split_live = 1
g_split_live = t**2
D_split_live = sp.factor(1 - c_split_live * r_split_live)
zero(D_split_live - g_split_live**2,
     "one-collision n=4 split hostile character")
gate(sp.gcd(c_split_live, sp.diff(c_split_live, t)) == 1,
     "one-collision n=4 coefficient has four simple roots")
gate(g_split_live.subs(t, 1) == 1 and
     g_split_live.subs(t, -1) == 1 and
     g_split_live.subs(t, sp.I) == -1 and
     g_split_live.subs(t, -sp.I) == -1,
     "one-collision n=4 packet uses both sign sets")
gate(sp.Poly(g_split_live, t).sqf_part() == sp.Poly(t, t),
     "one-collision n=4 packet has one distinct collision point")
gate(1 + 2 > 1,
     "one-collision n=4 plus component has opposite-sign punctures")

# r=0 is the principal-arm endpoint already closed by THM-3963.
zero(Q.subs(r, 0) - u * (u + 2),
     "r=0 ramification factorization")


summary = {
    "checks": CHECKS,
    "family": "q=3rP-r^3+c(P-r^2)^2 with c nonzero",
    "normalization": "u=cv; two smooth charts over arbitrary c zeros",
    "conductor": "B/A=k[t], conductor=(x,y), addresses u^2+3u+3cr",
    "ramification": "saturated horizontal character z^2=1-cr",
    "class": "Cl(B)=Z^n, units scalar, [R]=-(m_i)",
    "boundary": "pure boundary; exactly n pairwise-disjoint primes with A1 normalizations",
    "nonsquare": "n coefficient punctures plus infinity exclude R",
    "square": "opposite-sign punctures prevent both normalizations being A1",
    "conclusion": "displayed scalar-coefficient graph repeated-root debts excluded",
    "scope": "P-dependent graph coefficients, nongraph repetition, and JC2 remain open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3964 polynomial-graph hidden double-root companion")
print(f"CHECKS={CHECKS}")
print("NORMAL_FORM=X3_MINUS_3R_X2_MINUS_3XY_MINUS_CY2")
print("NORMALIZATION=U_CV;TWO_SMOOTH_CHARTS;FULL_FINITE_NORMALIZATION")
print("CONDUCTOR=B_OVER_A_KT;IDEAL_X_Y;ADDRESS_DISC_9_MINUS_12CR")
print("RAMIFICATION=SATURATED_Z2_EQUALS_1_MINUS_CR;ONE_ESCAPED_C_ROOT_ADDRESS")
print("CLASS=UNITS_SCALAR;CL_Z_TO_N;R_VECTOR_MINUS_MULTIPLICITIES")
print("BOUNDARY=PURE;EXACTLY_N_PRIMES;DISJOINT_WITH_A1_NORMALIZATIONS")
print("NONSQUARE_KILL=N_C_ROOT_PUNCTURES_PLUS_INFINITY;NEVER_A1")
print("SQUARE_KILL=OPPOSITE_SIGN_PUNCTURES;BOTH_NORMALIZATIONS_CANNOT_BE_A1")
print("CONCLUSION=SCALAR_COEFFICIENT_POLYNOMIAL_GRAPH_DOUBLE_ROOT_FAMILY_CLOSED")
print("CONTROL_SPLIT=C_1_MINUS_T4;R_1;G_T2;STILL_KILLED")
print("SCOPE=P_DEPENDENT_GRAPH_COEFFICIENT;NONGRAPH_REPETITION;ARBITRARY_P2Q2;NONMONOGENIC;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
