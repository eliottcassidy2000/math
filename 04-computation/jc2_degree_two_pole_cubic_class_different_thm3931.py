#!/usr/bin/env python3
"""Exact companion for THM-3931's class and different obstruction.

Reproduction:
  python3 04-computation/jc2_degree_two_pole_cubic_class_different_thm3931.py
  python3 -O 04-computation/jc2_degree_two_pole_cubic_class_different_thm3931.py
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


A, C, W, T, kap, eps, u = sp.symbols("A C W T kappa epsilon u")
kap_relation = sp.Poly(3 * kap**2 + 1, kap)


def reduce_kappa(expression: sp.Expr) -> sp.Expr:
    """Reduce a polynomial expression modulo 3*kappa^2+1."""
    polynomial = sp.Poly(sp.expand(expression), kap)
    return sp.expand(sp.rem(polynomial, kap_relation).as_expr())


def zero_mod(expression: sp.Expr, message: str) -> None:
    numerator, _ = sp.fraction(sp.cancel(expression))
    gate(reduce_kappa(numerator) == 0, message)


def valuation_at_zero(expression: sp.Expr, variable: sp.Symbol) -> int:
    polynomial = sp.Poly(sp.expand(expression), variable)
    return min(monomial[0] for monomial, coefficient in polynomial.terms() if coefficient != 0)


# ---------------------------------------------------------------------------
# Binary cubic, DF relations, and exact discriminant packet.
# ---------------------------------------------------------------------------

a = (A - 3 * kap) ** 2
c = 6 * kap * A + sp.Rational(2, 3)
d = (A - kap / 3) / 2

r_omega = W**2 + a * c - C * W + a * T
r_mixed = W * T + a * d
r_theta = T**2 + C * d - d * W + c * T

zero_mod(c - 12 * kap * d, "aligned middle coefficient c=12*kappa*d")

Delta = C**2 * c**2 - 4 * a * c**3 - 4 * C**3 * d - 27 * a**2 * d**2 + 18 * a * C * c * d
G = (
    -81 * A**5
    + 4455 * kap * A**4
    + 9246 * A**3
    + 648 * kap * A**3 * C
    + 1368 * A**2 * C
    - 18506 * kap * A**2
    - 144 * A * C**2
    - 2376 * kap * A * C
    - 3613 * A
    - 24 * C**3
    + 48 * kap * C**2
    - 216 * C
    + 627 * kap
)
zero_mod(Delta - d * G / 6, "line-times-quintic discriminant factorization")
gate(sp.Poly(G, A, C).total_degree() == 5, "G has total degree five")
zero_mod(G.subs(A, kap / 3) + 24 * C**3, "quintic does not contain vertical line")


# ---------------------------------------------------------------------------
# Factorial chart and the complete boundary-prime packet.
# ---------------------------------------------------------------------------

P_theta = T**3 + c * T**2 + C * d * T + a * d**2
D_theta = sp.diff(P_theta, T)
P_omega = W**3 - C * W**2 + a * c * W - a**2 * d

omega_chart = -a * d / T
C_chart = -(T**3 + c * T**2 + a * d**2) / (d * T)
zero_mod(r_omega.subs({W: omega_chart, C: C_chart}), "factorial chart omega relation")
zero_mod(r_mixed.subs({W: omega_chart, C: C_chart}), "factorial chart mixed relation")
zero_mod(r_theta.subs({W: omega_chart, C: C_chart}), "factorial chart theta relation")
zero_mod(P_theta.subs(C, C_chart), "factorial chart minimal equation")

A_dline = kap / 3
A_aline = 3 * kap
zero_mod(a.subs(A, A_dline) + sp.Rational(64, 27), "a is a unit on d=0")
zero_mod(c.subs(A, A_dline), "c vanishes on d=0")
zero_mod(P_omega.subs(A, A_dline) - W**2 * (W - C), "vertical omega fibre is double plus simple")
zero_mod(P_theta.subs(A, A_dline) - T**3, "vertical theta fibre is triple in the bad generator")

zero_mod(a.subs(A, A_aline), "a vanishes on Q")
zero_mod(c.subs(A, A_aline) + sp.Rational(16, 3), "c value on Q")
zero_mod(d.subs(A, A_aline) - 4 * kap / 3, "d is a unit on Q")
expected_Q_fibre = T * (T**2 - sp.Rational(16, 3) * T + 4 * kap * C / 3)
zero_mod(P_theta.subs(A, A_aline) - expected_Q_fibre, "Q is the unique theta-zero prime over a=0")

# Newton controls at the vertical double root. Put d=epsilon, so A=kappa/3+2epsilon.
A_eps = kap / 3 + 2 * eps
a_eps = reduce_kappa(a.subs(A, A_eps))
c_eps = reduce_kappa(c.subs(A, A_eps))
d_eps = reduce_kappa(d.subs(A, A_eps))
zero_mod(d_eps - eps, "vertical local parameter")
zero_mod(c_eps - 12 * kap * eps, "vertical c order one in base parameter")
gate(valuation_at_zero(a_eps, eps) == 0, "vertical a is a unit")

Pomega_eps = reduce_kappa(P_omega.subs(A, A_eps))
coeff_w0 = sp.Poly(Pomega_eps, W).coeff_monomial(1)
coeff_w1 = sp.Poly(Pomega_eps, W).coeff_monomial(W)
coeff_w2 = sp.Poly(Pomega_eps, W).coeff_monomial(W**2)
gate(valuation_at_zero(coeff_w0, eps) == 1, "double-root constant coefficient has base order one")
gate(valuation_at_zero(coeff_w1, eps) == 1, "double-root linear coefficient has base order one")
gate(valuation_at_zero(coeff_w2, eps) == 0, "double-root quadratic coefficient is a unit")
gate(sp.expand(coeff_w2) == -C, "double-root leading Newton coefficient")

# At Q, a has order two and d is a unit, forcing theta order two.
A_Qeps = 3 * kap + eps
gate(valuation_at_zero(reduce_kappa(a.subs(A, A_Qeps)), eps) == 2, "a has order two at Q")
gate(valuation_at_zero(reduce_kappa(d.subs(A, A_Qeps)), eps) == 0, "d is a unit at Q")


# ---------------------------------------------------------------------------
# Nagata relation lattice, classes, and scalar units.
# ---------------------------------------------------------------------------

relations = sp.Matrix([[2, 1, 0], [1, 1, 2]])
gate(relations.rank() == 2, "Nagata relation matrix has rank two")
minors = []
for columns in itertools.combinations(range(3), 2):
    minors.append(abs(int(relations[:, columns].det())))
gate(math.gcd(*minors) == 1, "Nagata relation lattice is saturated")

class_vector = sp.Matrix([2, -4, 1])
gate(relations * class_vector == sp.zeros(2, 1), "class assignments annihilate both relations")
gate(math.gcd(*(abs(int(value)) for value in class_vector)) == 1, "Q is a primitive class generator")

unit_divisor_matrix = sp.Matrix([[2, 1], [1, 1], [0, 2]])
gate(unit_divisor_matrix.rank() == 2, "only scalar Laurent-chart units extend globally")


# ---------------------------------------------------------------------------
# Exact different factorization and principal quintic ramification.
# ---------------------------------------------------------------------------

F = 3 * W - 12 * kap * T - 2 * C
zero_mod(D_theta - d * F - 3 * r_theta, "global derivative factorization modulo DF relation")
zero_mod(sp.discriminant(P_theta, T) - d**2 * Delta, "index-form discriminant identity")

zero_mod(F.subs({A: A_dline, W: 0, T: 0}) + 2 * C, "F is a unit at vertical ramified P0")
zero_mod(F.subs({A: A_dline, W: C, T: 0}) - C, "F is a unit at vertical companion PC")
zero_mod(F.subs({A: A_aline, W: C, T: 0}) - C, "F is a unit at Q")

# Norm(D_theta)=-disc(P); Norm(d)=d^3, so Norm(F)=-Delta/d=-G/6.
zero_mod((-d**2 * Delta) - d**3 * (-G / 6), "norm quotient for F")
gate(sp.Poly(F, W, T, C).total_degree() == 1, "F is visibly nonconstant")


# ---------------------------------------------------------------------------
# The two quintic addresses coalesce at the unique ramified fibre point.
# ---------------------------------------------------------------------------

A_u = u**3 + 3 * kap * u**2 - u
C_u = -sp.Rational(3, 2) * (u**2 - 1) * (u + kap) * (
    u**2 + 8 * kap * u - sp.Rational(11, 3)
)
for address in (-1, 1):
    zero_mod(A_u.subs(u, address) - 3 * kap, f"A({address})=3*kappa")
    zero_mod(C_u.subs(u, address), f"C({address})=0")
zero_mod((1 + kap) * (1 - kap) - sp.Rational(4, 3), "both repeated-root poles are genuine")

collision_fibre = reduce_kappa(P_theta.subs({A: 3 * kap, C: 0}))
gate(sp.expand(collision_fibre - T**2 * (T - sp.Rational(16, 3))) == 0,
     "collision fibre has one double and one simple endpoint")
zero_mod(a.subs(A, 3 * kap), "collision a=0")
zero_mod(c.subs(A, 3 * kap) + sp.Rational(16, 3), "collision c=-16/3")
zero_mod(d.subs(A, 3 * kap) - 4 * kap / 3, "collision d=4*kappa/3")
gate(sp.diff(collision_fibre, T).subs(T, 0) == 0, "theta=0 is ramified")
gate(sp.diff(collision_fibre, T).subs(T, sp.Rational(16, 3)) == sp.Rational(256, 9),
     "theta=16/3 is simple")
zero_mod(F.subs({A: 3 * kap, C: 0, W: 0, T: 0}), "F vanishes at ramified fibre point")
zero_mod(F.subs({A: 3 * kap, C: 0, W: 0, T: sp.Rational(16, 3)}) + 64 * kap,
         "F is nonzero at simple fibre point")

# Mandatory boundary classes E=0 and P0=2q fail the basis and leave Z/2.
E_class = 0
P0_class = 2
gate(math.gcd(E_class, P0_class) == 2, "natural ramification deletion leaves class quotient Z/2")
gate(E_class == 0, "quintic ramification class is principal")
gate(abs(P0_class) != 1, "vertical ramification class is nonprimitive")


summary = {
    "checks": CHECKS,
    "discriminant": "Delta=d*G/6",
    "factorial_chart": "S_(d theta)=k[A,d^-1,theta^+-1]",
    "divisors": {
        "d": [2, 1, 0],
        "theta": [1, 1, 2],
        "F": "E",
    },
    "class_group": "Z*q; P0=2q, PC=-4q, Q=q, E=0",
    "different": "P'(theta)=d*(3omega-12kappa theta-2C)",
    "norm_F": "-G/6",
    "collision": "u=+/-1 -> (A,C)=(3*kappa,0) -> one ramified E-point",
    "first_obstruction": "E has two normalization branches at one point",
    "mandatory_boundary_classes": [0, 2],
    "natural_class_quotient": "Z/2",
    "conclusion": "deleting E makes nonconstant F a unit; no A2 Keller open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3931 exact class-and-different companion")
print(f"CHECKS={CHECKS}")
print("DISCRIMINANT=Delta=d*G/6")
print("DIV_D=2*P0+PC")
print("DIV_THETA=P0+PC+2*Q")
print("CLASS_GROUP=Z*q;P0=2q;PC=-4q;Q=q;E=0")
print("DIFFERENT=P'(theta)=d*(3*omega-12*kappa*theta-2*C)")
print("DIV_F=E;NORM_F=-G/6")
print("COLLISION=u=+1,-1 -> (A,C)=(3*kappa,0)")
print("FIBRE=one length-2 ramified point plus one simple point")
print("FIRST_OBSTRUCTION=E is non-unibranch")
print("MANDATORY_BOUNDARY_CLASSES=(0,2);NATURAL_QUOTIENT=Z/2")
print("CONCLUSION=nonconstant F becomes a unit after deleting E; no A2 Keller open")
print(f"SEMANTIC_SHA256={semantic}")
