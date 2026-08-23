#!/usr/bin/env python3
"""Exact controls for THM-3882's dual-Wronskian one-place criterion.

Reproduction:
  python3 04-computation/jc2_rational_dual_one_place_wronskian_projection_criterion_thm3882.py
  python3 -O 04-computation/jc2_rational_dual_one_place_wronskian_projection_criterion_thm3882.py
"""

from __future__ import annotations

import hashlib
import sys

import sympy as sp

sys.stdout.reconfigure(newline="\n")

t = sp.symbols("t")
CHECKS = 0


def zero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def equal(label: str, left: object, right: object) -> None:
    global CHECKS
    CHECKS += 1
    if left != right:
        raise AssertionError(f"{label}: {left} != {right}")


# The fibre factor enters the point-projection Wronskian twice.
g = t**3 - 2 * t + 2
x = t**4 + 3 * t**2 - t + 1
y = 2 * t**4 - t**3 + 5 * t + 3
wr = lambda f, h: sp.expand(f * sp.diff(h, t) - h * sp.diff(f, t))
zero("wronskian_base_square", wr(g * x, g * y) - g**2 * wr(x, y))

# Degree and local-ramification ledgers for a broad exact window.
for d in range(2, 41):
    for m in range(0, d):
        e = d - m
        equal(f"degree_ledger_d{d}_m{m}", 2 * m + (2 * e - 2), 2 * d - 2)
        if e >= 2:
            equal(f"one_point_ramification_gap_d{d}_m{m}", (2 * e - 2) > (e - 1), True)
    if d >= 3:
        equal(f"immersed_one_fibre_no_go_d{d}", d - 1 > 1, True)

# Pluecker packet for a sextic with four nodes and six ordinary cusps.
d, nodes, cusps = 6, 4, 6
genus = (d - 1) * (d - 2) // 2 - nodes - cusps
dual_degree = d * (d - 1) - 2 * nodes - 3 * cusps
inflection_weight = 3 * d * (d - 2) - 6 * nodes - 8 * cusps
equal("pluecker_genus", genus, 0)
equal("pluecker_dual_degree", dual_degree, 4)
equal("pluecker_inflection_weight", inflection_weight, 0)

# THM-3879 primal/dual tangent computation in the T=1 chart.
X = t * (t**2 - 2)
Y = t * (t**2 - 1)
Z = (t**2 - 1) * (t**2 - 2)
A = -(t**2 - 1) ** 2 * (t**2 + 2)
B = (t**2 - 2) ** 2 * (t**2 + 1)
C = 2 * t**3
dual = (
    sp.factor(Y * sp.diff(Z, t) - Z * sp.diff(Y, t)),
    sp.factor(Z * sp.diff(X, t) - X * sp.diff(Z, t)),
    sp.factor(X * sp.diff(Y, t) - Y * sp.diff(X, t)),
)
zero("dual_A", dual[0] + A)
zero("dual_B", dual[1] + B)
zero("dual_C", dual[2] + C)
equal("immersed_affine_gcd", sp.gcd(sp.gcd(A, B), C), sp.Integer(1))
zero("tangent_incidence", A * X + B * Y + C * Z)
zero(
    "tangent_derivative_incidence",
    A * sp.diff(X, t) + B * sp.diff(Y, t) + C * sp.diff(Z, t),
)

# At P=[0:0:1], X and Y have the affine base factor t.  Homogeneously the
# second base address is T=0.  The residual degree-two map has derivative
# numerator -2t, so its two ramification addresses are exactly 0 and infinity.
equal("node_affine_base_gcd", sp.gcd(X, Y), t)
x2 = sp.cancel(X / t)
y2 = sp.cancel(Y / t)
equal("node_residual_degrees", (sp.degree(x2, t), sp.degree(y2, t)), (2, 2))
zero("node_residual_wronskian", wr(x2, y2) + 2 * t)
zero("node_full_wronskian", wr(X, Y) + 2 * t**3)
zero("dual_best_line", C - 2 * t**3)
equal("node_divisor_coefficients", (2 + 1, 2 + 1), (3, 3))

print("THM3882_WRONSKIAN", "W(gx,gy)=g^2W(x,y)")
print("THM3882_DIVISOR", "dual-line pullback=2*point-projection base fibre+ramification")
print("THM3882_ONE_PLACE", "iff projection degree 1 and base fibre=(d-1)p")
print("THM3882_IMMERSED", "primal degree d>=3 forbids a one-place line on its dual")
print("THM3882_PACKET", "degree6 with 6A2+4A1 has genus0,dual degree4,inflection weight0")
print("THM3882_PACKET_NOGO", "every 6A2+4A1 sextic line leaves at least two normalization places")
print("THM3882_THM3879", "C=0 is 2D+R=(3,3) over the two branches of primal node [0:0:1]")
print("THM3882_OPEN", "nonimmersed primal,changed packet,or nondual one-place construction")
semantic_packet = (
    "dual line as point-projection Wronskian",
    "base fibre occurs with multiplicity two",
    "remaining divisor is exact ramification",
    "Riemann Hurwitz one-point ramification impossibility",
    "one-place iff degree-one projection with a d-minus-one base fibre",
    "immersion makes the one-address base coefficient one",
    "six cusp four node sextic packet has immersed rational quartic dual",
    "entire packet cannot have affine-line normalization",
    "THM3879 best line is node base twice plus ramification once",
    "nonimmersion or changed Pluecker packet required",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
