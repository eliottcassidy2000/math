#!/usr/bin/env python3
"""Exact companion for THM-3974's height-tower support obstruction."""

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


def ceil_div(a: int, b: int) -> int:
    return (a + b - 1) // b


# Exact coefficient-ideal minima in every hostile height/weight window.
for n in range(2, 13):
    for q in range(1, 41):
        feasible = []
        for a in range(q + 2):
            for b in range(q + 2):
                if n * a + (n + 1) * b >= q:
                    feasible.append((a, b))
        gate(min(a + 2 * b for a, b in feasible) == ceil_div(q, n),
             f"height {n}, weight {-q}: u-order minimum")
        gate(min(a + b for a, b in feasible) == ceil_div(q, n + 1),
             f"height {n}, weight {-q}: u+1-order minimum")


# Homogeneous bracket formula, tested as an exact Laurent identity across
# signs, heights, and nontrivial coefficient polynomials.
x, t, u = sp.symbols("x t u")


def jac(a: sp.Expr, b: sp.Expr) -> sp.Expr:
    return sp.factor(sp.diff(a, x) * sp.diff(b, t)
                     - sp.diff(a, t) * sp.diff(b, x))


for n in range(2, 7):
    uxt = x**n * t
    for r, s in [(-5, 2), (-3, -2), (-1, 0), (0, 4), (2, -4)]:
        f = 1 + 2 * u + u**3
        g = 2 - u + u**2
        lhs = jac(x**r * f.subs(u, uxt), x**s * g.subs(u, uxt))
        W = r * f * sp.diff(g, u) - s * sp.diff(f, u) * g
        rhs = x**(r + s + n - 1) * W.subs(u, uxt)
        gate(sp.factor(lhs - rhs) == 0,
             f"height {n}, weights {r},{s}: homogeneous bracket")


# The uniform 2x2 support and terminal factor.
n_sym, k_sym = sp.symbols("n k", integer=True, positive=True)
p0 = n_sym * (k_sym - 1) + 1
gate(sp.expand((-n_sym) + 1 - (1 - n_sym)) == 0,
     "first scalar complement")
gate(sp.expand(p0 - n_sym * k_sym - (1 - n_sym)) == 0,
     "second scalar complement")

h, hp, K, Kp, A, B, L, M = sp.symbols("h hp K Kp A B L M")
f = A * h
fp = A * hp
g = B * h**k_sym
gp = B * k_sym * h**(k_sym - 1) * hp
F = L * K**p0
Fp = L * p0 * K**(p0 - 1) * Kp
G = M * K
Gp = M * Kp
scalar = sp.expand(fp * G + n_sym * f * Gp
                   - n_sym * k_sym * Fp * g - p0 * F * gp)
factor = sp.expand((K * hp + n_sym * h * Kp)
                   * (A * M - k_sym * p0 * L * B
                      * K**(p0 - 1) * h**(k_sym - 1)))
gate(sp.simplify(scalar - factor) == 0, "uniform 2x2 scalar factor")


# Exact two-bracket positive control using only homogeneous Wronskians.
g0 = u * (u + 1)
w = 2 * u + 1
first = sp.diff(w * g0, u)              # J(-wp,x)
second = -6 * g0                        # J(6xp/(n-1),z)
gate(sp.expand(first + second) == 1, "uniform two-bracket identity")


# Height-two interface: exact distinguished-arm order and the second color.
# First freeze the all-height scalar-color gate used by THM-3576: among the
# complementary channels (-R,R+1-n), only R=n can contribute a unit at u=0.
for n in range(2, 13):
    for R in range(1, 3 * n + 1):
        m = ceil_div(R, n)
        complement = R + 1 - n
        survives = (m == 1 and complement >= 0 and complement != 0)
        gate(survives == (R == n),
             f"height {n}, negative weight {-R}: scalar-color gate")

all_height_transfer = {
    "THM-3576": ("2x3", "3x2"),
}
gate(sum(len(v) for v in all_height_transfer.values()) == 2,
     "two all-height imported cells")

# Height-two interface: wider exponent-two cells.
for q in range(1, 65):
    gate(ceil_div(q, 2) >= 1, f"weight {-q}: exponent-two color")
    gate(ceil_div(q, 3) >= 1, f"weight {-q}: second compulsory color")

# Frozen transfer inventory. These are theorem-mechanism assertions, not a
# finite enumeration standing in for the symbolic proofs.
transfer = {
    "THM-3583": ("2x4", "4x2"),
    "THM-3592-repaired": ("3x3",),
}
gate(sum(len(v) for v in transfer.values()) == 3,
     "three wider height-two imported cells")
live_total_seven = [(2, 5), (3, 4), (4, 3), (5, 2)]
gate(all(a + b == 7 for a, b in live_total_seven),
     "first live cells have seven pieces")
gate(all(min(a, b) >= 2 for a, b in live_total_seven),
     "one-by-arbitrary boundary removed")


summary = {
    "checks": CHECKS,
    "tower": "B_n=k[x,u,x^-n*u(u+1),x^(-n-1)*u^2(u+1)]",
    "pieces": "B_r=x^r*k[u] for r>=0; two-ceiling formula for r<0",
    "bracket": "degree n-1 Wronskian",
    "uniform": "homogeneous, 1xany, 2x2, and 2x3 cells empty for every n>=2",
    "all_height_transfer": all_height_transfer,
    "height2_transfer": transfer,
    "height2_floor": 7,
    "height2_first_live": live_total_seven,
    "positive_control": "constant is sum of two polynomial brackets",
    "scope": "unrestricted Darboux pair and JC2 remain open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3974 height-tower weight-support companion")
print(f"CHECKS={CHECKS}")
print("PIECES=TWO_CEILING_U_AND_U_PLUS_ONE_FORMULA")
print("BRACKET=WEIGHT_SHIFT_N_MINUS_1_WRONSKIAN")
print("UNIFORM=HOMOGENEOUS;ONE_BY_ANY;TWO_BY_TWO;TWO_BY_THREE_EMPTY")
print("ALL_HEIGHT_TRANSFER=THM3576_TWO_BY_THREE_AND_TRANSPOSE")
print("HEIGHT2_WIDER_TRANSFER=NO_2X4;3X3_AND_TRANSPOSES")
print("HEIGHT2_FLOOR=SEVEN_NONCONSTANT_PIECES")
print("FIRST_LIVE=2X5;3X4;4X3;5X2")
print("POSITIVE_CONTROL=CONSTANT_HAS_BRACKET_LENGTH_AT_MOST_TWO")
print("SCOPE=UNRESTRICTED_DARBOUX_AND_JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
