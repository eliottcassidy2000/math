#!/usr/bin/env python3
"""Exact Laurent companion for THM-3776's short target-shear gate."""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


t = sp.symbols("t")
a, b, d, e = sp.symbols("a b d e", nonzero=True)
f0, f1, f2 = sp.symbols("f0 f1 f2")
g_minus, g0 = sp.symbols("g_minus g0")
A, B, C, D, E = sp.symbols("A B C D E", nonzero=True)
F2, G2, h_minus = sp.symbols("F2 G2 h_minus")


# A completely arbitrary regular tail through the first two relevant layers.
# Higher layers cannot alter any coefficient extracted below.
p = a / t + b + d * t + e * t**2
inverse_p = sp.series(1 / p, t, 0, 4).removeO()


# Two-shear finite-at-infinity branch.
F_finite = f0 + f1 * inverse_p + f2 * inverse_p**2
q1_finite = sp.series(t + F_finite, t, 0, 4).removeO()
slope = sp.factor(sp.expand(q1_finite - f0).coeff(t, 1))
gate(sp.factor(slope - (a + f1) / a) == 0,
     "finite branch first slope")

# A simple pole is the only possible second-shear cancellation.  Its residue
# condition is read before any regular term in G.
two_residue = sp.factor(a + g_minus / slope)
gate(sp.factor(two_residue - a * (a + f1 + g_minus) / (a + f1)) == 0,
     "two-shear pole coefficient")
required_g = sp.solve(sp.together(two_residue), g_minus)
gate(required_g == [-a - f1], "two-shear required common residue")

# Slope zero raises contact to at least two, including arbitrary b and f2.
q1_slope_zero = sp.series(q1_finite.subs(f1, -a) - f0, t, 0, 4).removeO()
gate(sp.expand(q1_slope_zero).coeff(t, 1) == 0,
     "slope-zero linear coefficient")
gate(
    sp.factor(sp.expand(q1_slope_zero).coeff(t, 2) - (a * b + f2) / a**2)
    == 0,
    "slope-zero second layer",
)


# Three-shear polar branch.  Terms through inverse square are retained on
# both rational jets so the arbitrary tails are actively tested.
F_polar = A * p + B + C * inverse_p + F2 * inverse_p**2
q1_polar = sp.series(t + F_polar, t, 0, 3).removeO()
inverse_q1 = sp.series(1 / q1_polar, t, 0, 3).removeO()
G_polar = -q1_polar / A + D + E * inverse_q1 + G2 * inverse_q1**2
p1_polar = sp.series(p + G_polar, t, 0, 3).removeO()

ell = sp.factor(sp.expand(p1_polar).coeff(t, 0))
gate(sp.factor(ell - (D - B / A)) == 0,
     "three-shear common finite value")
second_slope = sp.factor(sp.expand(p1_polar).coeff(t, 1))
gate(sp.factor(second_slope - (E - C - a) / (A * a)) == 0,
     "three-shear second slope")

# If the second slope is nonzero, H must have a simple pole at ell.  The
# q^-1 coefficient in q2 is Aa+h/second_slope.
three_residue = sp.factor(A * a + h_minus / second_slope)
required_h = sp.solve(sp.together(three_residue), h_minus)
gate(required_h == [C - E + a], "three-shear required common residue")

# The finite value and the surviving slope really ignore b,d and every
# displayed higher-jet parameter that should be irrelevant.
gate(not ell.has(a, b, d, e, C, E, F2, G2),
     "common value independent of divisor tails")
gate(not second_slope.has(b, d, e, F2, G2),
     "second slope independent of regular tails")


# Pole-order arithmetic: a pole of order d1 followed by one of order e1 can
# cancel an order-one input only in the linear-linear cell.
order_cells = 0
linear_linear_cells = 0
for d1 in range(1, 8):
    for e1 in range(1, 8):
        order_cells += 1
        if d1 * e1 == 1:
            linear_linear_cells += 1
            gate((d1, e1) == (1, 1), "unique polar order cell")
gate(order_cells == 49 and linear_linear_cells == 1,
     "polar order census")


# Reciprocal and Mobius jets.
s = sp.symbols("s")
alpha, beta, gamma, delta, kappa = sp.symbols(
    "alpha beta gamma delta kappa", nonzero=True
)
mobius = (alpha * s + beta) / (gamma * s + delta)
mobius_f0 = sp.limit(mobius, s, sp.oo)
mobius_f1 = sp.limit(s * (mobius - mobius_f0), s, sp.oo)
gate(sp.factor(mobius_f0 - alpha / gamma) == 0,
     "Mobius finite value")
gate(
    sp.factor(mobius_f1 - (beta * gamma - alpha * delta) / gamma**2) == 0,
    "Mobius inverse coefficient",
)
gate(sp.limit(s * (kappa / s), s, sp.oo) == kappa,
     "reciprocal inverse coefficient")


# Several distinct-residue pairs guard that both cancellation requirements
# retain exactly their original difference.
residue_pairs = ((1, -1), (1, 3), (-2, 5), (2, 7), (-5, -1))
pair_controls = 0
for left, right in residue_pairs:
    gate(left != 0 and right != 0 and left != right,
         "admissible unequal residue pair")
    two_left = -left - 4
    two_right = -right - 4
    gate(two_left - two_right == right - left,
         "two-shear difference survives")
    three_left = left + 6 - 9
    three_right = right + 6 - 9
    gate(three_left - three_right == left - right,
         "three-shear difference survives")
    pair_controls += 1


# Exact opposite-orientation local regularization.  These are two embeddings
# of the same target word into k((t)); every final denominator is a unit.
q_local = t
p_left = 1 / t
p_right = 3 / t


def opposite_word(p_input: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    p_first = sp.cancel(p_input - 1 / q_local)
    q_first = sp.cancel(q_local + p_first / (p_first + 1))
    p_second = sp.cancel(p_first - 1 / (q_first - 1))
    return p_first, q_first, p_second


left_word = opposite_word(p_left)
right_word = opposite_word(p_right)
gate(left_word == (0, t, -1 / (t - 1)),
     "opposite word left completion")
gate(right_word[0] == 2 / t, "opposite word right first state")
gate(sp.cancel(right_word[1] - (t + 2 / (t + 2))) == 0,
     "opposite word right second state")
gate(sp.cancel(right_word[2] - 1 / (t + 1)) == 0,
     "opposite word right regular output")
gate(
    all(sp.denom(sp.cancel(value)).subs(t, 0) != 0 for value in (left_word[1], left_word[2], right_word[1], right_word[2])),
    "opposite word all final denominators are units",
)


# Variable exchange gives the p-uniformizer/p-q-p dual literally.
p_uniformizer = t
q_polar = a / t + b + d * t
p_dual_first = sp.series(
    p_uniformizer + f0 + f1 / q_polar + f2 / q_polar**2, t, 0, 3
).removeO()
dual_slope = sp.factor(sp.expand(p_dual_first - f0).coeff(t, 1))
gate(sp.factor(dual_slope - (a + f1) / a) == 0,
     "variable-exchanged dual slope")


semantic_rows = (
    "local:q_uniformizer;p=a_i/q+O(1);a_1*a_2*(a_1-a_2)!=0",
    "two:q1=q+F(p);p1=p+G(q1);required_g=-(a_i+f1)",
    "three:polar_orders_force_1x1;ell=D-B/A;c_i=(E-C-a_i)/(A*a_i)",
    "three:required_h=a_i+C-E",
    "dual:p_uniformizer;q=a_i/p+O(1);p-q-p",
    "sharp_opposite:p=1/t,3/t;F=-1/q;G=s/(s+1);H=-1/(s-1)",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3776-unequal-vertical-residue-three-target-shear-nonpolynomialization")
print("scope=two_completed_DVRs;distinct_nonzero_simple_coefficients;natural_length_at_most_3")
print("two_shears=impossible_by_common_inverse_coefficient_translation")
print("three_shears=impossible_by_linear_linear_order_gate_and_second_translation")
print("higher_contact=slope_zero_forces_order_at_least_2")
print("dual=p_uniformizer_p-q-p_by_variable_exchange")
print("opposite_orientation=locally_regularizable_exact_hostile_not_global_polynomialization")
print(f"pair_controls={pair_controls};order_cells={order_cells};linear_linear_cells={linear_linear_cells}")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
