#!/usr/bin/env python3
"""Exact controls for THM-3995's scale-two parity-hole tariff."""

from fractions import Fraction

import sympy as sp


def require(label, condition):
    if not condition:
        raise RuntimeError(f"FAIL  {label}")
    print(f"PASS  {label}")


# THM-3878's exact obstruction for the scale-two (p,q)=(1,9) row.
boundaries = (
    Fraction(2, 21),
    Fraction(8, 63),
    Fraction(55, 63),
    Fraction(19, 21),
)
component_lengths = (
    boundaries[1] - boundaries[0],
    boundaries[3] - boundaries[2],
)
require("two obstruction components each have length 2/63", component_lengths == (Fraction(2, 63),) * 2)
require("obstruction measure is 4/63", sum(component_lengths) == Fraction(4, 63))

# Entering a safe owner interval uses 14a+1; leaving uses 14a-1.  After
# clearing denominators at the four oriented target boundaries, the event
# numerators are the following.  For odd t they are always odd, independently
# of a,u,m, so equality is impossible and the numerator has absolute value >=1.
def parity_rows(t, a, u, m):
    return (
        3 * t * (14 * a + 1) - 4 * u - 42 * u * m,
        9 * t * (14 * a - 1) - 16 * u - 126 * u * m,
        9 * t * (14 * a + 1) - 110 * u - 126 * u * m,
        3 * t * (14 * a - 1) - 38 * u - 42 * u * m,
    )


parity_checks = 0
for t_residue in (1, 3):
    for a_residue in (0, 1):
        for u_residue in (0, 1):
            for m_residue in (0, 1):
                values = parity_rows(t_residue, a_residue, u_residue, m_residue)
                require_label = all(value % 2 == 1 for value in values)
                if not require_label:
                    raise RuntimeError((t_residue, a_residue, u_residue, m_residue, values))
                parity_checks += 4
print(f"PASS  all four cleared event numerators are odd ({parity_checks} residue checks)")

hole_coefficients = (
    Fraction(1, 42),
    Fraction(1, 126),
    Fraction(1, 126),
    Fraction(1, 42),
)
require("four hole coefficients total 4/63", sum(hole_coefficients) == Fraction(4, 63))

t = sp.symbols("t", positive=True, integer=True)
s_t = sp.Rational(4, 63) - sp.Rational(4, 63) / t
require("support cap simplifies", sp.factor(s_t - 4 * (t - 1) / (63 * t)) == 0)
require(
    "real variance coefficient simplifies",
    sp.factor((1 - s_t) / s_t - (59 * t + 4) / (4 * (t - 1))) == 0,
)

# Independent finite dynamic-programming controls for the integer optimization:
# among q equally weighted cells, total mass M, and at most B positive cells,
# the minimum sum of squares is obtained by spreading M over all B allowed
# cells at the adjacent integers floor(M/B), ceil(M/B).
integer_checks = 0
for q in range(2, 10):
    max_total = 4 * q
    dp = {(0, 0): 0}
    for _ in range(q):
        nxt = {}
        for (positive, total), square_sum in dp.items():
            for value in range(max_total - total + 1):
                key = (positive + (value > 0), total + value)
                candidate = square_sum + value * value
                if candidate < nxt.get(key, 10**9):
                    nxt[key] = candidate
        dp = nxt
    for B in range(1, q + 1):
        for M in range(0, 4 * B + 1):
            actual = min(
                square_sum
                for (positive, total), square_sum in dp.items()
                if total == M and positive <= B
            )
            k, remainder = divmod(M, B)
            predicted = (B - remainder) * k * k + remainder * (k + 1) ** 2
            if actual != predicted:
                raise RuntimeError((q, B, M, actual, predicted))
            integer_checks += 1
print(f"PASS  integer support-envelope dynamic controls ({integer_checks} cases)")

# A hostile in which the integral correction is strictly positive.
B_hostile, M_hostile = 5, 7
k_hostile, rem_hostile = divmod(M_hostile, B_hostile)
require("integer correction hostile is positive", 0 < rem_hostile < B_hostile and k_hostile == 1)

print("RESULT obstruction_measure=4/63")
print("RESULT parity_hole_measure=4/(63t)")
print("RESULT support_cap=4(t-1)/(63t)")
print("RESULT sufficient_gate=t*mu(G)/r>sqrt(4(t-1)/(3(59t+4)))")
print("ALL EXACT CHECKS PASSED")
