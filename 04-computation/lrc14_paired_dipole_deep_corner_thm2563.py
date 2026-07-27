#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2563.

Checks the paired translated-tooth capacities, the two-zero-plane Fourier
corner, the missing-left-residue hostile, and the exact THM-2561 physical
interval table.  All arithmetic is integer or Fraction arithmetic.
"""

from collections import defaultdict
from fractions import Fraction


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def floor_q(value):
    return value.numerator // value.denominator


def fractional_part(value):
    return value - floor_q(value)


def circle_norm(value):
    residue = fractional_part(value)
    return min(residue, 1 - residue)


def danger(value, level=1):
    return circle_norm(value) < Fraction(level, 14)


print("== paired translated-tooth capacities ==")
count_values = {1: set(), 2: set()}
for level in (1, 2):
    for cell in range(182):
        y = Fraction(2 * cell + 1, 364)
        count = sum(danger(y - Fraction(q, P), level) for q in range(P))
        expected = 2 * level - danger(P * y, level)
        require(count == expected, "translated-tooth identity failed")
        count_values[level].add(count)
require(count_values == {1: {1, 2}, 2: {3, 4}}, "sharp tooth counts changed")

pair_floor = {level: P - 2 - 2 * level for level in (1, 2)}
require(pair_floor == {1: 9, 2: 7}, "paired safe capacity changed")
role_products = {
    (1, 1): pair_floor[1] ** 2,
    (1, 2): pair_floor[1] * pair_floor[2],
    (2, 1): pair_floor[2] * pair_floor[1],
}
require(role_products == {(1, 1): 81, (1, 2): 63, (2, 1): 63},
        "one-guard product floor changed")

full_ordinary = Fraction(81, P**3 * (P - 1) ** 2)
full_guard = Fraction(63, P**3 * (P - 1) ** 2)
marginal_ordinary = P * full_ordinary
marginal_guard = P * full_guard
require(full_ordinary == Fraction(9, 35152), "ordinary full floor changed")
require(full_guard == Fraction(7, 35152), "guard full floor changed")
require(marginal_ordinary == Fraction(9, 2704), "ordinary marginal floor changed")
require(marginal_guard == Fraction(7, 2704), "guard marginal floor changed")

print("  ordinary paired safe capacity: 9")
print("  guard paired safe capacity: 7")
print("  at-most-one-guard product floor: 63 (ordinary/ordinary: 81)")
print(f"  full 13^3 coefficient floors: {full_guard} or {full_ordinary}")
print(f"  marginalized coefficient floors: {marginal_guard} or {marginal_ordinary}")


print("\n== double-zero finite-character corner ==")


def nonzero_character_sum(index):
    return P - 1 if index % P == 0 else -1


basis_checks = 0
for r in range(P):
    for s in range(P):
        for t in range(P):
            if r == 0 or s == 0 or r == t:
                continue
            require(
                nonzero_character_sum(r) * nonzero_character_sum(s) == 1,
                "double-zero basis cell has wrong corner weight",
            )
            basis_checks += 1
require(basis_checks == 1728, "allowed basis-cell count changed")

# A nontrivial exact array with both zero planes and the diagonal zero.
table_total = 0
corner_numerator = 0
for r in range(P):
    for s in range(P):
        for t in range(P):
            value = 0 if r == 0 or s == 0 or r == t else 1 + (3 * r + 5 * s + 7 * t) % 11
            table_total += value
            corner_numerator += (
                value * nonzero_character_sum(r) * nonzero_character_sum(s)
            )
require(table_total > 0, "control table vanished")
require(corner_numerator == table_total, "double-zero corner identity failed")

print(f"  admissible delta-basis cells: {basis_checks}")
print(f"  nontrivial control total/corner numerator: {table_total}")
print("  sum_(A,B!=0) Rhat(A,B,0) = total/13^3 exactly")


print("\n== missing-left-residue hostile ==")


def cyclotomic_reduce(coefficients):
    """Reduce a length-13 coefficient vector modulo Phi_13."""
    require(len(coefficients) == P, "cyclotomic vector has wrong length")
    top = coefficients[-1]
    return tuple(coefficients[j] - top for j in range(P - 1))


c = [Fraction(-1, P)] * P
c[0] = Fraction(P - 1, P)
for shift in range(P):
    coefficients = [Fraction(0)] * P
    for residue, weight in enumerate(c):
        coefficients[(-residue * shift) % P] += weight
    value = cyclotomic_reduce(coefficients)
    expected = (Fraction(int(shift != 0)),) + (Fraction(0),) * (P - 2)
    require(value == expected, "one-sided hostile reconstruction failed")
require(sum(c, Fraction(0)) == 0, "canonical equal-residue sum did not cancel")

print("  K(s)=1_(s!=0) has nonzero one-sided characters")
print("  assigning equal left/right residues makes every target difference zero")


print("\n== THM-2561 exact physical interval control ==")
J_LO = Fraction(4117, 399854)
J_HI = Fraction(4129, 399854)
X_LO = (J_LO + 3) / P
X_HI = (J_HI + 3) / P
C = P**5
A = P**2
K_A = 95
K_B = 93


def factor_walls(coefficient, shifts):
    result = set()
    for shift in shifts:
        lower_value = coefficient * X_LO + shift
        upper_value = coefficient * X_HI + shift
        for integer in range(floor_q(lower_value) - 2, floor_q(upper_value) + 4):
            for sign in (-1, 1):
                wall = (
                    Fraction(integer) + sign * Fraction(1, 14) - shift
                ) / coefficient
                if X_LO < wall < X_HI:
                    result.add(wall)
    return result


walls = {X_LO, X_HI}
walls |= factor_walls(C, [-Fraction(r, P) for r in range(P)])
walls |= factor_walls(A, [-Fraction(s, P) for s in range(P)])
walls |= factor_walls(K_A, [Fraction(s, P) for s in range(P)])
walls |= factor_walls(C, [-Fraction(t, P) for t in range(P)])
walls |= factor_walls(K_B, [Fraction(t, P) for t in range(P)])
walls = sorted(walls)

physical = defaultdict(Fraction)
minimum_counts = [P, P, P]
count_histogram = defaultdict(Fraction)
for left, right in zip(walls, walls[1:]):
    midpoint = (left + right) / 2
    deep_roots = [
        r for r in range(P) if danger(C * midpoint - Fraction(r, P))
    ]
    first_pair = [
        s for s in range(P)
        if not danger(A * midpoint - Fraction(s, P))
        and not danger(K_A * midpoint + Fraction(s, P))
    ]
    second_pair = [
        t for t in range(P)
        if not danger(C * midpoint - Fraction(t, P))
        and not danger(K_B * midpoint + Fraction(t, P))
    ]
    counts = (len(deep_roots), len(first_pair), len(second_pair))
    minimum_counts = [min(old, new) for old, new in zip(minimum_counts, counts)]
    count_histogram[counts] += right - left
    for r in deep_roots:
        for s in first_pair:
            for t in second_pair:
                physical[r, s, t] += right - left

rho = X_HI - X_LO
total = sum(physical.values(), Fraction(0))
require(rho == Fraction(6, 2599051), "physical packet mass changed")
require(len(walls) - 2 == 22, "physical wall count changed")
require(len(physical) == 1330, "physical positive-cell count changed")
require(minimum_counts == [1, 10, 10], "physical capacity minima changed")
require(total == Fraction(2110, 4826809), "physical table total changed")
require(total / rho == Fraction(7385, 39), "physical relative total changed")
require(max(physical.values()) == rho / 6, "physical maximum cell changed")
require(total >= 81 * rho, "physical table missed the ordinary capacity floor")
require(
    all(r != 0 and s != 0 and r != t for r, s, t in physical),
    "physical zero plane or diagonal became positive",
)
require(
    sum(count_histogram.values(), Fraction(0)) == rho,
    "physical chamber histogram lost mass",
)

print(f"  old-head mass: {rho}")
print(f"  strict walls / open chambers: {len(walls)-2} / {len(walls)-1}")
print(f"  positive table cells: {len(physical)}")
print(f"  minimum (deep,s-pair,t-pair) capacities: {tuple(minimum_counts)}")
print(f"  total table mass: {total} = ({total/rho}) rho")
print(f"  largest cell: {max(physical.values())} = (1/6) rho")
print("  r=0 plane, s=0 plane, and r=t diagonal are exactly zero")

print("\nall exact checks passed")
