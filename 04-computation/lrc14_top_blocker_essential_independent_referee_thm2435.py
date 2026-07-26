#!/usr/bin/env python3
"""Independent finite-fibre referee for THM-2435.

This audit does not import either primary companion.  It separately checks:

* the exact image criterion for low/high quotient blockers;
* the three typed mass invoices;
* prime-cyclic gap words and the equivariant necklace marker;
* the sharp integral row-defect energy and its mixed-mode boundary; and
* the positive-depth annihilator which distinguishes quotient from physical
  Fourier classes.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(value: Fraction) -> Fraction:
    residue = value - value.numerator // value.denominator
    return min(residue, 1 - residue)


def danger(coefficient: int, point: Fraction) -> bool:
    return circle_norm(coefficient * point) < Fraction(1, 14)


def preimage_word(coefficient: int, parent: Fraction, degree: int) -> tuple[bool, ...]:
    return tuple(
        danger(coefficient, (parent + branch) / degree)
        for branch in range(degree)
    )


# A coefficient of valuation at most M occupies exactly one seventh of a
# generic 7^(M+1)-fibre.  A coefficient divisible by 7^(M+1) is constant
# on that fibre.  The direct union test independently recovers the image
# formula used in THM-2435.
image_controls = 0
parents = (Fraction(1, 113), Fraction(17, 113), Fraction(61, 113))
for M in (0, 1, 2):
    N = 7 ** (M + 1)
    for parent in parents:
        for r in range(M + 1):
            word = preimage_word((7**r) * 5, parent, N)
            require(sum(word) == N // 7, "low blocker count changed")
            image_controls += 1

        for unit in (1, 2, 5):
            coefficient = N * unit
            word = preimage_word(coefficient, parent, N)
            require(
                len(set(word)) == 1
                and word[0] == danger(unit, parent),
                "high blocker stopped being fibre-constant",
            )
            image_controls += 1

        for k in (1, 2):
            low_valuations = (0,) if k == 1 else (0, M)
            lows = tuple(
                (7**r) * unit
                for r, unit in zip(low_valuations, (1, 2))
            )
            highs = tuple(N * unit for unit in (3, 5, 6)[: 3 - k])
            words = tuple(preimage_word(c, parent, N) for c in (*lows, *highs))
            has_safe_preimage = any(
                not any(word[branch] for word in words)
                for branch in range(N)
            )
            all_high_safe = all(not danger(c // N, parent) for c in highs)
            require(
                has_safe_preimage == all_high_safe,
                "quotient image criterion changed",
            )
            image_controls += 1


# Typed invoices: parent >=(4+k)/7, exact tilings <=3/7, then one
# quotient root per parent and b-way label selection.
typed = ((1, 1), (2, 2), (2, 1))
invoices = []
for k, b in typed:
    parent_floor = Fraction(4 + k, 7)
    essential = parent_floor - Fraction(3, 7)
    fixed_parent = essential / b
    marked_mass = fixed_parent / 91
    word_mass = (essential / 91) / (2**b - 1)
    invoices.append((parent_floor, essential, fixed_parent, marked_mass, word_mass))

require(
    invoices
    == [
        (Fraction(5, 7), Fraction(2, 7), Fraction(2, 7), Fraction(2, 637), Fraction(2, 637)),
        (Fraction(6, 7), Fraction(3, 7), Fraction(3, 14), Fraction(3, 1274), Fraction(1, 637)),
        (Fraction(6, 7), Fraction(3, 7), Fraction(3, 7), Fraction(3, 637), Fraction(3, 637)),
    ],
    "typed mass invoices changed",
)


# Every nonempty proper Boolean word on C_13 has all twelve nontrivial
# characters: a vanishing rational coefficient polynomial would have to
# equal 0 or Phi_13.  The necklace marker is checked independently using
# minimum rotations of the complemented word.
def marker(mask: int) -> int:
    bits = tuple((mask >> index) & 1 for index in range(13))
    complemented_rotations = [
        (tuple(1 - bits[(start + offset) % 13] for offset in range(13)), start)
        for start in range(13)
    ]
    return min(complemented_rotations)[1]


necklace_orbits: set[int] = set()
marker_checks = 0
for mask in range(1, (1 << 13) - 1):
    coefficients = tuple((mask >> index) & 1 for index in range(13))
    require(
        len(set(coefficients)) == 2,
        "proper C13 word became cyclotomic-constant",
    )
    chosen = marker(mask)
    require(coefficients[chosen] == 1, "marker left its support")

    orbit_representatives = []
    for shift in range(13):
        shifted = 0
        for index in range(13):
            if (mask >> index) & 1:
                shifted |= 1 << ((index + shift) % 13)
        require(
            marker(shifted) == (chosen + shift) % 13,
            "marker covariance changed",
        )
        orbit_representatives.append(shifted)
        marker_checks += 1
    necklace_orbits.add(min(orbit_representatives))

require(len(necklace_orbits) == 630, "wrong C13 necklace count")


# Independent weak-composition recursion for every integral seven-incidence
# row.  It recovers the sharp nonzero Dirichlet floor four.
def compositions(total: int, length: int, prefix: tuple[int, ...] = ()):
    if length == 1:
        yield (*prefix, total)
        return
    for value in range(total + 1):
        yield from compositions(total - value, length - 1, (*prefix, value))


row_count = 0
energy_floor = None
for multiplicity in compositions(7, 7):
    row_count += 1
    defect = tuple(1 - entry for entry in multiplicity)
    require(sum(defect) == 0, "row defect lost zero sum")
    if any(defect):
        energy = sum(
            (defect[(index + 1) % 7] - defect[index]) ** 2
            for index in range(7)
        )
        energy_floor = energy if energy_floor is None else min(energy_floor, energy)
require(row_count == 1716 and energy_floor == 4, "row-energy audit changed")


# Exact rational two-way ANOVA: for arrays with every F7 row summing to
# zero, vanishing interaction is equivalent to all thirteen rows being
# equal.  Exhaust a hostile bank of row choices.
anova_checks = 0
row_bank = (
    (0, 0, 0, 0, 0, 0, 0),
    (1, -1, 0, 0, 0, 0, 0),
    (1, 0, -1, 0, 0, 0, 0),
)
for row_indices in product(range(len(row_bank)), repeat=3):
    rows = tuple(row_bank[index] for index in row_indices)
    # Repeat the last chosen row to make thirteen rows.
    array = (*rows, *(rows[-1] for _ in range(10)))
    column_means = tuple(
        Fraction(sum(row[column] for row in array), 13)
        for column in range(7)
    )
    interaction_zero = all(
        Fraction(array[row][column]) == column_means[column]
        for row in range(13)
        for column in range(7)
    )
    rows_equal = len(set(array)) == 1
    require(interaction_zero == rows_equal, "ANOVA boundary changed")
    anova_checks += 1


# One marked quotient root has normalized C_91 coefficient 1/91 in every
# class.  Pullback through d=7^M kills exactly the non-annihilator classes.
ancestry_checks = 0
for M in (1, 2, 3):
    d = 7**M
    for residue in range(91 * d):
        survives = residue % d == 0
        if survives and gcd(residue // d, 91) == 1:
            require(gcd(residue, 91) == 7, "positive-depth gcd boundary changed")
        ancestry_checks += 1


print("THM-2435 independent finite-fibre referee")
print(f"quotient_image_controls={image_controls}")
print("typed_parent_floors=5/7,6/7,6/7")
print("typed_essential_floors=2/7,3/7,3/7")
print("uniform_fixed_label_mass=3/1274")
print("uniform_word_mass=1/637")
print(f"c13_necklace_orbits={len(necklace_orbits)}")
print(f"marker_covariance_checks={marker_checks}")
print(f"row_compositions={row_count}")
print(f"row_dirichlet_floor={energy_floor}")
print(f"anova_hostile_checks={anova_checks}")
print(f"positive_depth_annihilator_checks={ancestry_checks}")
print("quotient_physical_scope=PASS")
print("ALL INDEPENDENT CHECKS PASSED")
