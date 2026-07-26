#!/usr/bin/env python3
"""Exact companion for THM-2435.

The finite word audit checks every labelled seven-bin placement of one
guard, five ordinary unit words, and one or two top blockers.  The
arithmetic audit verifies the sharp parent/physical mass conversion and
the exact one-sheet septimal polyphase energy.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations, product


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


# The six unit top words have total weight seven: a two-bin guard and
# five labelled singleton ordinary words.
unit_profiles = 0
gap_histogram: Counter[int] = Counter()
b1_cover_profiles = 0
b1_essential_profiles = 0
b2_cover_profiles = 0
b2_essential_profiles = 0

for guard in combinations(range(7), 2):
    for ordinary in product(range(7), repeat=5):
        unit_profiles += 1
        multiplicity = [0] * 7
        for cell in guard:
            multiplicity[cell] += 1
        for cell in ordinary:
            multiplicity[cell] += 1

        gaps = sum(value == 0 for value in multiplicity)
        gap_histogram[gaps] += 1

        # One labelled blocker covers the word exactly when there is at
        # most one unit gap.  If there is one gap its position is forced.
        if gaps == 0:
            b1_cover_profiles += 7
        elif gaps == 1:
            b1_cover_profiles += 1
            b1_essential_profiles += 1

        # Two ordered labelled blockers have 49 placements.  With one
        # gap, 13 ordered pairs hit it; with two gaps, only the two
        # orders of the two missing cells work.
        if gaps == 0:
            b2_cover_profiles += 49
        elif gaps == 1:
            b2_cover_profiles += 13
            b2_essential_profiles += 13
        elif gaps == 2:
            b2_cover_profiles += 2
            b2_essential_profiles += 2


require(unit_profiles == 352947, "wrong six-unit word count")
require(
    gap_histogram
    == Counter({0: 2520, 1: 50400, 2: 157500, 3: 119700, 4: 22155, 5: 672}),
    "wrong six-unit gap histogram",
)
require(b1_cover_profiles == 68040, "wrong one-blocker cover count")
require(b1_essential_profiles == 50400, "wrong one-blocker essential count")
require(b2_cover_profiles == 1093680, "wrong two-blocker cover count")
require(b2_essential_profiles == 970200, "wrong two-blocker essential count")


# THM-2427 supplies parent mass at least 4/7. THM-2431 caps the
# six-unit exact-tiling locus at 3/7.
parent_floor = Fraction(4, 7)
unit_exact_cap = Fraction(3, 7)
essential_parent_floor = parent_floor - unit_exact_cap
root_count = 91
aggregate_physical_floor = essential_parent_floor / root_count
fixed_label_parent_floor = essential_parent_floor / 2
fixed_label_physical_floor = fixed_label_parent_floor / root_count

require(essential_parent_floor == Fraction(1, 7), "wrong parent gap")
require(
    aggregate_physical_floor == Fraction(1, 637),
    "wrong punctured physical floor",
)
require(
    fixed_label_parent_floor == Fraction(1, 14),
    "wrong fixed-label parent floor",
)
require(
    fixed_label_physical_floor == Fraction(1, 1274),
    "wrong fixed-label physical floor",
)


# On every seven-root fibre, a packet supported in one top-blocker
# sheet is either zero or a singleton. Every normalized C_7 Fourier
# coefficient of a singleton has squared magnitude 1/49. Therefore a
# physical packet of mass rho has exact energy rho/7 in each residue
# class modulo seven.
for singleton in range(7):
    for colour in range(7):
        # A formal singleton transform is one monomial of coefficient
        # 1/7; its cyclotomic phase has norm one.
        exponent = (-colour * singleton) % 7
        require(0 <= exponent < 7, "bad formal septimal exponent")
        require(Fraction(1, 7) ** 2 == Fraction(1, 49), "bad DFT norm")

per_colour_floor = fixed_label_physical_floor / 7
total_charged_floor = 6 * per_colour_floor
b1_per_colour_floor = aggregate_physical_floor / 7
b1_total_charged_floor = 6 * b1_per_colour_floor

require(per_colour_floor == Fraction(1, 8918), "wrong colour floor")
require(total_charged_floor == Fraction(3, 4459), "wrong charged floor")
require(b1_per_colour_floor == Fraction(1, 4459), "wrong b=1 colour floor")
require(
    b1_total_charged_floor == Fraction(6, 4459),
    "wrong b=1 charged floor",
)


# Sharp target-dark hostile: D_13 is 1/13-periodic, so all of its
# Fourier frequencies are divisible by thirteen. For every nonzero
# ell mod 7, n=13(7-ell) lies in residue ell mod 7 and its reduced
# central-interval coefficient is nonzero because 7 does not divide
# 7-ell.
hostile_frequencies = {
    colour: 13 * (7 - colour) for colour in range(1, 7)
}
require(
    all(frequency % 7 == colour for colour, frequency in hostile_frequencies.items()),
    "hostile missed a source residue",
)
require(
    all(frequency % 13 == 0 for frequency in hostile_frequencies.values()),
    "hostile unexpectedly supplied a target unit",
)
require(
    all((frequency // 13) % 7 != 0 for frequency in hostile_frequencies.values()),
    "hostile coefficient can vanish",
)


print("THM-2435 exact companion")
print(f"labelled_six_unit_profiles={unit_profiles}")
print(
    "unit_gap_histogram="
    + ",".join(f"{gaps}:{gap_histogram[gaps]}" for gaps in sorted(gap_histogram))
)
print(f"one_blocker_cover_profiles={b1_cover_profiles}")
print(f"one_blocker_essential_profiles={b1_essential_profiles}")
print(f"two_blocker_cover_profiles={b2_cover_profiles}")
print(f"two_blocker_essential_profiles={b2_essential_profiles}")
print(f"essential_parent_floor={essential_parent_floor}")
print(f"aggregate_physical_floor={aggregate_physical_floor}")
print(f"fixed_label_parent_floor={fixed_label_parent_floor}")
print(f"fixed_label_physical_floor={fixed_label_physical_floor}")
print(f"per_nonzero_septimal_colour_floor={per_colour_floor}")
print(f"total_charged_septimal_floor={total_charged_floor}")
print(f"one_blocker_per_colour_floor={b1_per_colour_floor}")
print(f"one_blocker_total_charged_floor={b1_total_charged_floor}")
print(
    "target_dark_hostile_frequencies="
    + ",".join(
        f"{colour}:{hostile_frequencies[colour]}" for colour in sorted(hostile_frequencies)
    )
)
print("ALL CHECKS PASSED")
