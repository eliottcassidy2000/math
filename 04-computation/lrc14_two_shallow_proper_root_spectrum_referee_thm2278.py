#!/usr/bin/env python3
"""Independent exact referee for THM-2278.

This path checks all complement-reduced masks by cyclotomic resultants and
checks the unit-comb root counts by an exact endpoint-cell construction.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import ceil

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(x: Fraction) -> Fraction:
    residue = x % 1
    return min(residue, 1 - residue)


# Exact cyclotomic norms of all 4,095 masks of sizes one through six.
x = sp.symbols("x")
phi_13 = sum(x**j for j in range(13))
norm_histogram: Counter[int] = Counter()
mask_count = 0

for size in range(1, 7):
    for support in combinations(range(13), size):
        polynomial = sum(x**j for j in support)
        norm = abs(int(sp.resultant(phi_13, polynomial, x)))
        require(norm >= 1, f"zero cyclotomic norm on {support}")
        norm_histogram[norm] += 1
        mask_count += 1

require(mask_count == 4095, "reduced-mask count")
require(
    sorted(norm_histogram.items())
    == [
        (1, 754),
        (27, 364),
        (53, 936),
        (79, 624),
        (131, 312),
        (157, 312),
        (313, 156),
        (547, 156),
        (599, 156),
        (625, 39),
        (729, 130),
        (911, 156),
    ],
    "cyclotomic norm histogram",
)

# Independently build every endpoint-free parent-phase cell for
# D_u on a thirteen-root fibre, for all nonzero residues u mod 13.
root_cell_count = 0
root_count_histogram: Counter[int] = Counter()

for unit in range(1, 13):
    endpoints = {Fraction(0), Fraction(1)}
    for root in range(13):
        for integer_part in range(-2, unit + 3):
            for sign in (-1, 1):
                boundary = (
                    Fraction(13, unit)
                    * (integer_part + sign * Fraction(1, 14))
                    - root
                )
                if 0 <= boundary <= 1:
                    endpoints.add(boundary)

    ordered = sorted(endpoints)
    local_counts = set()
    for left, right in zip(ordered, ordered[1:]):
        if left == right:
            continue
        parent = (left + right) / 2
        count = sum(
            circle_norm(Fraction(unit, 13) * (parent + root))
            < Fraction(1, 14)
            for root in range(13)
        )
        require(
            count in (1, 2),
            f"bad root count unit={unit} parent={parent}",
        )
        local_counts.add(count)
        root_count_histogram[count] += 1
        root_cell_count += 1
    require(
        local_counts == {1, 2},
        f"incomplete root-count spectrum for unit {unit}",
    )

require(
    root_count_histogram == Counter({1: 90, 2: 78}),
    "unit root-count histogram",
)
require(root_cell_count == 168, "unit root-cell count")

# Referee arithmetic uses the already proved THM-2273 image floor directly.
mode_floor = Fraction(5, 21) ** 5
image_floor = Fraction(5696989, 76962600)
integral_floor = mode_floor * image_floor
residue_floor = integral_floor / 13**2

require(
    mode_floor == Fraction(3125, 4084101),
    "cyclotomic floor fraction",
)
require(
    mode_floor * Fraction(21, 5) ** 5 == 1,
    "norm/AM-GM certificate",
)
require(
    integral_floor == Fraction(712123625, 12572921264904),
    "vector-mode integral floor",
)
require(
    residue_floor == Fraction(712123625, 2124823693768776),
    "residue energy floor",
)

gap_counts = {
    offset: ceil(Fraction(7 * 13**offset, 6) * image_floor)
    for offset in range(4)
}
require(
    gap_counts == {0: 1, 1: 2, 2: 15, 3: 190},
    "minimum gap counts",
)

print("THM-2278 independent exact referee")
print(f"cyclotomic_masks_checked={mask_count}")
print(f"cyclotomic_norm_histogram={dict(sorted(norm_histogram.items()))}")
print(f"unit_root_cells_checked={root_cell_count}")
print(f"unit_root_count_histogram={dict(sorted(root_count_histogram.items()))}")
print(f"cyclotomic_mode_floor={mode_floor}")
print(f"common_time_image_floor={image_floor}")
print(f"vector_mode_integral_floor={integral_floor}")
print(f"every_residue_energy_floor={residue_floor}")
print(f"minimum_gap_counts_by_depth_offset={gap_counts}")
print("exact_frequency_landing=NOT_PROVED")
print("ALL_CHECKS_PASSED")
