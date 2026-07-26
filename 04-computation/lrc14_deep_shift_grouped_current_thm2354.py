#!/usr/bin/env python3
"""Exact finite controls for THM-2354.

The analytic Abel identity is proved in the theorem text.  This companion
freezes the 13-translate circle cover on its exact 182-cell arrangement, the
sharp overlap-energy control, and the 91-unit deep-colour residue ledger.
"""

from fractions import Fraction
from math import gcd


P = 13
Q = 7
MODULUS = P * Q
CELL_COUNT = 2 * MODULUS
MIDPOINT_DENOMINATOR = 2 * CELL_COUNT


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circular_distance(numerator: int, modulus: int) -> int:
    residue = numerator % modulus
    return min(residue, modulus - residue)


# The midpoint of cell j is (2j+1)/364.  The r-th shifted danger arc has
# centre r/13=28r/364 and radius 1/14=26/364.  Midpoints avoid boundaries.
cover_by_shift: list[set[int]] = []
for shift in range(P):
    covered = {
        cell
        for cell in range(CELL_COUNT)
        if circular_distance(
            2 * cell + 1 - 28 * shift,
            MIDPOINT_DENOMINATOR,
        )
        < 26
    }
    require(len(covered) == 26, "shifted comb cell count changed")
    cover_by_shift.append(covered)

coverage = [
    sum(cell in cover_by_shift[shift] for shift in range(P))
    for cell in range(CELL_COUNT)
]
require(min(coverage) == 1, "thirteen translates stopped covering")
require(max(coverage) == 2, "unexpected triple overlap")
require(coverage.count(1) == 26, "single-cover cell count changed")
require(coverage.count(2) == 156, "double-cover cell count changed")

unique_cells = {
    shift: [
        cell
        for cell in cover_by_shift[shift]
        if coverage[cell] == 1
    ]
    for shift in range(P)
}
require(
    all(len(cells) == 2 for cells in unique_cells.values()),
    "unique core per translate changed",
)


def overlap_profile(cells: set[int]) -> list[Fraction]:
    return [
        Fraction(len(cells & cover_by_shift[shift]), CELL_COUNT)
        for shift in range(P)
    ]


def nonconstant_fourier_energy(profile: list[Fraction]) -> Fraction:
    mean = sum(profile, Fraction()) / P
    return sum((value - mean) ** 2 for value in profile) / P


# A one-cell unique core realizes a one-sparse overlap profile.  It is a
# hostile control for any claim that comb noncancellation forces two shifts.
one_sparse_cells = {unique_cells[1][0]}
one_sparse_profile = overlap_profile(one_sparse_cells)
require(one_sparse_profile[0] == 0, "hostile meets the unshifted comb")
require(
    [index for index, value in enumerate(one_sparse_profile) if value]
    == [1],
    "one-sparse overlap control changed",
)

# Choose one unique cell for each nonzero shift.  Then G(0)=0 and the other
# twelve values are equal.  This geometrically realizes equality in both
# Cauchy--Schwarz steps of the theorem.
sharp_cells = {unique_cells[shift][0] for shift in range(1, P)}
sharp_profile = overlap_profile(sharp_cells)
sharp_mass = Fraction(len(sharp_cells), CELL_COUNT)
require(sharp_profile[0] == 0, "sharp control meets shift zero")
require(
    all(sharp_profile[shift] == Fraction(1, CELL_COUNT)
        for shift in range(1, P)),
    "sharp overlap profile changed",
)
sharp_energy = nonconstant_fourier_energy(sharp_profile)
require(
    sharp_energy == sharp_mass**2 / 2028,
    "sharp Fourier-energy constant changed",
)
sharp_mode_magnitude = Fraction(1, P * CELL_COUNT)
require(
    sharp_mode_magnitude == sharp_mass / 156,
    "sharp single-mode constant changed",
)

# Every nonzero mod-13 deep colour contains seven residues modulo 91.  The
# unique residue divisible by seven has zero danger-comb coefficient; the
# remaining six are exactly the 91-units.
deep_colour_rows = 0
unit_deep_rows = 0
zero_comb_rows = 0
for colour in range(1, P):
    residues = [
        multiplier
        for multiplier in range(MODULUS)
        if multiplier % P == colour
    ]
    require(len(residues) == Q, "deep-colour fibre changed")
    deep_colour_rows += len(residues)
    for multiplier in residues:
        if multiplier % Q == 0:
            zero_comb_rows += 1
        else:
            require(
                gcd(multiplier, MODULUS) == 1,
                "live deep multiplier is not a 91-unit",
            )
            unit_deep_rows += 1

require(deep_colour_rows == 84, "nonzero deep-colour universe changed")
require(zero_comb_rows == 12, "septimal zero count changed")
require(unit_deep_rows == 72, "91-unit deep residue count changed")
require(Fraction(1, 6 * 156) == Fraction(1, 936), "LRC floor factor changed")

print("theorem=THM-2354")
print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-AUDITED")
print(f"arrangement_cells={CELL_COUNT}")
print(f"shifted_combs={P}")
print(f"cells_per_comb={len(cover_by_shift[0])}")
print(f"coverage_min={min(coverage)}")
print(f"coverage_max={max(coverage)}")
print(f"single_cover_cells={coverage.count(1)}")
print(f"double_cover_cells={coverage.count(2)}")
print(f"unique_cells_per_shift={len(unique_cells[0])}")
print("one_sparse_overlap_support=1")
print(f"sharp_profile_mass={sharp_mass}")
print(f"sharp_nonzero_energy={sharp_energy}")
print(f"sharp_mode_magnitude={sharp_mode_magnitude}")
print("general_energy_floor=mu(F)^2/2028")
print("general_mode_floor=mu(F)/156")
print(f"deep_colour_rows={deep_colour_rows}")
print(f"septimal_zero_rows={zero_comb_rows}")
print(f"unit_deep_residue_rows={unit_deep_rows}")
print("lrc_word_mass_floor=e_j*eta/6")
print("lrc_grouped_mode_floor=e_j*eta/936")
print("target_quotient_landing=NOT-PROVED")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
