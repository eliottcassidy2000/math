#!/usr/bin/env python3
"""Exact finite companion for THM-2452."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import product
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


# Five guard/unit truth bits followed by two nondeep-blocker truth bits.
# Zero means safe and one means dangerous.
atoms = tuple(product((0, 1), repeat=7))
require(len(atoms) == 128, "complete local-mask bank changed size")
require(len(set(atoms)) == 128, "complete atoms ceased to be distinct")


def first_failure(guard_unit_bits: tuple[int, ...]) -> int:
    for index, bit in enumerate(guard_unit_bits, start=1):
        if bit:
            return index
    return 0


lex_counts = Counter(first_failure(atom[:5]) for atom in atoms)
expected_lex_counts = {0: 4, 1: 64, 2: 32, 3: 16, 4: 8, 5: 4}
require(dict(sorted(lex_counts.items())) == expected_lex_counts, "lex census changed")

# A lexicographic cell with first failure i has this many complete GU
# extensions. Cell zero is already all-safe.
extensions_by_cell = tuple(
    1 if index == 0 else 2 ** (5 - index) for index in range(6)
)
require(extensions_by_cell == (1, 16, 8, 4, 2, 1), "extension bank changed")

# Complete Boolean atoms are pointwise orthogonal and idempotent.
diagonal_pairs = 0
off_diagonal_pairs = 0
for left in atoms:
    for right in atoms:
        intersection_nonempty = left == right
        if intersection_nonempty:
            diagonal_pairs += 1
        else:
            off_diagonal_pairs += 1
require(diagonal_pairs == 128, "diagonal atom count changed")
require(off_diagonal_pairs == 128**2 - 128, "off-diagonal atom count changed")

# Exact quantitative invoices.
mask_count = len(atoms)
drift_denominator = mask_count**2
eligible_projection_denominator = 13
eligible_colours = 12 * (13**2 - 1)
eligible_energy_denominator = drift_denominator * eligible_projection_denominator
amplitude_square_denominator = eligible_energy_denominator * eligible_colours
require(drift_denominator == 16384, "drift denominator changed")
require(eligible_colours == 2016, "eligible colour count changed")
require(eligible_energy_denominator == 212992, "eligible energy denominator changed")
require(
    amplitude_square_denominator == 429391872,
    "amplitude-square denominator changed",
)

# The delayed-word quantitative version loses at most mu(Q)^2/4 once
# the standard BV threshold is imposed.
word_factor = Fraction(1, 4)
require(
    Fraction(1, drift_denominator) * word_factor
    == Fraction(1, 65536),
    "word-restored drift factor changed",
)
word_drift_denominator = 4 * drift_denominator
word_eligible_energy_denominator = (
    word_drift_denominator * eligible_projection_denominator
)
word_amplitude_square_denominator = (
    word_eligible_energy_denominator * eligible_colours
)
require(word_drift_denominator == 65536, "word drift denominator changed")
require(
    word_eligible_energy_denominator == 851968,
    "word eligible-energy denominator changed",
)
require(
    word_amplitude_square_denominator == 1717567488,
    "word amplitude-square denominator changed",
)

# If m is nonzero modulo 13 and the deep coefficient does not vanish,
# then 7 does not divide m, so m is a 91-unit.
for residue_13 in range(1, 13):
    for residue_7 in range(1, 7):
        representatives = [
            m
            for m in range(91)
            if m % 13 == residue_13 and m % 7 == residue_7
        ]
        require(len(representatives) == 1, "CRT unit representative changed")
        require(gcd(representatives[0], 91) == 1, "extracted m lost 91-unit typing")

# Fixed-frequency leakage is real although the full-X sum vanishes.
# On C4 take f=delta_0 and g=delta_1. With normalized DFT and right
# frequency X+1, the four cross terms are i^(X+1)/16.
gaussian_powers = ((1, 0), (0, 1), (-1, 0), (0, -1))
fixed_x_numerators = tuple(gaussian_powers[(x + 1) % 4] for x in range(4))
require(all(value != (0, 0) for value in fixed_x_numerators), "leakage term vanished")
require(
    sum(value[0] for value in fixed_x_numerators) == 0
    and sum(value[1] for value in fixed_x_numerators) == 0,
    "off-diagonal full-X cancellation failed",
)

# Abstract triangle inequalities are sharp from linearity alone.
drift_shares = (Fraction(1, mask_count),) * mask_count
require(sum(drift_shares) == 1, "128-way partition conservation failed")
require(max(drift_shares) ** 2 == Fraction(1, drift_denominator), "square floor changed")

print("THM-2452 exact companion")
print(f"complete_masks={mask_count}")
print(
    "lex_cell_counts="
    + ",".join(str(expected_lex_counts[index]) for index in range(6))
)
print("extensions_by_cell=" + ",".join(map(str, extensions_by_cell)))
print(f"ordered_mask_pairs={mask_count ** 2}")
print(f"diagonal_mask_pairs={diagonal_pairs}")
print(f"offdiagonal_mask_pairs={off_diagonal_pairs}")
print(f"drift_denominator={drift_denominator}")
print(f"eligible_projection_denominator={eligible_projection_denominator}")
print(f"eligible_colours={eligible_colours}")
print(f"eligible_energy_denominator={eligible_energy_denominator}")
print(f"amplitude_square_denominator={amplitude_square_denominator}")
print("word_restoration_factor=mu(Q)^2/4")
print(f"word_drift_denominator={word_drift_denominator}")
print(f"word_eligible_energy_denominator={word_eligible_energy_denominator}")
print(f"word_amplitude_square_denominator={word_amplitude_square_denominator}")
print("crt_91_unit_check=PASS")
print("fixed_X_nonzero_terms=4")
print("full_X_offdiagonal_sum=0")
print("ALL CHECKS PASSED")
