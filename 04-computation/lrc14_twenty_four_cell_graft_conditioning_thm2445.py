#!/usr/bin/env python3
"""Exact algebraic companion for THM-2445."""

from __future__ import annotations

from fractions import Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def first_failure_cell(bits: tuple[int, ...]) -> int:
    """0 means all safe; i=1..5 means the first danger occurs at i."""

    for index, bit in enumerate(bits, start=1):
        if bit == 0:
            return index
    return 0


# 1. Six lexicographic cells and four blocker states.
lex_counts = [0] * 6
partition_checks = 0
for mask in range(1 << 5):
    bits = tuple((mask >> index) & 1 for index in range(5))
    cell = first_failure_cell(bits)
    indicators = []
    indicators.append(int(all(bits)))
    for index in range(5):
        indicators.append(int(all(bits[:index]) and bits[index] == 0))
    require(sum(indicators) == 1, "lexicographic cells did not partition")
    require(indicators[cell] == 1, "first-failure label changed")
    lex_counts[cell] += 1
    partition_checks += 1

require(lex_counts == [1, 16, 8, 4, 2, 1], "lexicographic census changed")

product_cells = tuple(
    (cell, source_danger, target_danger)
    for cell in range(6)
    for source_danger in range(2)
    for target_danger in range(2)
)
ghosts = tuple(
    state
    for state in product_cells
    if state == (0, 0, 0)
)
typed_cells = tuple(state for state in product_cells if state not in ghosts)
require(len(product_cells) == 24, "wrong product-cell count")
require(len(typed_cells) == 23, "wrong typed-cell count")
require(ghosts == ((0, 0, 0),), "ghost ceased to be unique")


# 2. Exact quantitative ledger.
cell_count = 24
cell_energy_denominator = cell_count**2
eligible_energy_denominator = 13 * cell_energy_denominator
eligible_colours = 12 * (13**2 - 1)
coefficient_square_denominator = (
    eligible_energy_denominator * eligible_colours
)
source_variance_denominator = 294
source_energy_denominator = (
    coefficient_square_denominator * source_variance_denominator
)

require(cell_energy_denominator == 576, "cell energy denominator changed")
require(eligible_energy_denominator == 7488, "eligible energy denominator changed")
require(eligible_colours == 2016, "eligible colour count changed")
require(coefficient_square_denominator == 15095808, "coefficient floor changed")
require(source_energy_denominator == 4438167552, "source energy floor changed")
require(13 * eligible_colours == 26208, "amplitude radical changed")


# 3. Exact rational rank of the C7 cyclotomic reduction.
def matrix_rank(matrix: list[list[Fraction]]) -> int:
    rows = [row[:] for row in matrix]
    if not rows:
        return 0
    row = 0
    columns = len(rows[0])
    for column in range(columns):
        pivot = next(
            (candidate for candidate in range(row, len(rows)) if rows[candidate][column]),
            None,
        )
        if pivot is None:
            continue
        rows[row], rows[pivot] = rows[pivot], rows[row]
        scale = rows[row][column]
        rows[row] = [entry / scale for entry in rows[row]]
        for candidate in range(len(rows)):
            if candidate == row or not rows[candidate][column]:
                continue
            factor = rows[candidate][column]
            rows[candidate] = [
                left - factor * right
                for left, right in zip(rows[candidate], rows[row])
            ]
        row += 1
    return row


cyclotomic_ranks = []
for character in range(1, 7):
    # Coefficients after X^7=1 are permuted by character.  Reduction
    # modulo Phi_7 subtracts the X^6 coefficient from degrees 0..5.
    reduction = [[Fraction(0) for _ in range(7)] for _ in range(6)]
    permuted_positions = [(character * exponent) % 7 for exponent in range(7)]
    for original, position in enumerate(permuted_positions):
        if position < 6:
            reduction[position][original] += 1
        else:
            for degree in range(6):
                reduction[degree][original] -= 1
    rank = matrix_rank(reduction)
    require(rank == 6, "cyclotomic reduction lost rank")
    for row in reduction:
        require(sum(row) == 0, "constant vector left cyclotomic kernel")
    cyclotomic_ranks.append(rank)


# 4. Sharp anchored-zero energy: z_0=0, z_1=...=z_6=1/6.
z = (Fraction(0),) + (Fraction(1, 6),) * 6
total_source_energy = Fraction(1, 7) * sum(value * value for value in z)
trivial_source_energy = (sum(z) / 7) ** 2
nontrivial_source_energy = total_source_energy - trivial_source_energy
require(nontrivial_source_energy == Fraction(1, 294), "sharp source energy changed")


print("THM-2445 exact companion")
print(f"five_bit_partition_checks={partition_checks}")
print("lex_cell_counts=" + ",".join(map(str, lex_counts)))
print(f"product_cells={len(product_cells)}")
print(f"typed_cells={len(typed_cells)}")
print(f"ghost_cells={len(ghosts)}")
print(f"cell_energy_denominator={cell_energy_denominator}")
print(f"eligible_energy_denominator={eligible_energy_denominator}")
print(f"eligible_colours={eligible_colours}")
print(f"coefficient_square_denominator={coefficient_square_denominator}")
print(f"source_variance_denominator={source_variance_denominator}")
print(f"source_energy_denominator={source_energy_denominator}")
print("cyclotomic_reduction_ranks=" + ",".join(map(str, cyclotomic_ranks)))
print(f"sharp_nontrivial_source_energy={nontrivial_source_energy}")
print("ALL CHECKS PASSED")
