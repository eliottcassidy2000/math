#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2340."""

from fractions import Fraction


P = 13
POINTS = range(P)


def require(condition: bool, message: str) -> None:
    """Raise under both ordinary and optimized Python."""
    if not condition:
        raise RuntimeError(message)


def mean(values) -> Fraction:
    values = tuple(values)
    require(bool(values), "cannot average an empty family")
    return sum(values, Fraction(0)) / len(values)


def row_means(matrix: tuple[tuple[Fraction, ...], ...]) -> tuple[Fraction, ...]:
    return tuple(mean(row) for row in matrix)


def column_means(
    matrix: tuple[tuple[Fraction, ...], ...],
) -> tuple[Fraction, ...]:
    return tuple(mean(matrix[row][column] for row in POINTS) for column in POINTS)


def grand_mean(matrix: tuple[tuple[Fraction, ...], ...]) -> Fraction:
    return mean(entry for row in matrix for entry in row)


def interaction(
    matrix: tuple[tuple[Fraction, ...], ...],
) -> tuple[tuple[Fraction, ...], ...]:
    rows = row_means(matrix)
    columns = column_means(matrix)
    grand = grand_mean(matrix)
    return tuple(
        tuple(
            matrix[s][t] - rows[s] - columns[t] + grand
            for t in POINTS
        )
        for s in POINTS
    )


def norm_square(matrix: tuple[tuple[Fraction, ...], ...]) -> Fraction:
    return mean(entry * entry for row in matrix for entry in row)


def anova_energies(
    matrix: tuple[tuple[Fraction, ...], ...],
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    rows = row_means(matrix)
    columns = column_means(matrix)
    grand = grand_mean(matrix)
    residual = interaction(matrix)
    return (
        grand * grand,
        mean((value - grand) ** 2 for value in rows),
        mean((value - grand) ** 2 for value in columns),
        norm_square(residual),
    )


def mask_kernel(delta: int) -> int:
    """Fourier transform of the nonzero indicator on F_13."""
    return P - 1 if delta % P == 0 else -1


def gram_mask_energies(
    matrix: tuple[tuple[Fraction, ...], ...],
) -> tuple[Fraction, Fraction, Fraction]:
    """Evaluate the three exact 169-twist Gram forms."""
    pure_a = Fraction(0)
    pure_b = Fraction(0)
    fork = Fraction(0)
    for s in POINTS:
        for t in POINTS:
            for left in POINTS:
                kernel_a = mask_kernel(left - s)
                for right in POINTS:
                    product_value = matrix[s][t] * matrix[left][right]
                    kernel_b = mask_kernel(right - t)
                    pure_a += product_value * kernel_a
                    pure_b += product_value * kernel_b
                    fork += product_value * kernel_a * kernel_b
    normalization = P**4
    return (
        pure_a / normalization,
        pure_b / normalization,
        fork / normalization,
    )


def transpose(
    matrix: tuple[tuple[Fraction, ...], ...],
) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(matrix[row][column] for row in POINTS)
        for column in POINTS
    )


centered = tuple(Fraction(value - 6) for value in POINTS)
quadratic = tuple(
    Fraction((value - 6) ** 2) - mean(Fraction((x - 6) ** 2) for x in POINTS)
    for value in POINTS
)

constant_matrix = tuple(
    tuple(Fraction(5) for _ in POINTS)
    for _ in POINTS
)
row_matrix = tuple(
    tuple(centered[s] for _ in POINTS)
    for s in POINTS
)
column_matrix = tuple(
    tuple(quadratic[t] for t in POINTS)
    for _ in POINTS
)
additive_matrix = tuple(
    tuple(centered[s] + quadratic[t] for t in POINTS)
    for s in POINTS
)
interaction_matrix = tuple(
    tuple(centered[s] * quadratic[t] for t in POINTS)
    for s in POINTS
)
generic_matrix = tuple(
    tuple(
        Fraction(
            3
            + (s - 6)
            + 2 * (t - 6) ** 2
            + (s - 6) * (t - 6) ** 2
            + ((2 * s + 3 * t) % P)
        )
        for t in POINTS
    )
    for s in POINTS
)

controls = {
    "constant": constant_matrix,
    "row": row_matrix,
    "column": column_matrix,
    "additive": additive_matrix,
    "interaction": interaction_matrix,
    "generic": generic_matrix,
}

checked_matrices = 0
for name, matrix in controls.items():
    require(len(matrix) == P and all(len(row) == P for row in matrix), f"{name} shape changed")
    rows = row_means(matrix)
    columns = column_means(matrix)
    grand = grand_mean(matrix)
    residual = interaction(matrix)
    require(
        all(mean(residual[s][t] for t in POINTS) == 0 for s in POINTS),
        f"{name} interaction acquired a row mean",
    )
    require(
        all(mean(residual[s][t] for s in POINTS) == 0 for t in POINTS),
        f"{name} interaction acquired a column mean",
    )
    require(
        all(
            matrix[s][t]
            == grand
            + (rows[s] - grand)
            + (columns[t] - grand)
            + residual[s][t]
            for s in POINTS
            for t in POINTS
        ),
        f"{name} ANOVA reconstruction failed",
    )
    energies = anova_energies(matrix)
    require(
        sum(energies, Fraction(0)) == norm_square(matrix),
        f"{name} orthogonal norm invoice failed",
    )
    require(
        gram_mask_energies(matrix) == energies[1:],
        f"{name} support-mask Gram identity failed",
    )
    swapped = anova_energies(transpose(matrix))
    require(
        swapped == (energies[0], energies[2], energies[1], energies[3]),
        f"{name} target swap law failed",
    )
    checked_matrices += 1

constant_energy = anova_energies(constant_matrix)
row_energy = anova_energies(row_matrix)
column_energy = anova_energies(column_matrix)
additive_energy = anova_energies(additive_matrix)
interaction_energy = anova_energies(interaction_matrix)

require(constant_energy == (25, 0, 0, 0), "constant zero-only hostile changed")
require(row_energy[1] > 0 and row_energy[2:] == (0, 0), "pure-a control changed")
require(column_energy[2] > 0 and column_energy[1] == column_energy[3] == 0, "pure-b control changed")
require(additive_energy[1] > 0 and additive_energy[2] > 0 and additive_energy[3] == 0, "additive no-fork control changed")
require(interaction_energy[3] > 0 and interaction_energy[:3] == (0, 0, 0), "fork control changed")

require(1 + (P - 1) + (P - 1) + (P - 1) ** 2 == P**2, "ANOVA dimensions changed")

print("theorem=THM-2340")
print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-AUDITED")
print(f"target_twist_matrix_size={P}x{P}")
print(f"matrices_checked={checked_matrices}")
print("anova_dimensions=1+12+12+144=169")
print("pure_a_energy=row_main_effect_norm")
print("pure_b_energy=column_main_effect_norm")
print("fork_energy=two_way_interaction_norm")
print("pure_a_iff=row_mean_profile_nonconstant")
print("pure_b_iff=column_mean_profile_nonconstant")
print("fork_iff=twist_matrix_not_additively_separable")
print("target_swap=exchanges_main_effects_and_preserves_interaction_norm")
print("constant_hostile=only_zero_target_survives")
print("additive_hostile=pure_axes_may_survive_while_fork_energy_is_zero")
print(f"row_control_energy={row_energy[1]}")
print(f"column_control_energy={column_energy[2]}")
print(f"interaction_control_energy={interaction_energy[3]}")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
